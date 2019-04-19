// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include "ie-solver/ie_mat.h"
#include "ie-solver/initialization.h"
#include "ie-solver/tools/ie_solver_tools.h"
#include "ie-solver/quadtree.h"
#include "ie-solver/kernel.h"
#include "ie-solver/helpers.h"
#include "ie-solver/log.h"
#include "ie-solver/ie_solver_config.h"
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/boundaries/circle.h"
#include "ie-solver/boundaries/rounded_square.h"
#include "ie-solver/boundaries/rounded_square_with_bump.h"
#include "ie-solver/boundaries/squiggly.h"
#include "ie-solver/boundaries/annulus.h"
#include "ie-solver/boundaries/cubic_spline.h"


namespace ie_solver {

void initialize_U_mat(const std::vector<Hole>& holes,
                      const std::vector<double>& tgt_points, ie_Mat* U) {
  for (unsigned int i = 0; i < tgt_points.size(); i += 2) {
    for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
      Hole hole = holes[hole_idx];

      Vec2 r = Vec2(tgt_points[i], tgt_points[i + 1]) - hole.center;
      double scale = 1.0 / (4 * M_PI);
      // TODO(John) implement tensor product to improve readability?
      U->set(i, 3 * hole_idx, scale *
             (log(1 / r.norm()) +
              (1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[0]));
      U->set(i + 1, 3 * hole_idx, scale *
             (log(1 / r.norm()) +
              (1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[0]));
      U->set(i, 3 * hole_idx + 1, scale *
             (log(1 / r.norm()) +
              (1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[1]));
      U->set(i + 1, 3 * hole_idx + 1, scale *
             (log(1 / r.norm()) +
              (1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[1]));

      U->set(i, 3 * hole_idx + 2, r.a[1] * (scale / pow(r.norm(), 2)));
      U->set(i + 1, 3 * hole_idx + 2, -r.a[0] * (scale / pow(r.norm(), 2)));
    }
  }
}


void initialize_Psi_mat(const std::vector<Hole>& holes,
                        Boundary* boundary, ie_Mat* Psi) {
  for (unsigned int i = 0; i < boundary->points.size(); i += 2) {
    Vec2 x = Vec2(boundary->points[i], boundary->points[i + 1]);
    for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
      Hole hole = holes[hole_idx];
      if ((x - hole.center).norm() < hole.radius + 1e-8) {
        Psi->set(3 * hole_idx, i, boundary->weights[i / 2]);
        Psi->set(3 * hole_idx + 1, i + 1, boundary->weights[i / 2]);
        Psi->set(3 * hole_idx + 2, i, boundary->weights[i / 2]*x.a[1]);
        Psi->set(3 * hole_idx + 2, i + 1, -boundary->weights[i / 2]*x.a[0]);
        break;
      }
    }
  }
}


// TODO clean up this mess
void schur_solve(const ie_solver_config & config,
                 const Kernel& kernel, const IeSolverTools& ie_solver_tools,
                 const QuadTree& quadtree, const std::vector<double>& domain_points,
                 ie_Mat* solution) {
  int domain_dimension = 2;
  int solution_dimension = 1;
  if (config.pde == ie_solver_config::STOKES) {
    solution_dimension = 2;
  }
  ie_Mat f = config.boundary->boundary_values;
  ie_Mat mu(config.num_boundary_points * solution_dimension, 1);

  std::vector<Hole> holes = config.boundary->holes;
  int num_holes = holes.size();
  ie_Mat alpha(3 * num_holes, 1);
  ie_Mat U(config.num_boundary_points * solution_dimension, 3 * num_holes);
  ie_Mat Psi(3 * num_holes, config.num_boundary_points * solution_dimension);

  ie_Mat ident;
  ident.eye(3 * num_holes);

  initialize_U_mat(holes, config.boundary->points, &U);
  initialize_Psi_mat(holes, config.boundary.get(), &Psi);

  ie_Mat Dinv_U(config.num_boundary_points * solution_dimension, 3 * num_holes);

  ie_solver_tools.solve(kernel, quadtree, &Dinv_U, U);

  ie_Mat Psi_Dinv_U(3 * num_holes, 3 * num_holes);
  ie_Mat::gemm(NORMAL, NORMAL, 1., Psi, Dinv_U, 0., &Psi_Dinv_U);

  ie_Mat S = ident;
  S += Psi_Dinv_U;
  S *= -1.;

  ie_Mat Dinv_u(config.num_boundary_points * solution_dimension, 1);
  ie_solver_tools.solve(kernel, quadtree, &Dinv_u, f);

  ie_Mat Psi_Dinv_u(3 * num_holes, 1);
  ie_Mat::gemv(NORMAL, 1., Psi, Dinv_u, 0., &Psi_Dinv_u);

  ie_Mat Sinv_Psi_Dinv_u(3 * num_holes, 1);

  S.left_multiply_inverse(Psi_Dinv_u, &Sinv_Psi_Dinv_u);
  alpha = Sinv_Psi_Dinv_u;
  alpha *= -1.;

  ie_Mat right_vec(config.num_boundary_points * solution_dimension, 1);
  ie_Mat::gemv(NORMAL, 1., U, Sinv_Psi_Dinv_u, 0., &right_vec);

  right_vec += f;
  ie_solver_tools.solve(kernel, quadtree, &mu, right_vec);

  // This will be done as a sparse mat vec in the future, for now we do
  // dense matvec

  // QuadTree boundary_to_domain;
  // boundary_to_domain.initialize_tree(config.boundary.get(),
  //                                     domain_points, solution_dimension,
  //                                     domain_dimension);
  // ie_solver_tools.b2dskeletonize(kernel, &boundary_to_domain);
  // ie_solver_tools.b2dsparse_matvec(kernel, boundary_to_domain, mu, &domain);

  ie_Mat K_domain(TEST_SIZE * TEST_SIZE * solution_dimension,
                  config.num_boundary_points * solution_dimension);
  Initialization init;
  init.InitializeDomainKernel(&K_domain,
                              domain_points, TEST_SIZE, kernel,
                              solution_dimension);

  ie_Mat::gemv(NORMAL, 1., K_domain, mu, 0., solution);

  ie_Mat U_forward(solution->height(), 3 * num_holes);

  initialize_U_mat(holes, domain_points, &U_forward);

  for (unsigned int i = 0; i < solution->height(); i += 2) {
    Vec2 point = Vec2(domain_points[i], domain_points[i + 1]);
    if (!config.boundary->is_in_domain(point)) {
      for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
        U_forward.set(i, 3 * hole_idx, 0);
        U_forward.set(i + 1, 3 * hole_idx, 0);
        U_forward.set(i, 3 * hole_idx + 1, 0);
        U_forward.set(i + 1, 3 * hole_idx + 1, 0);
        U_forward.set(i, 3 * hole_idx + 2, 0);
        U_forward.set(i + 1, 3 * hole_idx + 2, 0);
      }
    }
  }
  ie_Mat U_forward_alpha(solution->height(), 1);
  ie_Mat::gemv(NORMAL, 1., U_forward, alpha, 0., &U_forward_alpha);

  (*solution) += U_forward_alpha;
}

// TODO so config should be constant, currently the animation is hacked together,
// some mild refactoring should be done to handle the perturbation of the boundary
// which belongs to this constant class. The likely scenario is going to be to
// remove the boundary from the config altogether. THIS IS URGENT.
double boundary_integral_solve(ie_solver_config & config,
                               std::vector<double>* skel_times) {
  int domain_dimension = 2;
  int solution_dimension = 1;
  if (config.pde == ie_solver_config::STOKES) {
    solution_dimension = 2;
  }

  bool is_time_trial = (skel_times != nullptr);
  double id_tol = config.id_tol;

  if (!is_time_trial && !config.testing) {
    write_boundary_to_file(config.boundary->points);
  }


  QuadTree quadtree;
  quadtree.initialize_tree(config.boundary.get(), std::vector<double>(),
                           solution_dimension,
                           domain_dimension);

  if (!is_time_trial && !config.testing) {
    quadtree.write_quadtree_to_file();
  }

  // Consider making init instead of constructor for readability
  bool strong_admissibility =
    (config.admissibility == ie_solver_config::STRONG);
  IeSolverTools ie_solver_tools(id_tol, strong_admissibility,
                                solution_dimension, domain_dimension);

  std::vector<double> domain_points;
  get_domain_points(&domain_points, quadtree.min, quadtree.max);

  Kernel kernel;
  kernel.load(config.boundary.get(), domain_points, config.pde,
              solution_dimension,
              domain_dimension);

  if (is_time_trial) {
    double elapsed = 0;
    for (int i = 0; i < TIMING_ITERATIONS; i++) {
      double start = omp_get_wtime();
      ie_solver_tools.skeletonize(kernel, &quadtree);
      double end = omp_get_wtime();
      elapsed += (end - start);
      quadtree.reset();
    }
    skel_times->push_back(elapsed / TIMING_ITERATIONS);
    return -1;
  }


  for (int frame = 0; frame < 60; frame++) {

    double x = 0.35 + 0.3 * (frame / 60);
    double y = 0.35 ;

    config.boundary->holes[0].center = Vec2(x, y);


    double x1 = 0.65 - 0.3 * (frame / 60);
    double y1 = 0.65;

    config.boundary->holes[1].center = Vec2(x1, y1);



    config.boundary->initialize(1000,
                                config.boundary_condition);
    config.num_boundary_points = config.boundary->weights.size();
    quadtree.reset();
    ie_solver_tools.skeletonize(kernel, &quadtree);
    ie_Mat domain(TEST_SIZE * TEST_SIZE * solution_dimension, 1);

    schur_solve(config, kernel, ie_solver_tools, quadtree, domain_points, &domain);
    std::string filename = "output/bake/sol/" + std::to_string(frame)
                           + ".txt";
    write_solution_to_file(filename, domain, domain_points, solution_dimension);
  }



  // double error;
  // switch (config.pde) {
  //   case ie_solver_config::LAPLACE:
  //     error = laplace_error(domain, config.id_tol, domain_points,
  //                           config.boundary.get());
  //     break;
  //   case ie_solver_config::STOKES:
  //     error = stokes_error(domain, config.id_tol, domain_points,
  //                          config.boundary.get());
  //     break;
  // }

  // if (!config.testing) {
  //   write_solution_to_file("output/data/ie_solver_solution.txt", domain,
  //                          domain_points, solution_dimension);
  // }
  // return error;
  return 0;
}


// TODO(John) put write functions in another file or something
void write_boundary_to_file(const std::vector<double>& points) {
  assert(points.size() % 2 == 0);
  std::ofstream output;
  output.open("output/data/ie_solver_boundary.txt");
  if (output.is_open()) {
    for (unsigned int i = 0; i < points.size(); i += 2) {
      output << points[i] << "," << points[i + 1] << std::endl;
    }
    output.close();
  } else {
    LOG::ERROR("Failed to open boundary output file!");
  }
}

void write_potential_to_file() {
  std::cout << "unimplemented" << std::endl;
}


void write_times_to_files(int* scale_n, const std::vector<double>& n_times,
                          double * scale_eps,
                          const std::vector<double>& eps_times) {
  std::ofstream n_output, e_output;
  n_output.open("output/data/ie_solver_n_scaling.txt");
  e_output.open("output/data/ie_solver_e_scaling.txt");

  if (n_output.is_open()) {
    for (unsigned int i = 0; i < n_times.size(); i++) {
      n_output << scale_n[i] << "," << n_times[i] << std::endl;
    }
    n_output.close();
  } else {
    printf("Failed to open n output file!\n");
  }

  if (e_output.is_open()) {
    for (unsigned int i = 0; i < eps_times.size(); i++) {
      e_output << scale_eps[i] << "," << eps_times[i] << std::endl;
    }
    e_output.close();
  } else {
    printf("Failed to open e output file!\n");
  }
}


void write_solution_to_file(const std::string & filename, const ie_Mat & domain,
                            const std::vector<double>&
                            domain_points, int solution_dimension) {
  assert(domain.height() > 0 && domain.width() == 1);
  std::ofstream output;
  output.open(filename);

  int points_index = 0;

  if (output.is_open()) {
    for (unsigned int i = 0; i < domain.height(); i += solution_dimension) {
      output << domain_points[points_index] << "," <<
             domain_points[points_index + 1] << ",";
      points_index += 2;  // depends on domain dimension
      for (int dim = 0; dim < solution_dimension; dim++) {
        output << domain.get(i + dim, 0);
        if (dim < solution_dimension - 1) {
          output << ",";
        }
      }
      output << std::endl;
    }
    output.close();
  } else {
    printf("Failed to open solution output file!\n");
  }
}


void get_domain_points(std::vector<double>* points, double min, double max) {
  for (int i = 0; i < TEST_SIZE; i++) {
    double x = min + ((i + 0.0) / TEST_SIZE) * (max - min);
    for (int j = 0; j < TEST_SIZE; j++) {
      double y = min + ((j + 0.0) / TEST_SIZE) * (max - min);
      points->push_back(x);
      points->push_back(y);
    }
  }
}


double laplace_error(const ie_Mat & domain, double id_tol,
                     const std::vector<double>& domain_points,
                     Boundary * boundary) {
  double max = 0;
  double diff_norm = 0;
  double avg = 0;
  double norm_of_true = 0;
  for (unsigned int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    Vec2 x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    double potential;
    switch (boundary->boundary_condition) {
      case Boundary::BoundaryCondition::SINGLE_ELECTRON:
        potential = log(sqrt(pow(x0 + 2, 2) + pow(x1 + 2, 2))) / (2 * M_PI);
        break;
      case Boundary::BoundaryCondition::ALL_ONES:
        potential = 1.0;
        break;
      case Boundary::BoundaryCondition::BUMP_FUNCTION: {
        std::cout << "Error: check Laplace called on Bump BC;"
                  << " no analytic solution known to check against."
                  << std::endl;
        break;
      }
      case Boundary::BoundaryCondition::STOKES:
        std::cout << "Error: check Laplace called on Stokes BC." << std::endl;
        break;
    }
    if (std::isnan(domain.get(i / 2, 0))) {
      continue;
    }
    double diff = std::abs(potential - domain.get(i / 2, 0));
    avg += diff / potential;
    diff_norm += pow(diff, 2);
    norm_of_true += pow(potential, 2);
    max = std::max(max, diff / potential);
  }
  avg /= domain_points.size();
  diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
  return diff_norm;
}


double stokes_error(const ie_Mat & domain_solution, double id_tol,
                    const std::vector<double>& domain_points,
                    Boundary * boundary) {
  if (boundary->boundary_shape != Boundary::CIRCLE) {
    std::cout << "Error: cannot check stokes error on non-circle" << std::endl;
    return -1;
  }
  double truth_size = 0;
  double total_diff = 0;
  for (unsigned int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    Vec2 x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    Vec2 center(0.5, 0.5);



    Vec2 r = x - center;
    Vec2 sol = Vec2(domain_solution.get(i, 0), domain_solution.get(i + 1, 0));

    Vec2 truth = Vec2(-r.a[1], r.a[0]);
    truth = truth * (1 / truth.norm());
    double om1 = -30;
    double om2 = 4;
    double r1 = 0.05;
    double r2 = 0.25;
    double c1 = (om2 * pow(r2, 2) - om1 * pow(r1, 2))
                / (pow(r2, 2) - pow(r1, 2));
    double c2 = ((om1 - om2) * pow(r2, 2) * pow(r1, 2))
                / (pow(r2, 2) - pow(r1, 2));

    double truth_length = c1 * r.norm() + (c2 / r.norm());

    truth = truth * truth_length;
    // domain_solution.set(i, 0, truth.a[0]);
    // domain_solution.set(i + 1, 0, truth.a[1]);

    double diff = (truth - sol).norm();
    truth_size += fabs(truth_length);
    total_diff += diff;
  }

  return total_diff / truth_size;
}

// TODO(John) this is missing the input of a Stokes Boundary Condition
int parse_input_into_config(int argc, char** argv, ie_solver_config * config) {
  std::string usage = "\n\tusage: ./ie-solver "
                      "-pde {LAPLACE|STOKES} "
                      "-boundary {CIRCLE|ROUNDED_SQUARE|"
                      "ROUNDED_SQUARE_WITH_BUMP|SQUIGGLY|CUBIC_SPLINE "
                      "-boundary_condition {SINGLE_ELECTRON|ALL_ONES|"
                      "BUMP_FUNCTION} "
                      "-N {number of nodes} "
                      "-e {ID error tolerance} "
                      "{-scaling} {-strong}"
                      "\nOmitting an arg triggers a default value.";
  std::string boundary_name = "CIRCLE";
  std::string pde_name = "LAPLACE";
  std::string boundary_condition_name = "SINGLE_ELECTRON";
  // The default boundary is Circle, set it here.
  config->boundary.reset(new Circle());

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-pde")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "STOKES")) {
          config->pde = ie_solver_config::STOKES;
        } else if (!strcmp(argv[i + 1], "LAPLACE")) {
          config->pde = ie_solver_config::LAPLACE;
        } else {
          LOG::ERROR("Unrecognized pde: " + std::string(argv[i + 1])
                     + "\n Acceptable pdes: STOKES, LAPLACE");
          return 0;
        }
        pde_name = argv[i + 1];
      }
      i++;
    } else if (!strcmp(argv[i], "-N")) {
      if (i < argc - 1) {
        config->num_boundary_points = std::stoi(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-strong")) {
      config->admissibility = ie_solver_config::STRONG;
    } else if (!strcmp(argv[i], "-scaling")) {
      config->scaling = true;
    } else if (!strcmp(argv[i], "-h")) {
      LOG::log_level_ = LOG::LOG_LEVEL::INFO_;
      LOG::INFO(usage);
      return 0;
    } else if (!strcmp(argv[i], "-e")) {
      if (i < argc - 1) {
        config->id_tol = std::stof(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-boundary")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "CIRCLE")) {
          config->boundary.reset(new Circle());
        } else if (!strcmp(argv[i + 1], "ROUNDED_SQUARE")) {
          config->boundary.reset(new RoundedSquare());
        } else if (!strcmp(argv[i + 1], "ROUNDED_SQUARE_WITH_BUMP")) {
          config->boundary.reset(new RoundedSquareWithBump());
        } else if (!strcmp(argv[i + 1], "SQUIGGLY")) {
          config->boundary.reset(new Squiggly());
        }  else if (!strcmp(argv[i + 1], "ANNULUS")) {
          config->boundary.reset(new Annulus());
        } else if (!strcmp(argv[i + 1], "CUBIC_SPLINE")) {
          config->boundary.reset(new CubicSpline());
        } else {
          LOG::ERROR("Unrecognized boundary: " + std::string(argv[i + 1])
                     + usage);
          return 0;
        }
        boundary_name = argv[i + 1];
      }
      i++;
    } else if (!strcmp(argv[i], "-boundary_condition")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "SINGLE_ELECTRON")) {
          config->boundary_condition =
            Boundary::BoundaryCondition::SINGLE_ELECTRON;
        } else if (!strcmp(argv[i + 1], "ALL_ONES")) {
          config->boundary_condition = Boundary::BoundaryCondition::ALL_ONES;
        } else if (!strcmp(argv[i + 1], "BUMP_FUNCTION")) {
          config->boundary_condition =
            Boundary::BoundaryCondition::BUMP_FUNCTION;
        } else {
          LOG::ERROR("Unrecognized boundary_condition: " +
                     std::string(argv[i + 1]) + "\n Acceptable "
                     "boundary_conditions: SINGLE_ELECTRON, ALL_ONES,"
                     "BUMP_FUNCTION");
          return 0;
        }
        boundary_condition_name = argv[i + 1];
      }
      i++;
    } else {
      LOG::ERROR("Unrecognized argument: " + std::string(argv[i]) + usage);
      return 0;
    }
  }

  LOG::INFO("PDE: " + pde_name);
  LOG::INFO("Boundary: " + boundary_name);
  LOG::INFO("Boundary condition: " + boundary_condition_name);
  LOG::INFO("Number of nodes: " + std::to_string(config->num_boundary_points));
  LOG::INFO("ID error tolerance: " + std::to_string(config->id_tol));
  if (config->scaling) {
    LOG::INFO("Scaling run");
  }
  if (config->admissibility == ie_solver_config::STRONG) {
    LOG::INFO("Strong admissibility");
  } else {
    LOG::INFO("Weak admissibility");
  }

  return 1;
  // TODO(John) after parsing, print out input configuration
}

}  // namespace ie_solver
