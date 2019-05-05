// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
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


ie_Mat initialize_U_mat(const ie_solver_config::Pde pde,
                        const std::vector<Hole>& holes,
                        const std::vector<double>& tgt_points) {
  ie_Mat U;
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE: {
      U = ie_Mat(tgt_points.size() / 2, holes.size());
      for (unsigned int i = 0; i < tgt_points.size(); i += 2) {
        for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
          Hole hole = holes[hole_idx];
          Vec2 r = Vec2(tgt_points[i], tgt_points[i + 1]) - hole.center;
          U.set(i / 2, hole_idx, log(r.norm()));
        }
      }
      break;
    }
    case ie_solver_config::Pde::STOKES: {
      U = ie_Mat(tgt_points.size(), 3 * holes.size());
      for (unsigned int i = 0; i < tgt_points.size(); i += 2) {
        for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
          Hole hole = holes[hole_idx];
          Vec2 r = Vec2(tgt_points[i], tgt_points[i + 1]) - hole.center;
          double scale = 1.0 / (4 * M_PI);
          // TODO(John) implement tensor product to improve readability?
          U.set(i, 3 * hole_idx, scale *
                (log(1 / r.norm()) +
                 (1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[0]));
          U.set(i + 1, 3 * hole_idx, scale *
                (
                  (1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[0]));
          U.set(i, 3 * hole_idx + 1, scale *
                (
                  (1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[1]));
          U.set(i + 1, 3 * hole_idx + 1, scale *
                (log(1 / r.norm()) +
                 (1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[1]));

          U.set(i, 3 * hole_idx + 2, r.a[1] * (scale / pow(r.norm(), 2)));
          U.set(i + 1, 3 * hole_idx + 2, -r.a[0] * (scale / pow(r.norm(), 2)));
        }
      }
    }
  }
  return U;
}


ie_Mat initialize_Psi_mat(const ie_solver_config::Pde pde,
                          const std::vector<Hole>& holes, Boundary * boundary) {
  ie_Mat Psi;
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE: {
      Psi = ie_Mat(holes.size(), boundary->points.size() / 2);
      for (unsigned int i = 0; i < boundary->points.size(); i += 2) {
        Vec2 x = Vec2(boundary->points[i], boundary->points[i + 1]);
        for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
          Hole hole = holes[hole_idx];
          if ((x - hole.center).norm() < hole.radius + 1e-8) {
            Psi.set(hole_idx, i / 2, boundary->weights[i / 2]);
            break;
          }
        }
      }
      break;
    }
    case ie_solver_config::Pde::STOKES: {
      Psi = ie_Mat(3 * holes.size(), boundary->points.size());
      for (unsigned int i = 0; i < boundary->points.size(); i += 2) {
        Vec2 x = Vec2(boundary->points[i], boundary->points[i + 1]);
        for (unsigned int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
          Hole hole = holes[hole_idx];
          if ((x - hole.center).norm() < hole.radius + 1e-8) {
            Psi.set(3 * hole_idx, i, boundary->weights[i / 2]);
            Psi.set(3 * hole_idx + 1, i + 1, boundary->weights[i / 2]);
            Psi.set(3 * hole_idx + 2, i, boundary->weights[i / 2]*x.a[1]);
            Psi.set(3 * hole_idx + 2, i + 1, -boundary->weights[i / 2]*x.a[0]);
            break;
          }
        }
      }
    }
  }
  return Psi;
}


void schur_solve(const QuadTree & quadtree, const ie_Mat & U,
                 const ie_Mat & Psi,
                 const ie_Mat & f, const ie_Mat & K_domain,
                 const ie_Mat & U_forward,  ie_Mat * solution) {
  ie_Mat ident(U.width(), U.width()), alpha(U.width(), 1),
         Dinv_U(U.height(), U.width()), Psi_Dinv_U(U.width(), U.width()),
         Dinv_u(U.height(), 1), Psi_Dinv_u(U.width(), 1),
         U_alpha(U.height(), 1), mu(K_domain.width(),  1),
         U_forward_alpha(solution->height(), 1);
  if (U.width() == 0) {
    quadtree.solve(&mu, f);
    ie_Mat::gemv(NORMAL, 1., K_domain, mu, 0., solution);
    return;
  }
  quadtree.multiply_connected_solve(&mu, &alpha, f);
  ie_Mat::gemv(NORMAL, 1., K_domain, mu, 0., solution);
  ie_Mat::gemv(NORMAL, 1., U_forward, alpha, 0., &U_forward_alpha);
  (*solution) += U_forward_alpha;
  return;
}


ie_Mat boundary_integral_solve(const ie_solver_config & config,
                               Boundary * boundary, QuadTree * quadtree,
                               const std::vector<double>& domain_points,
                               bool is_time_trial) {
  // Consider making init instead of constructor for readability
  IeSolverTools ie_solver_tools(config.id_tol, config.is_strong_admissibility,
                                config.solution_dimension,
                                config.domain_dimension);

  Kernel kernel;
  kernel.load(boundary, domain_points, config.pde,
              config.solution_dimension, config.domain_dimension);
  ie_solver_tools.skeletonize(kernel, quadtree);
  if (is_time_trial) {
    return ie_Mat(0, 0);
  }
  ie_Mat f = boundary->boundary_values;
  int num_holes = boundary->holes.size();
  ie_Mat U = initialize_U_mat(config.pde, boundary->holes, boundary->points);
  ie_Mat Psi = initialize_Psi_mat(config.pde, boundary->holes,
                                  boundary);
  ie_Mat U_forward = initialize_U_mat(config.pde, boundary->holes,
                                      domain_points);
  ie_Mat K_domain(config.domain_size * config.domain_size *
                  config.solution_dimension,
                  boundary->weights.size() * config.solution_dimension);
  Initialization init;
  init.InitializeDomainKernel(&K_domain, domain_points,
                              config.domain_size, kernel,
                              config.solution_dimension);

  for (unsigned int i = 0; i < domain_points.size(); i += 2) {
    Vec2 point = Vec2(domain_points[i], domain_points[i + 1]);
    if (!boundary->is_in_domain(point)) {
      for (unsigned int hole_idx = 0; hole_idx < num_holes; hole_idx++) {
        switch (config.pde) {
          case ie_solver_config::Pde::LAPLACE: {
            U_forward.set(i / 2, hole_idx, 0);
            break;
          }
          case ie_solver_config::Pde::STOKES: {
            U_forward.set(i, 3 * hole_idx, 0);
            U_forward.set(i + 1, 3 * hole_idx, 0);
            U_forward.set(i, 3 * hole_idx + 1, 0);
            U_forward.set(i + 1, 3 * hole_idx + 1, 0);
            U_forward.set(i, 3 * hole_idx + 2, 0);
            U_forward.set(i + 1, 3 * hole_idx + 2, 0);
            break;
          }
        }
      }
    }
  }

  ie_Mat domain_solution(config.domain_size * config.domain_size *
                         config.solution_dimension, 1);
  quadtree->U = U;
  quadtree->Psi = Psi;
  schur_solve(*quadtree, U, Psi, f, K_domain, U_forward, &domain_solution);
  return domain_solution;
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


void get_domain_points(unsigned int domain_size, std::vector<double>* points,
                       double min, double max) {
  for (unsigned int i = 0; i < domain_size; i++) {
    double x = min + ((i + 0.0) / domain_size) * (max - min);
    for (int j = 0; j < domain_size; j++) {
      double y = min + ((j + 0.0) / domain_size) * (max - min);
      points->push_back(x);
      points->push_back(y);
    }
  }
}


double laplace_error(const ie_Mat & domain, double id_tol,
                     const std::vector<double>& domain_points,
                     Boundary * boundary) {
  if (boundary->holes.size() > 0) {
    std::cout <<
              "Error: laplace error not currently calculated for multiply"
              << " connected domain." << std::endl;
    return -1;
  }
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
  if (boundary->boundary_shape != Boundary::ANNULUS) {
    std::cout << "Error: cannot currently check stokes error on non-annulus" <<
              std::endl;
    return -1;
  }
  if (boundary->holes.size() != 1) {
    std::cout << "Error: can only check error on boundary with one hole" <<
              std::endl;
    return -1;
  } else if (boundary->holes[0].center.a[0] != 0.5
             || boundary->holes[0].center.a[1] != 0.5
             || boundary->holes[0].radius != 0.05) {
    std::cout << "Error: can only check error on boundary with hole at center "
              << "and radius 0.05" << std::endl;
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
                      "ROUNDED_SQUARE_WITH_BUMP|SQUIGGLY|CUBIC_SPLINE|ANNULUS "
                      "-boundary_condition {SINGLE_ELECTRON|ALL_ONES|"
                      "BUMP_FUNCTION|STOKES} "
                      "-N {number of nodes} "
                      "-D {side length of square domain grid} "
                      "-e {ID error tolerance} "
                      "{-scaling} {-strong} {-animation}"
                      "\nOmitting an arg triggers a default value.";
  std::string boundary_name = "CIRCLE";
  std::string pde_name = "LAPLACE";
  std::string boundary_condition_name = "SINGLE_ELECTRON";
  // The default boundary is Circle, set it here.

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-pde")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "STOKES")) {
          config->pde = ie_solver_config::STOKES;
        } else if (!strcmp(argv[i + 1], "LAPLACE")) {
          config->pde = ie_solver_config::LAPLACE;
        } else {
          LOG::ERROR("Unrecognized pde: " + std::string(argv[i + 1])
                     + "\n Acceptable pdes: LAPLACE, STOKES");
          return -1;
        }
        pde_name = argv[i + 1];
      }
      i++;
    } else if (!strcmp(argv[i], "-N")) {
      if (i < argc - 1) {
        config->num_boundary_points = std::stoi(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-D")) {
      if (i < argc - 1) {
        config->domain_size = std::stoi(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-strong")) {
      config->is_strong_admissibility = true;
    } else if (!strcmp(argv[i], "-scaling")) {
      config->scaling = true;
    } else if (!strcmp(argv[i], "-animation")) {
      config->animation = true;
    } else if (!strcmp(argv[i], "-h")) {
      LOG::log_level_ = LOG::LOG_LEVEL::INFO_;
      LOG::INFO(usage);
      return -1;
    } else if (!strcmp(argv[i], "-e")) {
      if (i < argc - 1) {
        config->id_tol = std::stof(argv[i + 1]);
      }
      i++;
    } else if (!strcmp(argv[i], "-boundary")) {
      if (i < argc - 1) {
        if (!strcmp(argv[i + 1], "CIRCLE")) {
          config->boundary_shape = Boundary::BoundaryShape::CIRCLE;
        } else if (!strcmp(argv[i + 1], "ROUNDED_SQUARE")) {
          config->boundary_shape = Boundary::BoundaryShape::ROUNDED_SQUARE;
        } else if (!strcmp(argv[i + 1], "ROUNDED_SQUARE_WITH_BUMP")) {
          config->boundary_shape =
            Boundary::BoundaryShape::ROUNDED_SQUARE_WITH_BUMP;
        } else if (!strcmp(argv[i + 1], "SQUIGGLY")) {
          config->boundary_shape = Boundary::BoundaryShape::SQUIGGLY;
        }  else if (!strcmp(argv[i + 1], "ANNULUS")) {
          config->boundary_shape = Boundary::BoundaryShape::ANNULUS;
        } else if (!strcmp(argv[i + 1], "CUBIC_SPLINE")) {
          config->boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
        } else {
          LOG::ERROR("Unrecognized boundary: " + std::string(argv[i + 1])
                     + usage);
          return -1;
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
        } else if (!strcmp(argv[i + 1], "STOKES")) {
          config->boundary_condition =
            Boundary::BoundaryCondition::STOKES;
        } else {
          LOG::ERROR("Unrecognized boundary_condition: " +
                     std::string(argv[i + 1]) + "\n Acceptable "
                     "boundary_conditions: SINGLE_ELECTRON, ALL_ONES,"
                     "BUMP_FUNCTION, STOKES");
          return -1;
        }
        boundary_condition_name = argv[i + 1];
      }
      i++;
    } else {
      LOG::ERROR("Unrecognized argument: " + std::string(argv[i]) + usage);
      return -1;
    }
  }

  if (config->boundary_condition
      ==  ie_solver::Boundary::BoundaryCondition::STOKES
      || config->pde == ie_solver::ie_solver_config::STOKES) {
    config->boundary_condition =
      ie_solver::Boundary::BoundaryCondition::STOKES;
    config->pde = ie_solver::ie_solver_config::STOKES;
    boundary_condition_name = "STOKES";
    pde_name = "STOKES";
    config->solution_dimension = 2;
  }

  LOG::INFO("PDE: " + pde_name);
  LOG::INFO("Boundary: " + boundary_name);
  LOG::INFO("Boundary condition: " + boundary_condition_name);
  LOG::INFO("Number of nodes: " + std::to_string(config->num_boundary_points));
  LOG::INFO("Side length of square domain grid: " + std::to_string(
              config->domain_size));
  LOG::INFO("ID error tolerance: " + std::to_string(config->id_tol));
  if (config->scaling) {
    LOG::INFO("Scaling run");
  }
  if (config->animation) {
    LOG::INFO("Animation run");
  }
  if (config->is_strong_admissibility) {
    LOG::INFO("Strong admissibility");
  } else {
    LOG::INFO("Weak admissibility");
  }

  return 0;
}

}  // namespace ie_solver
