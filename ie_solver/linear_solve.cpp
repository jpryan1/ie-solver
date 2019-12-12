// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include "ie_solver/ie_mat.h"
#include "ie_solver/initialization.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/linear_solve.h"
#include "ie_solver/log.h"


namespace ie_solver {


ie_Mat initialize_U_mat(const ie_solver_config::Pde pde,
                        const std::vector<Hole>& holes,
                        const std::vector<double>& tgt_points) {
  ie_Mat U;
  if (holes.size() == 0) return ie_Mat(0, 0);
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
                ((1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[0]));
          U.set(i, 3 * hole_idx + 1, scale *
                ((1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[1]));
          U.set(i + 1, 3 * hole_idx + 1, scale *
                (log(1 / r.norm()) +
                 (1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[1]));
          U.set(i, 3 * hole_idx + 2, r.a[1] * (scale / pow(r.norm(), 2)));
          U.set(i + 1, 3 * hole_idx + 2, -r.a[0] * (scale / pow(r.norm(), 2)));
        }
      }
      break;
    }
    case ie_solver_config::Pde::LAPLACE_NEUMANN: {
      U = ie_Mat(0, 0);
      break;
    }
  }
  return U;
}


ie_Mat initialize_Psi_mat(const ie_solver_config::Pde pde,
                          const std::vector<Hole>& holes, Boundary * boundary) {
  if (holes.size() == 0) return ie_Mat(0, 0);
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
      break;
    }
    case ie_solver_config::Pde::LAPLACE_NEUMANN: {
      Psi = ie_Mat(0, 0);
      break;
    }
  }
  return Psi;
}


void linear_solve(const SkelFactorization& skel_factorization,
                  const QuadTree& quadtree, const ie_Mat& f, ie_Mat* mu,
                  ie_Mat* alpha) {
  double start = omp_get_wtime();
  *mu = ie_Mat(quadtree.boundary->weights.size() *
               quadtree.solution_dimension, 1);
  if (alpha == nullptr) {
    skel_factorization.solve(quadtree, mu, f);
  } else {
    *alpha = ie_Mat(quadtree.U.width(), 1);
    skel_factorization.multiply_connected_solve(quadtree, mu, alpha, f);
  }
  double end = omp_get_wtime();
  std::cout << "timing: linear_solve " << (end - start) << std::endl;
}

// void linear_solve(const SkelFactorization& skel_factorization,
//                   const QuadTree& quadtree, const ie_Mat& f, ie_Mat* mu,
//                   double* c) {
//   *mu = ie_Mat(quadtree.boundary->weights.size(), 1);
//   ie_Mat C(1, 1);
//   skel_factorization.multiply_connected_solve(quadtree, mu, &C, f);
//   *c = C.get(0, 0);
// }

void schur_solve(const SkelFactorization & skel_factorization,
                 const QuadTree & quadtree, const ie_Mat & U,
                 const ie_Mat & Psi,
                 const ie_Mat & f, const ie_Mat & K_domain,
                 const ie_Mat & U_forward,  ie_Mat * solution) {
  ie_Mat mu;
  if (U.width() == 0) {
    linear_solve(skel_factorization, quadtree, f, &mu);
    *solution = K_domain*mu;
  } else {
    ie_Mat alpha;
    linear_solve(skel_factorization, quadtree, f, &mu, &alpha);
    *solution = (K_domain*mu) + (U_forward*alpha);
  }
}

void bie_time_trial(const ie_solver_config& config,
                    QuadTree * quadtree, double * avg_skel_time,
                    double * avg_solve_time) {
  Boundary* boundary = quadtree->boundary;
  // Consider making init instead of constructor for readability
  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       DEFAULT_NUM_THREADS);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary, std::vector<double>());
  ie_Mat f = boundary->boundary_values;
  ie_Mat mu(f.height(), 1);
  double skel_time = 0;
  double solve_time = 0;

  for (int i = 0; i < 10; i++) {
    double start = omp_get_wtime();
    skel_factorization.skeletonize(kernel, quadtree);
    double end = omp_get_wtime();
    skel_time += (end - start);
    start = omp_get_wtime();
    skel_factorization.solve(*quadtree, &mu, f);
    end = omp_get_wtime();
    solve_time += (end - start);
    quadtree->reset();
  }
  *avg_skel_time = (skel_time) / 10.0;
  *avg_solve_time = (solve_time) / 10.0;
}


ie_Mat boundary_integral_solve(const ie_solver_config & config,
                               QuadTree * quadtree,
                               const std::vector<double>& domain_points) {
  Boundary* boundary = quadtree->boundary;
  // Consider making init instead of constructor for readability
  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.num_threads);
  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary, domain_points);

  // Domain kernel init takes less time than skel, init in background
  ie_Mat K_domain((domain_points.size() / 2)*
                  config.solution_dimension,
                  boundary->weights.size() * config.solution_dimension);
  // std::thread init_domain_kernel(&Initialization::InitializeDomainKernel,
  //                                &K_domain, domain_points,
  //                                kernel, config.solution_dimension);
  double init_start = omp_get_wtime();
  Initialization::InitializeDomainKernel(
    &K_domain, domain_points,
    kernel, config.solution_dimension);
  double init_end = omp_get_wtime();
  std::cout << "timing: init " << init_end - init_start << std::endl;
  // std::vector<unsigned int> all_inds;
  // for (unsigned int i = 0;
  //      i < boundary->points.size() / (3 - config.solution_dimension); i++) {
  //   all_inds.push_back(i);
  // }
  // ie_Mat all = kernel(all_inds, all_inds);
  // std::cout << "Condition number: " << all.condition_number() << std::endl;

  ie_Mat f = boundary->boundary_values;
  int num_holes = boundary->holes.size();

  ie_Mat U = initialize_U_mat(config.pde, boundary->holes, boundary->points);
  ie_Mat Psi = initialize_Psi_mat(config.pde, boundary->holes,
                                  boundary);
  ie_Mat U_forward = initialize_U_mat(config.pde, boundary->holes,
                                      domain_points);

  // Zero out points outside the domain
  if (config.pde == ie_solver_config::Pde::LAPLACE
      || config.pde == ie_solver_config::Pde::STOKES) {
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
            case ie_solver_config::Pde::LAPLACE_NEUMANN:
              break;  // just suppressing a compiler warning.
          }
        }
      }
    }
  }
  ie_Mat domain_solution((domain_points.size() / 2)*
                         config.solution_dimension, 1);
  quadtree->U = U;
  quadtree->Psi = Psi;
  skel_factorization.skeletonize(kernel, quadtree);

  // init_domain_kernel.join();

  schur_solve(skel_factorization, *quadtree, U, Psi, f, K_domain,
              U_forward, &domain_solution);

  return domain_solution;
}


void get_domain_points(unsigned int domain_size, std::vector<double>* points,
                       double min, double max) {
  for (unsigned int i = 0; i < domain_size; i++) {
    double x = min + ((i + 0.0) / (domain_size - 1)) * (max - min);
    for (int j = 0; j < domain_size; j++) {
      double y = min + ((j + 0.0) / (domain_size - 1)) * (max - min);
      points->push_back(x);
      points->push_back(y);
    }
  }
}



}  // namespace ie_solver
