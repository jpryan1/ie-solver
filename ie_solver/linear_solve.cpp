// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "ie_solver/ie_mat.h"
#include "ie_solver/initialization.h"
#include "ie_solver/skel_factorization/skel_factorization.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/kernel/kernel.h"
#include "ie_solver/linear_solve.h"
#include "ie_solver/log.h"
#include "ie_solver/boundaries/boundary.h"
#include "ie_solver/boundaries/circle.h"
#include "ie_solver/boundaries/rounded_square.h"
#include "ie_solver/boundaries/rounded_square_with_bump.h"
#include "ie_solver/boundaries/squiggly.h"
#include "ie_solver/boundaries/annulus.h"
#include "ie_solver/boundaries/cubic_spline.h"
#include "ie_solver/boundaries/ex1boundary.h"
#include "ie_solver/boundaries/ex3boundary.h"


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

void check_factorization_against_kernel(const Kernel& kernel, const SkelFactorization&
    skel_factorization, QuadTree* tree) {
  int check_size = 100;
  // This ensures that operations know what the remaining skels are.

  tree->remove_inactive_dofs_at_all_boxes();

  int dofs = tree->boundary->points.size() / 2;
  // Take a random 100x100 submatrix of A-A^hat, estimate 2-norm by power method
  std::vector<unsigned int> rand_x_indices, rand_y_indices;
  for (int i = 0; i < check_size; i++) {
    rand_x_indices.push_back(rand() % dofs);
    rand_y_indices.push_back(rand() % dofs);
  }

  ie_Mat A = kernel(rand_x_indices, rand_y_indices);
  ie_Mat A_hat(check_size, check_size);
  ie_Mat basis(dofs, 1);
  ie_Mat b(dofs, 1);

  for (unsigned int y = 0; y < rand_y_indices.size(); y++) {
    int rand_y_idx = rand_y_indices[y];
    for (int i = 0; i < dofs; i++) {
      basis.set(i, 0, 0);
    }
    basis.set(rand_y_idx, 0, 1.0);
    skel_factorization.sparse_matvec(*tree, basis, &b);
    for (unsigned int x = 0; x < rand_x_indices.size(); x++) {
      int rand_x_idx = rand_x_indices[x];
      A_hat.set(x, y, b.get(rand_x_idx, 0));
    }
  }
  double truenorm = A.frob_norm();
  A -= A_hat;
  std::cout << "Total error: " << 100.0 * A.frob_norm() / truenorm << "%" <<
            std::endl;
}


void schur_solve(const SkelFactorization& skel_factorization,
                 const QuadTree & quadtree, const ie_Mat & U,
                 const ie_Mat & Psi,
                 const ie_Mat & f, const ie_Mat & K_domain,
                 const ie_Mat & U_forward,  ie_Mat * solution) {
  ie_Mat ident(U.width(), U.width()), alpha(U.width(), 1),
         Dinv_U(U.height(), U.width()), Psi_Dinv_U(U.width(), U.width()),
         Dinv_u(U.height(), 1), Psi_Dinv_u(U.width(), 1),
         U_alpha(U.height(), 1), mu(K_domain.width(),  1),
         U_forward_alpha(solution->height(), 1);
  if (U.width() == 0) {
    skel_factorization.solve(quadtree, &mu, f);
    ie_Mat::gemv(NORMAL, 1., K_domain, mu, 0., solution);
    return;
  }
  skel_factorization.multiply_connected_solve(quadtree, &mu, &alpha, f);

  ie_Mat::gemv(NORMAL, 1., K_domain, mu, 0., solution);
  ie_Mat::gemv(NORMAL, 1., U_forward, alpha, 0., &U_forward_alpha);
  (*solution) += U_forward_alpha;
  return;
}

void bie_time_trial(const ie_solver_config & config,
                    QuadTree * quadtree, double* avg_skel_time,
                    double* avg_solve_time) {
  Boundary* boundary = quadtree->boundary;
  // Consider making init instead of constructor for readability
  SkelFactorization skel_factorization(config.id_tol, config.is_strong_admissibility,
                                config.solution_dimension,
                                config.domain_dimension);

  Kernel kernel;
  kernel.load(boundary, std::vector<double>(), config.pde,
              config.solution_dimension, config.domain_dimension);
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
  SkelFactorization skel_factorization(config.id_tol, config.is_strong_admissibility,
                                config.solution_dimension,
                                config.domain_dimension);

  Kernel kernel;
  kernel.load(boundary, domain_points, config.pde,
              config.solution_dimension, config.domain_dimension);

  skel_factorization.skeletonize(kernel, quadtree);

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
  skel_factorization.U = U;
  skel_factorization.Psi = Psi;
  schur_solve(skel_factorization, *quadtree, U, Psi, f, K_domain, U_forward,
              &domain_solution);

  return domain_solution;
  
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



}  // namespace ie_solver
