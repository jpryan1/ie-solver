// Copyright John Paul Ryan 2019
#include <omp.h>
#include <string.h>
#include <fstream>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include "ie_solver/ie_mat.h"
#include "ie_solver/initialization.h"
#include "ie_solver/boundaries/circle.h"
#include "ie_solver/boundaries/annulus.h"
#include "ie_solver/boundaries/donut.h"
#include "ie_solver/boundaries/ex1boundary.h"
#include "ie_solver/boundaries/ex2boundary.h"
#include "ie_solver/boundaries/ex3boundary.h"
#include "ie_solver/boundaries/cubic_spline.h"
#include "ie_solver/skel_factorization/skel_factorization.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/kernel/kernel.h"
#include "ie_solver/log.h"
#include "ie_solver/linear_solve.h"
#include "gtest/gtest.h"

namespace ie_solver {



void check_solve_err(const ie_solver_config& config, Boundary* boundary) {
  QuadTree quadtree;
  quadtree.initialize_tree(boundary, std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.num_threads);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary, std::vector<double>());
  double max = 0.;
  ie_Mat dense;

  if (boundary->holes.size() > 0
      && config.pde != ie_solver_config::Pde::LAPLACE_NEUMANN) {
    ie_Mat U = initialize_U_mat(config.pde, boundary->holes, boundary->points);
    ie_Mat Psi = initialize_Psi_mat(config.pde, boundary->holes, boundary);
    quadtree.U = U;
    quadtree.Psi = Psi;
    skel_factorization.skeletonize(kernel, &quadtree);

    ie_Mat mu, alpha;
    linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu,
                 &alpha);
    ie_Mat stacked(mu.height() + alpha.height(), 1);
    stacked.set_submatrix(0, mu.height(), 0, 1, mu);
    stacked.set_submatrix(mu.height(), stacked.height(), 0, 1, alpha);

    std::vector<int> all_dofs;
    for (int i = 0; i < config.solution_dimension * boundary->weights.size();
         i++) {
      all_dofs.push_back(i);
    }
    ie_Mat kern = kernel(all_dofs, all_dofs);
    int hole_factor = config.solution_dimension == 2 ? 3 : 1;
    int added = hole_factor * boundary->holes.size();
    dense = ie_Mat(all_dofs.size() + added, all_dofs.size() + added);
    ie_Mat ident(added, added);
    ident.eye(added);
    dense.set_submatrix(0, all_dofs.size(), 0, all_dofs.size(), kern);
    dense.set_submatrix(0, all_dofs.size(), all_dofs.size(), dense.width(), U);
    dense.set_submatrix(all_dofs.size(), dense.height(),
                        0, all_dofs.size(), Psi);
    dense.set_submatrix(all_dofs.size(), dense.height(),
                        all_dofs.size(), dense.width(), -ident);

    ie_Mat fzero_prime = dense * stacked;

    ie_Mat err1 = (fzero_prime(0, mu.height(), 0, 1)
                   - boundary->boundary_values);
    ie_Mat err2 = (fzero_prime(mu.height(), fzero_prime.height(), 0, 1));
    max = err1.max_entry_magnitude();
    EXPECT_LE(std::max(max, err2.max_entry_magnitude()), 50 * config.id_tol);
  } else {
    skel_factorization.skeletonize(kernel, &quadtree);
    ie_Mat mu;
    linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);
    std::vector<int> all_dofs;
    for (int i = 0; i < config.solution_dimension * boundary->weights.size();
         i++) {
      all_dofs.push_back(i);
    }

    ie_Mat err = (kernel(all_dofs, all_dofs) * mu) - boundary->boundary_values;
    EXPECT_LE(err.max_entry_magnitude(), 50 * config.id_tol);
  }
}


ie_Mat get_dom_sol(const ie_solver_config& config,
                   std::vector<double>* domain_points,
                   Boundary* boundary) {
  QuadTree quadtree;
  quadtree.initialize_tree(boundary, std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  get_domain_points(config.domain_size, domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  return boundary_integral_solve(config, &quadtree, *domain_points);
}


TEST(IeSolverTest, LaplaceCircleBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, LaplaceNeumannCircleBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::ALL_ONES;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, StokesCircleBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, LaplaceStarfishBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, LaplaceNeumannStarfishBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::ALL_ONES;
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, StokesStarfishBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}



TEST(IeSolverTest, LaplaceAnnulusBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, LaplaceNeumannAnnulusBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, StokesAnnulusBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}



TEST(IeSolverTest, BigLaplaceCircleBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigLaplaceNeumannCircleBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::ALL_ONES;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigStokesCircleBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigLaplaceStarfishBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigLaplaceNeumannStarfishBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::ALL_ONES;
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigStokesStarfishBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}



TEST(IeSolverTest, BigLaplaceAnnulusBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigLaplaceNeumannAnnulusBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel::IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}


TEST(IeSolverTest, BigStokesAnnulusBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.num_boundary_points *= 8;
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  Kernel:: IMPROVE_CONDITION = 0;
  check_solve_err(config, boundary.get());
}





double laplace_error(const ie_Mat& domain,
                     const std::vector<double>& domain_points,
                     Boundary * boundary) {
  double diff_norm = 0;
  double norm_of_true = 0;
  for (int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    Vec2 x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    double potential = log(sqrt(pow(x0 + 3, 2) + pow(x1 + 2, 2))) / (2 * M_PI);

    if (std::isnan(domain.get(i / 2, 0))) {
      continue;
    }
    double diff = std::abs(potential - domain.get(i / 2, 0));
    diff_norm += pow(diff, 2);
    norm_of_true += pow(potential, 2);
  }
  diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
  return diff_norm;
}


TEST(IeSolverTest, LaplaceCircleAnalyticAgreementElectron) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  std::vector<double> domain_points;
  ie_Mat sol = get_dom_sol(config, &domain_points,
                           boundary.get());
  double err = laplace_error(sol, domain_points, boundary.get());
  EXPECT_LE(std::abs(err), 10 * config.id_tol);
}


TEST(IeSolverTest, LaplaceAnnulusAnalyticAgreementElectron) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  std::vector<double> domain_points;
  ie_Mat sol = get_dom_sol(config, &domain_points,
                           boundary.get());
  double err = laplace_error(sol, domain_points, boundary.get());
  EXPECT_LE(std::abs(err), 10 * config.id_tol);
}

double laplace_neumann_error(const ie_Mat& domain,
                             const std::vector<double>& domain_points,
                             Boundary * boundary) {
  ie_Mat no_ln = domain;
  for (int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    Vec2 x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      no_ln.set(i / 2, 0, 0);
      continue;
    }
    if (std::isnan(domain.get(i / 2, 0))) {
      no_ln.set(i / 2, 0, 0);
      continue;
    }
    double r = sqrt(pow(x0 - 0.5, 2) + pow(x1 - 0.5, 2));
    no_ln.set(i / 2, 0, domain.get(i / 2, 0) - log(r));
  }

  double avg = 0;
  int avg_count = 0;
  // Now we take the avg of the nonzeros, subtract the nonzeros by it,
  // and expect a smalllll norm
  for (int i = 0; i < no_ln.height(); i++) {
    if (no_ln.get(i, 0) != 0 && !std::isnan(no_ln.get(i, 0))) {
      avg_count++;
      avg += no_ln.get(i, 0);
    }
  }
  avg /= avg_count;
  double max = 0;
  for (int i = 0; i < no_ln.height(); i++) {
    if (no_ln.get(i, 0) != 0 && !std::isnan(no_ln.get(i, 0))) {
      max = std::max(max, no_ln.get(i, 0) - avg);
    }
  }
  return max;
}


double stokes_error(const ie_Mat& domain,
                    const std::vector<double>& domain_points,
                    Boundary * boundary) {
  double max_err = 0;
  for (int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    Vec2 x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    Vec2 center(0.5, 0.5);

    Vec2 r = x - center;

    Vec2 sol = Vec2(domain.get(i, 0),
                    domain.get(i + 1, 0));
    Vec2 truth = Vec2(-r.a[1], r.a[0]);
    switch (boundary->boundary_shape) {
      case Boundary::BoundaryShape::CIRCLE:
        truth = truth * (r.norm() / truth.norm()) * 4.;
        break;
      case Boundary::BoundaryShape::DONUT:
      default:
        double p = (-1. / r.norm()) + 2 * r.norm();
        truth = truth * (1. / truth.norm()) * p;
        break;
    }

    max_err = std::max(max_err, (truth - sol).norm());
  }

  return max_err;
}


TEST(IeSolverTest, LaplaceNeumannDonutAnalyticAgreement) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_shape = Boundary::BoundaryShape::DONUT;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Donut());

  boundary->initialize(config.num_boundary_points, DEFAULT);
  std::vector<double> domain_points;
  ie_Mat sol = get_dom_sol(config, &domain_points,
                           boundary.get());
  double err = laplace_neumann_error(sol, domain_points, boundary.get());
  EXPECT_LE(std::abs(err), 10 * config.id_tol);
}


TEST(IeSolverTest, StokesCircleAnalyticAgreementTangent) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.solution_dimension = 2;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  std::vector<double> domain_points;
  ie_Mat sol = get_dom_sol(config, &domain_points,
                           boundary.get());
  double err = stokes_error(sol, domain_points, boundary.get());
  EXPECT_LE(std::abs(err), 10 * config.id_tol);
}

TEST(IeSolverTest, StokesDonutAnalyticAgreementTangent) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.solution_dimension = 2;
  // Need more nodes for this one.
  config.num_boundary_points *= 4;
  config.pde = ie_solver_config::Pde::STOKES;
  config.boundary_condition = BoundaryCondition::TANGENT_VEC;
  config.boundary_shape = Boundary::BoundaryShape::DONUT;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Donut());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);
  std::vector<double> domain_points;
  ie_Mat sol = get_dom_sol(config, &domain_points,
                           boundary.get());
  double err = stokes_error(sol, domain_points, boundary.get());
  EXPECT_LE(std::abs(err), 10 * config.id_tol);
}

ie_solver_config get_experiment_one_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = pow(2, 12);
  config.domain_size = 20;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX1;
  return config;
}

TEST(IeSolverTest, Ex1UpdateLosesNoAcc) {
  ie_solver_config config = get_experiment_one_config();
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex1Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex1Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 BoundaryCondition::DEFAULT);

  int FRAME_CAP = 10;
  for (int frame = 0; frame < FRAME_CAP; frame++) {
    double ang = (frame / (0.0 + FRAME_CAP)) * 2 * M_PI;

    perturbed_boundary->perturbation_parameters[0] = ang;
    perturbed_boundary->initialize(config.num_boundary_points,
                                   config.boundary_condition);
    quadtree.perturb(*perturbed_boundary.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree,
                      domain_points);
    QuadTree fresh;
    fresh.initialize_tree(perturbed_boundary.get(), std::vector<double>(), 2,  2);
    ie_Mat new_sol = boundary_integral_solve(config, &fresh,
                     domain_points);
    ASSERT_LE((new_sol - solution).max_entry_magnitude(), 100 * config.id_tol);
  }
}


ie_solver_config get_experiment_two_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = pow(2, 10);
  config.domain_size = 20;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX2;
  return config;
}


// Note that increasing N causes this test to fail, likely due to conditioning
TEST(IeSolverTest, Ex2UpdateLosesNoAcc) {
  ie_solver_config config = get_experiment_two_config();
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 BoundaryCondition::DEFAULT);
  int FRAME_CAP = 10;
  for (int frame = 0; frame < FRAME_CAP; frame++) {
    int rand_idx = floor(8 * (rand() / (0. + RAND_MAX)));
    perturbed_boundary->perturbation_parameters[rand_idx] = 0.35
        + 0.3 * (rand() / (0. + RAND_MAX));
    perturbed_boundary->initialize(config.num_boundary_points,
                                   config.boundary_condition);

    quadtree.perturb(*perturbed_boundary.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree,
                      domain_points);
    QuadTree fresh;
    fresh.initialize_tree(perturbed_boundary.get(), std::vector<double>(), 2,  2);
    ie_Mat new_sol = boundary_integral_solve(config, &fresh,
                     domain_points);
    // Allow 10* more error due to conditioning
    ASSERT_LE((new_sol - solution).max_entry_magnitude(), 300 * config.id_tol);
  }
}



ie_solver_config get_experiment_three_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.num_boundary_points = pow(2, 10);
  config.domain_size = 20;
  config.domain_dimension = 2;
  config.solution_dimension = 1;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX3;
  return config;
}


TEST(IeSolverTest, Ex3UpdateLosesNoAcc) {
  ie_solver_config config = get_experiment_three_config();
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);

  int FRAME_CAP = 5;
  for (int frame1 = 0; frame1 < FRAME_CAP; frame1++) {
    double ang1 = (frame1 / (0.0 + FRAME_CAP)) * 4 * M_PI / 5.;
    for (int frame2 = 0; frame2 < FRAME_CAP; frame2++) {
      double ang2 = ((frame2 / (0.0 + FRAME_CAP)) * 4 * M_PI / 5.) + M_PI;
      perturbed_boundary->perturbation_parameters[0] = ang1;
      perturbed_boundary->perturbation_parameters[1] = ang2;

      perturbed_boundary->initialize(config.num_boundary_points,
                                     config.boundary_condition);

      quadtree.perturb(*perturbed_boundary.get());
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);
      QuadTree fresh;
      fresh.initialize_tree(perturbed_boundary.get(), std::vector<double>(), 1,  2);
      ie_Mat new_sol = boundary_integral_solve(config, &fresh,
                       domain_points);
      ASSERT_LE((new_sol - solution).max_entry_magnitude(), 100 * config.id_tol);
    }
  }
}


TEST(IeSolverTest, TreeCopyGivesSameAnswer) {
  ie_solver_config config = get_experiment_two_config();
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);

  QuadTree fresh;

  quadtree.copy_into(&fresh);
  ie_Mat new_sol = boundary_integral_solve(config, &fresh,
                   domain_points);
  ASSERT_LE((new_sol - solution).max_entry_magnitude(), 1e-15);
}



}  // namespace ie_solver
