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
#include "ie_solver/boundaries/cubic_spline.h"
#include "ie_solver/skel_factorization/skel_factorization.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/kernel/kernel.h"
#include "ie_solver/log.h"
#include "ie_solver/linear_solve.h"
#include "gtest/gtest.h"

namespace ie_solver {

// struct ie_solver_config {
//   enum Pde {
//     LAPLACE,
//     LAPLACE_NEUMANN,
//     STOKES
//   };
//   int num_boundary_points = DEFAULT_NUM_DISCRETIZATION_POINTS;
//   int domain_size = DEFAULT_DOMAIN_SIZE;
//   int domain_dimension = 2;
//   int solution_dimension = 1;
//   double id_tol = DEFAULT_ID_TOL;
//   Pde pde = LAPLACE;
//   bool is_strong_admissibility = false;
//   BoundaryCondition boundary_condition =
//     BoundaryCondition::SINGLE_ELECTRON;
//   Boundary::BoundaryShape boundary_shape =
//     Boundary::BoundaryShape::CIRCLE;
//   bool scaling = false;
//   bool animation = false;
// };

TEST(IeSolverTest, LaplaceCircleBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }

  ie_Mat dense = kernel(all_dofs, all_dofs);

  ie_Mat f_prime(mu.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
  ie_Mat err = (f_prime - boundary->boundary_values);
  double max = 0;
  for (int i = 0; i < err.height(); i++) {
    max = std::max(max, std::abs(err.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 100000);
  EXPECT_LE(max, 100 * config.id_tol);
}


TEST(IeSolverTest, LaplaceNeumannCircleBackwardError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.boundary_condition = BoundaryCondition::ALL_ONES;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat dense = kernel(all_dofs, all_dofs);
  ie_Mat f_prime(mu.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
  ie_Mat err = (f_prime - boundary->boundary_values);
  double max = 0;
  for (int i = 0; i < err.height(); i++) {
    max = std::max(max, std::abs(err.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 100000);
  EXPECT_LE(max, 100 * config.id_tol);
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

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < 2 * boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat dense = kernel(all_dofs, all_dofs);
  ie_Mat f_prime(mu.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
  ie_Mat err = (f_prime - boundary->boundary_values);
  double max = 0;
  for (int i = 0; i < err.height(); i++) {
    max = std::max(max, std::abs(err.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 100000);
  EXPECT_LE(max, 100 * config.id_tol);
}


TEST(IeSolverTest, LaplaceStarfishBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.boundary_shape = Boundary::BoundaryShape::CUBIC_SPLINE;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat dense = kernel(all_dofs, all_dofs);
  ie_Mat f_prime(mu.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
  ie_Mat err = (f_prime - boundary->boundary_values);
  double max = 0;
  for (int i = 0; i < err.height(); i++) {
    max = std::max(max, std::abs(err.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 100000);
  EXPECT_LE(max, 100 * config.id_tol);
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

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat dense = kernel(all_dofs, all_dofs);
  ie_Mat f_prime(mu.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
  ie_Mat err = (f_prime - boundary->boundary_values);
  double max = 0;
  for (int i = 0; i < err.height(); i++) {
    max = std::max(max, std::abs(err.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 100000);
  EXPECT_LE(max, 100 * config.id_tol);
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

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < 2 * boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat dense = kernel(all_dofs, all_dofs);
  ie_Mat f_prime(mu.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
  ie_Mat err = (f_prime - boundary->boundary_values);
  double max = 0;
  for (int i = 0; i < err.height(); i++) {
    max = std::max(max, std::abs(err.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 100000);
  EXPECT_LE(max, 100 * config.id_tol);
}



TEST(IeSolverTest, LaplaceAnnulusBackwardError) {
  srand(0);

  ie_solver_config config = ie_solver_config();
  config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());

  ie_Mat U = initialize_U_mat(config.pde, boundary->holes, boundary->points);
  ie_Mat Psi = initialize_Psi_mat(config.pde, boundary->holes,
                                  boundary.get());
  skel_factorization.U = U;
  skel_factorization.Psi = Psi;
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu, alpha;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu,
               &alpha);
  ie_Mat stacked(mu.height() + alpha.height(), 1);
  stacked.set_submatrix(0, mu.height(), 0, 1, mu);
  stacked.set_submatrix(mu.height(), stacked.height(), 0, 1, alpha);

  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat kern = kernel(all_dofs, all_dofs);

  ie_Mat dense(all_dofs.size() + boundary->holes.size(),
               all_dofs.size() + boundary->holes.size());
  ie_Mat ident(boundary->holes.size(), boundary->holes.size());
  ident.eye(boundary->holes.size());
  dense.set_submatrix(0, all_dofs.size(), 0, all_dofs.size(), kern);
  dense.set_submatrix(0, all_dofs.size(), all_dofs.size(), dense.width(), U);
  dense.set_submatrix(all_dofs.size(), dense.height(), 0, all_dofs.size(), Psi);
  dense.set_submatrix(all_dofs.size(), dense.height(),
                      all_dofs.size(), dense.width(), -ident);

  ie_Mat fzero_prime(mu.height() + alpha.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, stacked, 0., &fzero_prime);

  ie_Mat err1 = (fzero_prime(0, mu.height(), 0, 1) - boundary->boundary_values);
  ie_Mat err2 = (fzero_prime(mu.height(), fzero_prime.height(), 0, 1));
  double max = 0;
  for (int i = 0; i < err1.height(); i++) {
    max = std::max(max, std::abs(err1.get(i, 0)));
  }
  for (int i = 0; i < err2.height(); i++) {
    max = std::max(max, std::abs(err2.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 10000);
  EXPECT_LE(max, 100 * config.id_tol);
}


// TEST(IeSolverTest, LaplaceNeumannAnnulusBackwardError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
//   config.boundary_condition = BoundaryCondition::ALL_ONES;
//   config.boundary_shape = Boundary::BoundaryShape::ANNULUS;
//   std::unique_ptr<Boundary> boundary =
//     std::unique_ptr<Boundary>(new Annulus());
//   boundary->initialize(config.num_boundary_points, config.boundary_condition);

//   QuadTree quadtree;
//   quadtree.initialize_tree(boundary.get(), std::vector<double>(),
//                            config.solution_dimension, config.domain_dimension);

//   SkelFactorization skel_factorization(config.id_tol,
//                                        config.is_strong_admissibility,
//                                        config.solution_dimension,
//                                        config.domain_dimension);

//   Kernel kernel(config.solution_dimension, config.domain_dimension,
//                 config.pde, boundary.get(), std::vector<double>());
//   skel_factorization.skeletonize(kernel, &quadtree);

//   ie_Mat mu;
//   linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);


//   std::vector<unsigned int> all_dofs;
//   for (int i = 0; i < boundary->weights.size(); i++) {
//     all_dofs.push_back(i);
//   }
//   ie_Mat dense = kernel(all_dofs, all_dofs);
//   ie_Mat f_prime(mu.height(), 1);
//   ie_Mat::gemm(NORMAL, NORMAL, 1., dense, mu, 0., &f_prime);
//   ie_Mat err = (f_prime - boundary->boundary_values);
//   double max = 0;
//   for (int i = 0; i < err.height(); i++) {
//     max = std::max(max, std::abs(err.get(i, 0)));
//   }
//   EXPECT_LE(dense.condition_number(), 100000);
//   EXPECT_LE(max, 100 * config.id_tol);
// }


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

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  SkelFactorization skel_factorization(config.id_tol,
                                       config.is_strong_admissibility,
                                       config.solution_dimension,
                                       config.domain_dimension);

  Kernel kernel(config.solution_dimension, config.domain_dimension,
                config.pde, boundary.get(), std::vector<double>());

  ie_Mat U = initialize_U_mat(config.pde, boundary->holes, boundary->points);
  ie_Mat Psi = initialize_Psi_mat(config.pde, boundary->holes,
                                  boundary.get());
  skel_factorization.U = U;
  skel_factorization.Psi = Psi;
  skel_factorization.skeletonize(kernel, &quadtree);

  ie_Mat mu, alpha;
  linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu,
               &alpha);
  ie_Mat stacked(mu.height() + alpha.height(), 1);
  stacked.set_submatrix(0, mu.height(), 0, 1, mu);
  stacked.set_submatrix(mu.height(), stacked.height(), 0, 1, alpha);

  std::vector<unsigned int> all_dofs;
  for (int i = 0; i < 2 * boundary->weights.size(); i++) {
    all_dofs.push_back(i);
  }
  ie_Mat kern = kernel(all_dofs, all_dofs);

  ie_Mat dense(all_dofs.size() + 3 * boundary->holes.size(),
               all_dofs.size() + 3 * boundary->holes.size());
  ie_Mat ident(3 * boundary->holes.size(), 3 * boundary->holes.size());
  ident.eye(3 * boundary->holes.size());
  dense.set_submatrix(0, all_dofs.size(), 0, all_dofs.size(), kern);
  dense.set_submatrix(0, all_dofs.size(), all_dofs.size(), dense.width(), U);
  dense.set_submatrix(all_dofs.size(), dense.height(), 0, all_dofs.size(), Psi);
  dense.set_submatrix(all_dofs.size(), dense.height(),
                      all_dofs.size(), dense.width(), -ident);


  // ie_Mat cond(dense.height(), dense.width());
  // cond.eye(dense.height());
  // dense += cond;
  // dense += cond;
  // dense += cond;

  ie_Mat fzero_prime(mu.height() + alpha.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., dense, stacked, 0., &fzero_prime);

  ie_Mat err1 = (fzero_prime(0, mu.height(), 0, 1) - boundary->boundary_values);
  ie_Mat err2 = (fzero_prime(mu.height(), fzero_prime.height(), 0, 1));
  double max = 0;
  for (int i = 0; i < err1.height(); i++) {
    max = std::max(max, std::abs(err1.get(i, 0)));
  }
  for (int i = 0; i < err2.height(); i++) {
    max = std::max(max, std::abs(err2.get(i, 0)));
  }
  EXPECT_LE(dense.condition_number(), 10000);
  EXPECT_LE(max, 100 * config.id_tol);
}

// TEST(IeSolverTest, LaplaceCircleOnesSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new Circle());
//   config.boundary_condition = BoundaryCondition::ALL_ONES;
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-4);
// }


// TEST(IeSolverTest, LaplaceCircleElectronTinyError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new Circle());
//   config.id_tol = 1e-10;
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }


// TEST(IeSolverTest, StokesCircleSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.pde = ie_solver_config::STOKES;
//   config.boundary_condition =  BoundaryCondition::STOKES;

//   config.boundary.reset(new Circle());
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }


// TEST(IeSolverTest, LaplaceRoundedSquareElectronSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new RoundedSquare());
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-4);
// }


// TEST(IeSolverTest, LaplaceRoundedSquareWithBumpElectronSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new RoundedSquareWithBump());
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-4);
// }



// TEST(IeSolverTest, LaplaceSquigglyElectronSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new Squiggly());
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-4);
// }


// TEST(IeSolverTest, LaplaceCubicSplineElectronSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new CubicSpline());
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-4);
// }


// TEST(IeSolverTest, LaplaceEllipsesElectronSmallError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new Ellipses());
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-4);
// }


// The below tests await a fast quadrature rule
// TEST(IeSolverTest, LaplaceRoundedSquareElectronTinyError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new RoundedSquare());
//   config.testing = true;
//   config.id_tol = 1e-10;
//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }


// TEST(IeSolverTest, LaplaceRoundedSquareWithBumpElectronTinyError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new RoundedSquareWithBump());
//   config.testing = true;
//   config.id_tol = 1e-10;

//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }



// TEST(IeSolverTest, LaplaceSquigglyElectronTinyError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new Squiggly());
//   config.testing = true;
//   config.id_tol = 1e-10;

//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }


// TEST(IeSolverTest, LaplaceCubicSplineElectronTinyError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new CubicSpline());
//   config.testing = true;
//   config.id_tol = 1e-10;

//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }


// TEST(IeSolverTest, LaplaceEllipsesElectronTinyError) {
//   srand(0);
//   ie_solver_config config = ie_solver_config();
//   config.boundary.reset(new Ellipses());
//   config.testing = true;
//   config.id_tol = 1e-10;

//   // First we init the boundary so we can correct the num_boundary_points
//   config.boundary->initialize(config.num_boundary_points,
//                               config.boundary_condition);
//   config.num_boundary_points = config.boundary->weights.size();
//   double error = boundary_integral_solve(config, nullptr);
//   EXPECT_GE(error, 0);
//   EXPECT_LE(error, 1e-8);
// }

}  // namespace ie_solver
