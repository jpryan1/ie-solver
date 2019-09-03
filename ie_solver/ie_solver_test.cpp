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
#include "ie_solver/skel_factorization/skel_factorization.h"
#include "ie_solver/quadtree.h"
#include "ie_solver/kernel.h"
#include "ie_solver/log.h"
#include "ie_solver/linear_solve.h"
#include "gtest/gtest.h"

namespace ie_solver {


TEST(IeSolverTest, LaplaceCircleElectronSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new Circle());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}


TEST(IeSolverTest, LaplaceCircleOnesSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new Circle());
  config.testing = true;
  config.boundary_condition = BoundaryCondition::ALL_ONES;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}


TEST(IeSolverTest, LaplaceCircleElectronTinyError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new Circle());
  config.testing = true;
  config.id_tol = 1e-10;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-8);
}


TEST(IeSolverTest, StokesCircleSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.pde = ie_solver_config::STOKES;
  config.boundary_condition =  BoundaryCondition::STOKES;

  config.boundary.reset(new Circle());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-8);
}


TEST(IeSolverTest, LaplaceRoundedSquareElectronSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new RoundedSquare());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}


TEST(IeSolverTest, LaplaceRoundedSquareWithBumpElectronSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new RoundedSquareWithBump());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}



TEST(IeSolverTest, LaplaceSquigglyElectronSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new Squiggly());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}


TEST(IeSolverTest, LaplaceCubicSplineElectronSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new CubicSpline());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}


TEST(IeSolverTest, LaplaceEllipsesElectronSmallError) {
  srand(0);
  ie_solver_config config = ie_solver_config();
  config.boundary.reset(new Ellipses());
  config.testing = true;
  // First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();
  double error = boundary_integral_solve(config, nullptr);
  EXPECT_GE(error, 0);
  EXPECT_LE(error, 1e-4);
}


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
