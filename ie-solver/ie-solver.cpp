// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <string.h>
#include <fstream>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include "ie-solver/ie_mat.h"
#include "ie-solver/initialization.h"
#include "ie-solver/tools/ie_solver_tools.h"
#include "ie-solver/quadtree.h"
#include "ie-solver/kernel.h"
#include "ie-solver/log.h"
#include "ie-solver/helpers.h"
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/boundaries/circle.h"
#include "ie-solver/boundaries/rounded_square.h"
#include "ie-solver/boundaries/rounded_square_with_bump.h"
#include "ie-solver/boundaries/squiggly.h"
#include "ie-solver/boundaries/annulus.h"
#include "ie-solver/boundaries/cubic_spline.h"

namespace ie_solver {


void run_animation(const ie_solver_config& config) {
  std::unique_ptr<Boundary> boundary;
  switch (config.boundary_shape) {
    case Boundary::BoundaryShape::CIRCLE:
      boundary.reset(new Circle());
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE:
      boundary.reset(new RoundedSquare());
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE_WITH_BUMP:
      boundary.reset(new RoundedSquareWithBump());
      break;
    case Boundary::BoundaryShape::SQUIGGLY:
      boundary.reset(new Squiggly());
      break;
    case Boundary::BoundaryShape::ANNULUS:
      boundary.reset(new Annulus());
      break;
    case Boundary::BoundaryShape::CUBIC_SPLINE:
      boundary.reset(new CubicSpline());
      break;
  }

  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);
  for (int frame = 0; frame < 15; frame++) {
    double x = 0.4 + 0.2 * (frame / 15.0);
    double y = 0.5;

    boundary->holes[0].center = Vec2(x, y);
    boundary->initialize(500, config.boundary_condition);
    quadtree.reset(boundary.get());
    ie_Mat solution = boundary_integral_solve(config, boundary.get(), &quadtree,
                      domain_points);
    std::string filename = "output/bake/sol/" + std::to_string(frame)  + ".txt";
    write_solution_to_file(filename, solution, domain_points,
                           config.solution_dimension);
  }
}


void run_time_trial(const ie_solver_config & config) {
  // std::vector<double> n_times, eps_times;
  // int scale_n[] = {5000, 6000, 7000, 8000, 9000};
  // double scale_eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};

  // for (int n : scale_n) {
  //   config.num_boundary_points = n;
  //   config.id_tol = DEFAULT_ID_TOL;
  //   ie_solver::boundary_integral_solve(config, &n_times);
  // }

  // for (double eps : scale_eps) {
  //   config.num_boundary_points = DEFAULT_NUM_DISCRETIZATION_POINTS * 10;
  //   config.id_tol = eps;
  //   ie_solver::boundary_integral_solve(config, &eps_times);
  // }

  // ie_solver::write_times_to_files(scale_n, n_times, scale_eps, eps_times);
}


void run_single_solve(const ie_solver_config & config) {
  std::unique_ptr<Boundary> boundary;
  switch (config.boundary_shape) {
    case Boundary::BoundaryShape::CIRCLE:
      boundary.reset(new Circle());
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE:
      boundary.reset(new RoundedSquare());
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE_WITH_BUMP:
      boundary.reset(new RoundedSquareWithBump());
      break;
    case Boundary::BoundaryShape::SQUIGGLY:
      boundary.reset(new Squiggly());
      break;
    case Boundary::BoundaryShape::ANNULUS:
      boundary.reset(new Annulus());
      break;
    case Boundary::BoundaryShape::CUBIC_SPLINE:
      boundary.reset(new CubicSpline());
      break;
  }
  boundary->boundary_shape = config.boundary_shape;
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);

  ie_Mat solution = boundary_integral_solve(config, boundary.get(), &quadtree,
                    domain_points);
  double error = stokes_error(solution, config.id_tol, domain_points,
                              boundary.get());
  std::cout << "Error: " << error << std::endl;
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  // TODO(John) allow for command line args for setting parameters
  srand(0);  // omp_get_wtime());

  ie_solver::ie_solver_config config;
  if (ie_solver::parse_input_into_config(argc, argv, &config)) {
    return 1;
  }

// First we init the boundary so we can correct the num_boundary_points

  if (config.animation) {
    ie_solver::run_animation(config);
  } else if (config.scaling) {
    ie_solver::run_time_trial(config);

  } else {
    ie_solver::run_single_solve(config);
  }
  return 0;
}

