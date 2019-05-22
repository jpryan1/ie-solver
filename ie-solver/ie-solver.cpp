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
#include "ie-solver/boundaries/ex1boundary.h"
#include "ie-solver/boundaries/ex3boundary.h"

namespace ie_solver {


void run_experiment1() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = 1000;
  config.domain_size = 49;
  int domain_dimension = 2;
  int solution_dimension = 2;
  config.domain_dimension = 2;
  config.solution_dimension = 2;

  config.boundary_condition = Boundary::BoundaryCondition::STOKES;
  config.boundary_shape = Boundary::BoundaryShape::EX1;



  std::unique_ptr<Boundary> boundary, perturbed_boundary;
  boundary.reset(new Ex1Boundary());
  perturbed_boundary.reset(new Ex1Boundary());
  boundary->initialize(config.num_boundary_points,
                       Boundary::BoundaryCondition::STOKES);
  perturbed_boundary->initialize(config.num_boundary_points,
                                 Boundary::BoundaryCondition::STOKES);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           solution_dimension, domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);


  for (int frame = 0; frame < 30; frame++) {
    double ang = (frame / 60.0) * 2 * M_PI;

    perturbed_boundary->holes[0].center = Vec2(0.5 + 0.3 * cos(M_PI + ang),
                                          0.5 + 0.3 * sin(M_PI + ang));
    perturbed_boundary->holes[3].center = Vec2(0.5 + 0.3 * cos(ang),
                                          0.5 + 0.3 * sin(ang));

    perturbed_boundary->initialize(config.num_boundary_points,
                                   config.boundary_condition);

    quadtree.perturb(*perturbed_boundary.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);
    std::string filename = "output/bake/sol/" + std::to_string(frame)  + ".txt";
    write_solution_to_file(filename, solution, domain_points,
                           config.solution_dimension);
  }
  // ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);

  // write_solution_to_file("output/data/ie_solver_solution.txt", solution,
  //                        domain_points, config.solution_dimension);
  // write_boundary_to_file(boundary->points);
  // quadtree.write_quadtree_to_file();
}


void run_experiment2() {
}


void run_experiment3() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = 1000;
  config.domain_size = 49;
  int domain_dimension = 2;
  int solution_dimension = 2;
  config.domain_dimension = 2;
  config.solution_dimension = 2;

  config.boundary_condition = Boundary::BoundaryCondition::STOKES;
  config.boundary_shape = Boundary::BoundaryShape::EX3;



  std::unique_ptr<Boundary> boundary, perturbed_boundary;
  boundary.reset(new Ex3Boundary());
  boundary->initialize(config.num_boundary_points,
                       Boundary::BoundaryCondition::STOKES);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           solution_dimension, domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);

  ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);

  write_solution_to_file("output/data/ie_solver_solution.txt", solution,
                         domain_points, config.solution_dimension);
  write_boundary_to_file(boundary->points);
  quadtree.write_quadtree_to_file();
}


void run_animation(const ie_solver_config& config) {
  std::unique_ptr<Boundary> boundary, perturbed_boundary;
  switch (config.boundary_shape) {
    case Boundary::BoundaryShape::CIRCLE:
      boundary.reset(new Circle());
      perturbed_boundary.reset(new Circle());
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE:
      boundary.reset(new RoundedSquare());
      perturbed_boundary.reset(new RoundedSquare());
      break;
    case Boundary::BoundaryShape::ROUNDED_SQUARE_WITH_BUMP:
      boundary.reset(new RoundedSquareWithBump());
      perturbed_boundary.reset(new RoundedSquareWithBump());
      break;
    case Boundary::BoundaryShape::SQUIGGLY:
      boundary.reset(new Squiggly());
      perturbed_boundary.reset(new Squiggly());
      break;
    case Boundary::BoundaryShape::ANNULUS:
      boundary.reset(new Annulus());
      perturbed_boundary.reset(new Annulus());
      break;
    case Boundary::BoundaryShape::CUBIC_SPLINE:
      boundary.reset(new CubicSpline());
      perturbed_boundary.reset(new CubicSpline());
      break;
  }

  // Currently the animation is of moving holes inside an annulus
  assert(config.boundary_shape == Boundary::BoundaryShape::ANNULUS);

  boundary->initialize(1000, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);
  for (int frame = 0; frame < 20; frame++) {
    double ang = (frame / 60.0) * 2 * M_PI;
    perturbed_boundary->holes.clear();
    Hole hole;
    hole.center = Vec2(0.5 + 0.1 * cos(ang), 0.5 + 0.1 * sin(ang));
    hole.radius = 0.025;
    perturbed_boundary->holes.push_back(hole);
    hole.center = Vec2(0.5 + 0.1 * cos(M_PI + ang),
                       0.5 + 0.1 * sin(M_PI + ang));
    perturbed_boundary->holes.push_back(hole);
    perturbed_boundary->initialize(1000, config.boundary_condition);
    quadtree.perturb(*perturbed_boundary.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);
    std::string filename = "output/bake/sol/" + std::to_string(frame)  + ".txt";
    write_solution_to_file(filename, solution, domain_points,
                           config.solution_dimension);
  }
}


void run_time_trial(const ie_solver_config & config) {
  std::unique_ptr<Boundary> boundary;
  boundary.reset(new CubicSpline());
  boundary->boundary_shape = config.boundary_shape;
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);
  double avg_skel_time = 0;
  double avg_solve_time = 0;
  bie_time_trial(config, &quadtree, &avg_skel_time, &avg_solve_time);
  std::cout << "Avg skel time: " << avg_skel_time << " Avg solve time: " <<
            avg_solve_time << std::endl;
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

  ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);

  write_solution_to_file("output/data/ie_solver_solution.txt", solution,
                         domain_points, config.solution_dimension);
  write_boundary_to_file(boundary->points);
  quadtree.write_quadtree_to_file();

  double error;
  switch (config.pde) {
    case ie_solver_config::Pde::LAPLACE:
      error = laplace_error(solution, config.id_tol, domain_points,
                            boundary.get());
      break;
    case ie_solver_config::Pde::STOKES:
      error = stokes_error(solution, config.id_tol, domain_points,
                           boundary.get());
      break;
  }
  std::cout << "Error: " << error << std::endl;
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  // TODO(John) allow for command line args for setting parameters
  srand(0);  // omp_get_wtime());

  ie_solver::run_experiment3();


//   ie_solver::ie_solver_config config;



//   if (ie_solver::parse_input_into_config(argc, argv, &config)) {
//     return 1;
//   }

// // First we init the boundary so we can correct the num_boundary_points

//   if (config.animation) {
//     ie_solver::run_animation(config);
//   } else if (config.scaling) {
//     ie_solver::run_time_trial(config);

//   } else {
//     ie_solver::run_single_solve(config);
//   }
  return 0;
}

