// Copyright 2019 John Paul Ryan
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
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/kernel/kernel.h"
#include "ie_solver/log.h"
#include "ie_solver/linear_solve.h"
#include "ie_solver/boundaries/ex2boundary.h"

namespace ie_solver {


ie_solver_config get_experiment_two_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = pow(2, 12);
  config.domain_size = 50;  // 375 for image
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX2;
  return config;
}


void run_experiment2() {
  ie_solver_config config = get_experiment_two_config();
  std::unique_ptr<Boundary> boundary
    =  std::unique_ptr<Boundary>(new Ex2Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  double x_min = boundary->points[0], x_max = boundary->points[0],
         y_min = boundary->points[1], y_max = boundary->points[1];
  for (int i = 0; i < boundary->points.size(); i += 2) {
    x_min = std::min(x_min, boundary->points[i]);
    x_max = std::max(x_max, boundary->points[i]);
    y_min = std::min(y_min, boundary->points[i + 1]);
    y_max = std::max(y_max, boundary->points[i + 1]);
  }
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, x_min, x_max, y_min,
                    y_max);

  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);

  int FRAME_CAP = 20;
  for (int frame = 0; frame < FRAME_CAP; frame++) {
    int rand_idx = floor(8 * (rand() / (0. + RAND_MAX)));
    boundary->perturbation_parameters[rand_idx] = 0.35
        + 0.3 * (rand() / (0. + RAND_MAX));
    boundary->initialize(config.num_boundary_points,
                         config.boundary_condition);
    quadtree.perturb(*boundary.get());
    solution = boundary_integral_solve(config, &quadtree,
                                       domain_points);
    io::write_solution_to_file("output/bake/sol/" + std::to_string(
                                 frame)  + ".txt", solution, domain_points,
                               config.solution_dimension);
    io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                 frame) + ".txt", boundary->points);
  }

  // io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
  //                            domain_points,
  //                            config.solution_dimension);
  // io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
  //                            boundary->points);
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment2();
  return 0;
}

