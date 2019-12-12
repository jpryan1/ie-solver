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
#include "ie_solver/boundaries/ex1boundary.h"

namespace ie_solver {

ie_solver_config get_experiment_one_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = pow(2, 14);
  config.domain_size = 200;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX1;
  return config;
}


void run_experiment1() {
  // double start = omp_get_wtime();
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
                    quadtree.max);

  // // We'll iteratively reinitialized another Boundary and use that
  // // to update the quadtree's Boundary.
  // std::unique_ptr<Boundary> perturbed_boundary =
  //   std::unique_ptr<Boundary>(new Ex1Boundary());

  // perturbed_boundary->initialize(config.num_boundary_points,
  //                                config.boundary_condition);

  // int FRAME_CAP = 30;
  // for (int frame = 0; frame < FRAME_CAP; frame++) {
  //   double ang = (frame / (0.0 + FRAME_CAP)) * 2 * M_PI;

  //   perturbed_boundary->perturbation_parameters[0] = ang;
  //   perturbed_boundary->initialize(config.num_boundary_points,
  //                                  config.boundary_condition);

  //   double pstr = omp_get_wtime();
  //   quadtree.perturb(*perturbed_boundary.get());
  //   double pend = omp_get_wtime();
  //   std::cout << "timing: qtree_perturb " << (pend - pstr) << std::endl;

  //   ie_Mat solution = boundary_integral_solve(config, &quadtree,
  //                     domain_points);
  //   io::write_solution_to_file("output/bake/sol/" + std::to_string(frame)
  //                              + ".txt", solution, domain_points,
  //                              config.solution_dimension);
  //   io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
  //                                frame)  + ".txt",
  //                              perturbed_boundary->points);
  //   io::write_quadtree_to_file("output/bake/tree/" + std::to_string(
  //                                frame)  + ".txt", quadtree);
  // }

  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);
  io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
                             domain_points,
                             config.solution_dimension);
  io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
                             boundary->points);
  io::write_quadtree_to_file("output/bake/tree/ie_solver_tree.txt", quadtree);
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);
  ie_solver::run_experiment1();
  return 0;
}

