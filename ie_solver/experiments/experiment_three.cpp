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
#include "ie_solver/boundaries/ex3boundary.h"

namespace ie_solver {


ie_solver_config get_experiment_three_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = pow(2, 13);
  config.domain_size = 70;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX3;
  return config;
}


void run_experiment3() {
  ie_solver_config config = get_experiment_three_config();

  std::unique_ptr<Boundary> boundary;
  boundary.reset(new Ex3Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);

  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);

  int FRAME_CAP = 10;
  for (int frame = 0; frame < FRAME_CAP; frame++) {
    double ang = (frame / (FRAME_CAP + 0.)) * 2 * M_PI;

    perturbed_boundary->perturbation_parameters[0] = ang;
    perturbed_boundary->initialize(config.num_boundary_points,
                                   config.boundary_condition);

    quadtree.perturb(*perturbed_boundary.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree,
                      domain_points);
    io::write_solution_to_file("output/bake/sol/" + std::to_string(frame)
                               + ".txt", solution, domain_points,
                               config.solution_dimension);
    io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                 frame)  + ".txt",
                               perturbed_boundary->points);
    io::write_quadtree_to_file("output/bake/tree/" + std::to_string(
                                 frame)  + ".txt", quadtree);
  }
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment3();
  return 0;
}

