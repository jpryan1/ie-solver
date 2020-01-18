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
  config.domain_size = 10;  // 200;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX1;
  return config;
}


void run_experiment1(int N) {
  // double start = omp_get_wtime();
  ie_solver_config config = get_experiment_one_config();
  config.num_boundary_points = N;
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
  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);

  boundary->perturbation_parameters[0] = (M_PI / 2.) + 0.1;
  boundary->initialize(config.num_boundary_points, config.boundary_condition);

  quadtree.perturb(*boundary.get());
  solution = boundary_integral_solve(config, &quadtree, domain_points);
}

}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);
  std::cout << "\nEXPERIMENT ONE" << std::endl;
  for (int i = 15; i < 20 ; i++) {
    std::cout << "i = " << i << std::endl;
    for (int k = 0; k < 3; k++) {
      ie_solver::run_experiment1(pow(2, i));
    }
  }
  return 0;
}

