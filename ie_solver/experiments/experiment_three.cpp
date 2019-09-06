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


ie_solver_config get_experiment_three_config(){
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = 3000;
  config.domain_size = 150;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::STOKES;
  config.boundary_shape = Boundary::BoundaryShape::EX3;
  return config;
}


void run_experiment3() {
  ie_solver_config config = get_experiment_three_config();
  int domain_dimension = 2;
  int solution_dimension = 2;

  std::unique_ptr<Boundary> boundary, perturbed_boundary;
  boundary.reset(new Ex3Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::STOKES);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           solution_dimension, domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);

  ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);

  io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
                         domain_points, config.solution_dimension);
  io::write_boundary_to_file("output/data/ie_solver_boundary.txt", boundary->points);
  io::write_quadtree_to_file("output/data/ie_solver_tree.txt", quadtree);
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment3();
  return 0;
}

