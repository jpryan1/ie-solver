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

ie_solver_config get_experiment_one_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = 3000;
  config.domain_size = 49;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::STOKES;
  config.boundary_shape = Boundary::BoundaryShape::EX1;
  return config;
}


void run_experiment1() {

  double start = omp_get_wtime();
  ie_solver_config config = get_experiment_one_config();
  int domain_dimension = 2;
  int solution_dimension = 2;

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex1Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::STOKES);

  QuadTree quadtree;

  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           solution_dimension, domain_dimension);

  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);

  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex1Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);
  // for (int frame = 0; frame < 30; frame++) {
  //   double ang = (frame / 120.0) * 2 * M_PI;

  //   perturbed_boundary->holes[0].center = Vec2(0.5 + 0.3 * cos(M_PI + ang),
  //                                         0.5 + 0.3 * sin(M_PI + ang));
  //   perturbed_boundary->holes[3].center = Vec2(0.5 + 0.3 * cos(ang),
  //                                         0.5 + 0.3 * sin(ang));

  //   perturbed_boundary->initialize(config.num_boundary_points,
  //                                 config.boundary_condition);
  //   quadtree.perturb(*perturbed_boundary.get());
  //   ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);
  //   io::write_solution_to_file("output/bake/sol/" + std::to_string(frame)  + ".txt",
  //                             solution, domain_points, config.solution_dimension);
  //   io::write_boundary_to_file("output/bake/boundary/" + std::to_string(frame)  + ".txt",
  //                             perturbed_boundary->points);
  //   io::write_quadtree_to_file("output/bake/tree/ie_solver_tree.txt", quadtree);

  // }
  ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);

  io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
                             domain_points, config.solution_dimension);
  io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
                             boundary->points);
  io::write_quadtree_to_file("output/data/ie_solver_tree.txt", quadtree);

  double end = omp_get_wtime();
  std::cout << "timing: experiment_one_total " << (end - start) << std::endl;

}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment1();
  return 0;
}

