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
#include "ie_solver/helpers.h"
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
    io::write_solution_to_file(filename, solution, domain_points,
                           config.solution_dimension);
  }
  // ie_Mat solution = boundary_integral_solve(config, &quadtree, domain_points);

  // io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
  //                        domain_points, config.solution_dimension);
  // io::write_boundary_to_file(boundary->points);
  //   io::write_quadtree_to_file(quadtree);

}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment1();
  return 0;
}

