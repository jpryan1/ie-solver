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
#include "ie_solver/boundaries/ex2boundary.h"
#include "ie_solver/boundaries/ex3boundary.h"

namespace ie_solver {

ie_solver_config get_experiment_two_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.num_boundary_points = 2000;
  config.domain_size = 100;
  config.domain_dimension = 2;
  config.solution_dimension = 1;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX2;
  return config;
}

void run_experiment2() {
  double start = omp_get_wtime();
  ie_solver_config config = get_experiment_two_config();
  int domain_dimension = 2;
  int solution_dimension = 1;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());

  double init_start = omp_get_wtime();
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);
  double init_end = omp_get_wtime();
  std::cout << "timing: boundary_init " << (init_end - init_start) << std::endl;

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           solution_dimension, domain_dimension);

  std::vector<double> domain_points;
  get_domain_points(config.domain_size, &domain_points, quadtree.min,
                    quadtree.max);


  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);
  for (int frame = 0; frame < 30; frame++) {
    double ang = (frame / 30.0) * 2 * M_PI;

    perturbed_boundary->perturbation_parameters[0] = ang;
    perturbed_boundary->perturbation_parameters[1] = ang + M_PI;

    perturbed_boundary->initialize(config.num_boundary_points,
                                   config.boundary_condition);
    quadtree.perturb(*perturbed_boundary.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree,
                      domain_points);

    // Make solution zero mean.
    double avg = 0.;
    int total = 0;
    for (int i = 0; i < solution.height(); i++) {
      if (solution.get(i, 0) != 0.0) {
        avg += solution.get(i, 0);
        total++;
      }
    }
    avg /= total;
    for (int i = 0; i < solution.height(); i++) {
      if (solution.get(i, 0) != 0.0) {
        solution.set(i, 0, solution.get(i, 0) - avg);
      }
    }
    io::write_solution_to_file("output/bake/sol/" + std::to_string(frame)
                               + ".txt", solution, domain_points,
                               config.solution_dimension);
    io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                 frame)  + ".txt",
                               perturbed_boundary->points);
    io::write_quadtree_to_file("output/bake/tree/ie_solver_tree.txt",
                               quadtree);
  }

  // ie_Mat solution = boundary_integral_solve(config, &quadtree,
  //                   domain_points);

  // io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
  //                            domain_points, config.solution_dimension);
  // io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
  //                            boundary->points);
  // io::write_quadtree_to_file("output/data/ie_solver_tree.txt", quadtree);

  // double end = omp_get_wtime();
  // std::cout << "timing: experiment_two_total " << (end - start) << std::endl;
}

}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment2();
  return 0;
}

