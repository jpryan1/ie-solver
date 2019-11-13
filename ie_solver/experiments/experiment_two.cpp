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
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.num_boundary_points = pow(2, 13);
  config.domain_size = 50;
  config.solution_dimension = 1;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX2;
  return config;
}


void run_experiment2() {
  // double start = omp_get_wtime();
  ie_solver_config config = get_experiment_two_config();

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  std::vector<double> domain_points;
  // get_domain_points(config.domain_size, &domain_points, quadtree.min,
  //                   quadtree.max);
  domain_points.push_back(0.4999);
  domain_points.push_back(0.5001);

  domain_points.push_back(0.5001);
  domain_points.push_back(0.5001);
  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);

  double current_ang1 = 0.2;
  double current_ang2 = 2;
  double alpha = 0.1;
  double h = 1e-4;
  int FRAME_CAP = 10;
  for (int frame1 = 0; frame1 < FRAME_CAP; frame1++) {
    double findiff1[4];
    double samples1[4] = {current_ang1 - 2 * h, current_ang1 - h,
                          current_ang1 + h, current_ang1 + 2 * h
                         };
    for (int i = 0; i < 4; i++) {
      double temp_ang = samples1[i];
      perturbed_boundary->perturbation_parameters[0] = temp_ang;
      perturbed_boundary->initialize(config.num_boundary_points,
                                     config.boundary_condition);
      quadtree.perturb(*perturbed_boundary.get());
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);

      int sh = solution.height();
      double gradient = (solution.get(sh - 1, 0) - solution.get(sh - 2, 0))
                        / 0.0002;
      findiff1[i] = gradient;
      if (i == 2) std::cout << gradient << std::endl;
    }
    perturbed_boundary->perturbation_parameters[0] = current_ang1;

    double findiff2[4];
    double samples2[4] = {current_ang2 - 2 * h, current_ang2 - h,
                          current_ang2 + h, current_ang2 + 2 * h
                         };
    for (int i = 0; i < 4; i++) {
      double temp_ang = samples2[i];
      perturbed_boundary->perturbation_parameters[2] = temp_ang;
      perturbed_boundary->initialize(config.num_boundary_points,
                                     config.boundary_condition);
      quadtree.perturb(*perturbed_boundary.get());
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);

      int sh = solution.height();
      double gradient = (solution.get(sh - 1, 0) - solution.get(sh - 2, 0))
                        / (0.0002);
      findiff2[i] = gradient;
      if (i == 2) {
        io::write_solution_to_file("output/bake/sol/" + std::to_string(
                                     frame1)  + ".txt", solution, domain_points,
                                   config.solution_dimension);
        io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                     frame1) + ".txt", perturbed_boundary->points);
        io::write_quadtree_to_file("output/bake/tree/" + std::to_string(
                                     frame1)  + ".txt", quadtree);
      }
    }

    double grad1 = (findiff1[0] - 8 * findiff1[1] + 8 * findiff1[2] - findiff1[3]) /
                   (12 * h);
    double grad2 = (findiff2[0] - 8 * findiff2[1] + 8 * findiff2[2] - findiff2[3]) /
                   (12 * h);

    current_ang1 += alpha * grad1;
    current_ang2 += alpha * grad2;
    std::cout << grad1 << " " << grad2 << std::endl;
    perturbed_boundary->perturbation_parameters[0] = current_ang1;
    perturbed_boundary->perturbation_parameters[1] = current_ang2;


  }
}

}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment2();
  return 0;
}

