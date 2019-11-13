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
  config.num_boundary_points = pow(2, 12);
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
  get_domain_points(20, &domain_points, quadtree.min,
                    quadtree.max);


  domain_points.push_back(0.6);
  domain_points.push_back(0.325);

  domain_points.push_back(0.7);
  domain_points.push_back(0.325);

  domain_points.push_back(0.8);
  domain_points.push_back(0.325);

  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);

  int FRAME_CAP = 15;
  double current_ang = 0.2;
  double alpha = 0.1;
  double h = 1e-4;
  for (int frame = 0; frame < FRAME_CAP; frame++) {
    // Estimate gradient
    double findiff[4];
    double samples[4] = {current_ang - 2 * h, current_ang - h, current_ang + h,
                         current_ang + 2 * h
                        };
    // #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 4; i++) {
      double temp_ang = samples[i];
      perturbed_boundary->perturbation_parameters[0] = temp_ang;
      perturbed_boundary->initialize(config.num_boundary_points,
                                     config.boundary_condition);
      quadtree.perturb(*perturbed_boundary.get());
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);
      double flow = 0.;
      for (int i = solution.height() - 12; i < solution.height(); i += 2) {
        flow += solution.get(i, 0);
      }
      findiff[i] = flow;
      if (i == 2) {
        io::write_solution_to_file("output/bake/sol/" + std::to_string(
                                     frame)  + ".txt", solution, domain_points,
                                   config.solution_dimension);
        io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                     frame) + ".txt", perturbed_boundary->points);
        io::write_quadtree_to_file("output/bake/tree/" + std::to_string(
                                     frame)  + ".txt", quadtree);
      }
    }
    double grad = (findiff[0] - 8 * findiff[1] + 8 * findiff[2] - findiff[3]) /
                  (12 * h);
    std::cout << "Current flow and ang " << findiff[0] << " "
              << (current_ang - h) / M_PI  << "pi, grad " << grad <<
              std::endl;
    current_ang -= alpha * grad;
  }
}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment3();
  return 0;
}

