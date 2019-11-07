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
#include "ie_solver/boundaries/ex2boundary.h"
#include "ie_solver/boundaries/ex3boundary.h"

namespace ie_solver {

ie_solver_config get_landscape_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.num_boundary_points = pow(2, 13);
  config.domain_size = 100;
  config.domain_dimension = 2;
  config.solution_dimension = 1;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX2;
  return config;
}

void print_landscapes() {
  ie_solver_config config = get_landscape_config();
  std::unique_ptr<Boundary> boundary2 =
    std::unique_ptr<Boundary>(new Ex2Boundary());
  boundary2->initialize(config.num_boundary_points,
                        BoundaryCondition::DEFAULT);
  QuadTree quadtree2;
  quadtree2.initialize_tree(boundary2.get(), std::vector<double>(),
                            config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points2;

  domain_points2.push_back(0.4999);
  domain_points2.push_back(0.4999);

  domain_points2.push_back(0.4999);
  domain_points2.push_back(0.5001);

  domain_points2.push_back(0.5001);
  domain_points2.push_back(0.4999);

  domain_points2.push_back(0.5001);
  domain_points2.push_back(0.5001);

  std::unique_ptr<Boundary> perturbed_boundary2 =
    std::unique_ptr<Boundary>(new Ex2Boundary());
  perturbed_boundary2->initialize(config.num_boundary_points,
                                  config.boundary_condition);

  int FRAME_CAP = 50;
  std::vector<double> angs_and_grads;
  for (int frame1 = 0; frame1 < FRAME_CAP; frame1++) {
    double ang1 = (frame1 / (0.0 + FRAME_CAP)) * 2 * M_PI;

    for (int frame2 = 0; frame2 < FRAME_CAP; frame2++) {
      double ang2 = ((frame2 / (0.0 + FRAME_CAP)) * (M_PI)) +
                    (M_PI / 2.) + ang1;
      perturbed_boundary2->perturbation_parameters[0] = ang1;
      perturbed_boundary2->perturbation_parameters[1] = ang2;

      perturbed_boundary2->initialize(config.num_boundary_points,
                                      config.boundary_condition);
      quadtree2.perturb(*perturbed_boundary2.get());
      ie_Mat solution = boundary_integral_solve(config, &quadtree2,
                        domain_points2);

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

      int sh = solution.height();
      double gradient = solution.get(sh - 1, 0) + solution.get(sh - 2,
                        0) - solution.get(sh - 3, 0) - solution.get(sh - 4, 0);
      angs_and_grads.push_back(ang1);
      angs_and_grads.push_back(ang2);
      angs_and_grads.push_back(gradient);
    }
  }
  io::write_ex2_gradients_to_file("output/ex2grads.txt", angs_and_grads);

  /////////////////////////////////////////////////////////////////////////


  config.pde =  ie_solver_config::Pde::STOKES;
  config.solution_dimension = 2;
  config.boundary_shape = Boundary::BoundaryShape::EX3;

  std::unique_ptr<Boundary> boundary3 =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  boundary3->initialize(config.num_boundary_points,
                        BoundaryCondition::DEFAULT);
  QuadTree quadtree3;
  quadtree3.initialize_tree(boundary3.get(), std::vector<double>(),
                            config.solution_dimension, config.domain_dimension);
  std::vector<double> domain_points3;
  // get_domain_points(config.domain_size, &domain_points3, quadtree3.min,
  //                   quadtree3.max);

  domain_points3.push_back(0.2);
  domain_points3.push_back(0.675);

  domain_points3.push_back(0.3);
  domain_points3.push_back(0.675);

  domain_points3.push_back(0.4);
  domain_points3.push_back(0.675);

  domain_points3.push_back(0.6);
  domain_points3.push_back(0.325);

  domain_points3.push_back(0.7);
  domain_points3.push_back(0.325);

  domain_points3.push_back(0.8);
  domain_points3.push_back(0.325);

  std::unique_ptr<Boundary> perturbed_boundary3 =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  perturbed_boundary3->initialize(config.num_boundary_points,
                                  config.boundary_condition);
  FRAME_CAP = 50;
  std::vector<double> ang_and_flow;
  for (int frame1 = 0; frame1 < FRAME_CAP; frame1++) {
    double ang = (frame1 / (0.0 + FRAME_CAP)) * M_PI;

    perturbed_boundary3->perturbation_parameters[0] = ang;

    perturbed_boundary3->initialize(config.num_boundary_points,
                                    config.boundary_condition);
    quadtree3.perturb(*perturbed_boundary3.get());
    ie_Mat solution = boundary_integral_solve(config, &quadtree3,
                      domain_points3);

    int sh = solution.height();
    double flow = 0.;

    for (int i = sh - 12; i < sh; i += 2) {
      flow += solution.get(i, 0);
    }
    ang_and_flow.push_back(ang);
    ang_and_flow.push_back(flow);
  }

  io::write_ex3_flows_to_file("output/ex3flows.txt", ang_and_flow);
}

}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::print_landscapes();
  return 0;
}

