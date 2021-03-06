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

#define THREE_B 0
#define DELTA_X 0.0001

namespace ie_solver {

ie_solver_config get_experiment_three_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.num_boundary_points = pow(2, 15);  // 15 for timing, 12 for conv, img
  config.domain_size =  200;
  config.num_threads = 4;
  if (THREE_B) {
    config.boundary_condition = BoundaryCondition::EX3B;
    config.solution_dimension = 2;
    config.pde = ie_solver_config::Pde::STOKES;
  } else {
    config.boundary_condition = BoundaryCondition::EX3A;
    config.solution_dimension = 1;
    config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  }
  config.boundary_shape = Boundary::BoundaryShape::EX3;
  return config;
}


void get_sample_vals(const ie_solver_config& config, double* samples,
                     Boundary* boundary,
                     int perturbed_param,
                     const QuadTree& quadtree,
                     const std::vector<double>& domain_points,
                     double* findiff) {
  // TODO(John) try copying quadtree and running this in parallel
  QuadTree trees[4];
  for (int i = 0; i < 4; i++) {
    quadtree.copy_into(&(trees[i]));
  }

  for (int i = 0; i < 4; i++) {
    boundary->perturbation_parameters[perturbed_param] = samples[i];
    boundary->initialize(config.num_boundary_points,
                         config.boundary_condition);
    trees[i].perturb(*boundary);
  }

  int num_outer_threads;
  switch (config.num_threads) {
    case 1:
      num_outer_threads = 4;
      break;
    case 2:
      num_outer_threads = 2;
      break;
    default:
      num_outer_threads = 1;
      break;
  }
  #pragma omp parallel for num_threads(num_outer_threads)
  for (int i = 0; i < 4; i++) {
    ie_Mat solution = boundary_integral_solve(config, &trees[i],
                      domain_points);

    double gradient;
    if (THREE_B) {
      gradient = -solution.get(0, 0);
    } else {
      gradient = (solution.get(1, 0) - solution.get(0, 0))
                 / (2 * DELTA_X);
    }
    findiff[i] = gradient;
  }
}


void enforce_separation(double* ang1, double* ang2) {
  double current_ang1 = *ang1;
  double current_ang2 = *ang2;

  while (current_ang1 > 2 * M_PI) current_ang1 -= 2 * M_PI;
  while (current_ang1 < 0) current_ang1 += 2 * M_PI;
  while (current_ang2 > 2 * M_PI) current_ang2 -= 2 * M_PI;
  while (current_ang2 < 0) current_ang2 += 2 * M_PI;

  double* lowerang;
  double* upperang;
  if (current_ang1 < current_ang2) {
    lowerang = &current_ang1;
    upperang = &current_ang2;
  } else {
    lowerang = &current_ang2;
    upperang = &current_ang1;
  }
  double dist = std::min(*upperang - *lowerang,
                         *lowerang + 2 * M_PI - *upperang);
  if (dist < M_PI / 4.) {
    double prob = (M_PI / 4.) - dist;
    if (*upperang - *lowerang < *lowerang + 2 * M_PI - *upperang) {
      *upperang += prob;
      *lowerang -= prob;
    } else {
      *upperang -= prob;
      *lowerang += prob;
    }
  }
  *ang1 = current_ang1;
  *ang2 = current_ang2;
}


void run_experiment3a(int num_inner_threads) {
  // double start = omp_get_wtime();
  ie_solver_config config = get_experiment_three_config();
  config.num_threads = num_inner_threads;
  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());

  boundary->initialize(config.num_boundary_points,
                       config.boundary_condition);

  double current_ang1 = -2.879;
  double current_ang2 = 0.010;
  if (THREE_B) {
    current_ang1 = 1.4;
    current_ang2 = 3;
  }


  boundary->perturbation_parameters[0] = current_ang1;
  boundary->perturbation_parameters[1] = current_ang2;
  boundary->initialize(config.num_boundary_points,
                       config.boundary_condition);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  std::vector<double> domain_points;

  // TODO(John) the fact that round numbers screw things up is a problem -
  // out of domain should salt these maybe?
  domain_points.push_back(0.5 - DELTA_X);
  domain_points.push_back(0.5 + DELTA_X);

  if (!THREE_B) {
    domain_points.push_back(0.5 + DELTA_X);
    domain_points.push_back(0.5 + DELTA_X);
  }
  // get_domain_points(config.domain_size, &domain_points,
  //                   quadtree.min, quadtree.max, quadtree.min, quadtree.max);
  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);

  double prev_gradient;
  if (THREE_B) {
    prev_gradient = -solution.get(0, 0);
  } else {
    prev_gradient = (solution.get(1, 0) - solution.get(0, 0))
                    / (2 * DELTA_X);
  }
  double start_alpha = 1;
  double alpha_decay = 0.8;
  double h = 1e-4;
  int FRAME_CAP = 30;

  for (int step = 0; step < FRAME_CAP; step++) {
    // First, find gradient.
    double prev_time = omp_get_wtime();


    double findiff1[4];
    double samples1[4] = {current_ang1 - 2 * h, current_ang1 - h,
                          current_ang1 + h, current_ang1 + 2 * h
                         };
    double findiff2[4];
    double samples2[4] = {current_ang2 - 2 * h, current_ang2 - h,
                          current_ang2 + h, current_ang2 + 2 * h
                         };
    get_sample_vals(config, samples1, boundary.get(), 0,
                    quadtree, domain_points, findiff1);
    boundary->perturbation_parameters[0] = current_ang1;

    get_sample_vals(config, samples2, boundary.get(), 1,
                    quadtree, domain_points, findiff2);

    double grad1 = (findiff1[0] - 8 * findiff1[1] + 8 * findiff1[2]
                    - findiff1[3]) / (12 * h);
    double grad2 = (findiff2[0] - 8 * findiff2[1] + 8 * findiff2[2]
                    - findiff2[3]) / (12 * h);
    double curr_time = omp_get_wtime();
    // std::cout << curr_time - prev_time << std::endl;
    // Now, perform line search
    double alpha = start_alpha;
    while (alpha > 0.01) {
      double trial_ang1 = current_ang1 + alpha * grad1;
      double trial_ang2 = current_ang2 + alpha * grad2;
      enforce_separation(&trial_ang1, &trial_ang2);

      // Calculate new obj val, check wolfe cond satisfaction,
      // else update param, repeat.

      boundary->perturbation_parameters[0] = trial_ang1;
      boundary->perturbation_parameters[1] = trial_ang2;
      boundary->initialize(config.num_boundary_points,
                           config.boundary_condition);
      quadtree.perturb(*boundary);
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);
      double gradient;
      if (THREE_B) {
        gradient = -solution.get(0, 0);
      } else {
        gradient = (solution.get(1, 0) - solution.get(0, 0))
                   / (2 * DELTA_X);
      }
      if (prev_gradient < gradient) {
        prev_gradient = gradient;
        current_ang1 = trial_ang1;
        current_ang2 = trial_ang2;
        std::cout << "theta1: " << trial_ang1 << " theta2: " << trial_ang2
                  << " obj: " << gradient << std::endl;
        boundary->perturbation_parameters[0] = current_ang1;
        boundary->perturbation_parameters[1] = current_ang2;

        io::write_solution_to_file("output/bake/sol/" + std::to_string(
                                     step)  + ".txt", solution, domain_points,
                                   config.solution_dimension);
        io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                     step) + ".txt",
                                   boundary->points);
        io::write_quadtree_to_file("output/bake/tree/" + std::to_string(
                                     step)  + ".txt", quadtree);
        break;
      } else {
        alpha *= alpha_decay;
      }
    }
    if (alpha <= 0.01) {
      std::cout << "Line search did not terminate" << std::endl;
      exit(0);
    }
  }
  // io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
  //                            domain_points,
  //                            config.solution_dimension);
  // io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
  //                            boundary->points);
}

}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  omp_set_nested(1);
  ie_solver::run_experiment3a(1);

  // std::cout << "inner_threads " << 1 << std::endl;
  // for (int k = 0; k < 3; k++) {
  //   std::cout << "k " << k << std::endl;
  //   ie_solver::run_experiment3a(1);
  // }
  // std::cout << "inner_threads " << 2 << std::endl;

  // for (int k = 0; k < 3; k++) {
  //   std::cout << "k " << k << std::endl;

  //   ie_solver::run_experiment3a(2);
  // }
  // std::cout << "inner_threads " << 4 << std::endl;

  // for (int k = 0; k < 3; k++) {
  //   std::cout << "k " << k << std::endl;

  //   ie_solver::run_experiment3a(4);
  // }
  return 0;
}

