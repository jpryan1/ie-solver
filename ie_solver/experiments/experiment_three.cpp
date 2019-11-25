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
  config.pde = ie_solver_config::Pde::LAPLACE_NEUMANN;
  config.num_boundary_points = pow(2, 13);
  config.domain_size = 200;
  config.solution_dimension = 1;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX3;
  return config;
}


void get_sample_vals(const ie_solver_config& config, double* samples,
                     Boundary* boundary,
                     int perturbed_param,
                     QuadTree* quadtree,
                     const std::vector<double>& domain_points,
                     double* findiff) {
  // TODO(John) try copying quadtree and running this in parallel
  QuadTree trees[4];
  for (int i = 0; i < 4; i++) {
    quadtree->copy_into(&(trees[i]));
  }

  for (int i = 0; i < 4; i++) {
    boundary->perturbation_parameters[perturbed_param] = samples[i];
    boundary->initialize(config.num_boundary_points,
                         config.boundary_condition);
    trees[i].perturb(*boundary);
  }

  #pragma omp parallel for num_threads(4)
  for (int i = 0; i < 4; i++) {
    ie_Mat solution = boundary_integral_solve(config, &trees[i],
                      domain_points);

    double gradient = (solution.get(1, 0) - solution.get(0, 0))
                      / 0.0002;
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


void run_experiment3() {
  // double start = omp_get_wtime();
  ie_solver_config config = get_experiment_three_config();

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  boundary->initialize(config.num_boundary_points,
                       BoundaryCondition::DEFAULT);

  std::vector<double> domain_points;

  // TODO(John) the fact that round numbers screw things up is a problem -
  // out of domain should salt these maybe?
  domain_points.push_back(0.4999);
  domain_points.push_back(0.5001);

  domain_points.push_back(0.5001);
  domain_points.push_back(0.5001);

  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);

  double current_ang1 = 1;
  double current_ang2 = -1;

  perturbed_boundary->perturbation_parameters[0] = current_ang1;
  perturbed_boundary->perturbation_parameters[1] = current_ang2;
  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);
  QuadTree quadtree;
  quadtree.initialize_tree(perturbed_boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);
  double prev_gradient = (solution.get(1, 0) - solution.get(0, 0))
                         / 0.0002;
  double start_alpha = 1;
  double alpha_decay = 0.8;
  double h = 1e-4;
  int FRAME_CAP = 20;


  for (int step = 0; step < FRAME_CAP; step++) {
    // First, find gradiant.
    double findiff1[4];
    double samples1[4] = {current_ang1 - 2 * h, current_ang1 - h,
                          current_ang1 + h, current_ang1 + 2 * h
                         };
    // TODO(John) try copying quadtree and running this in parallel
    get_sample_vals(config, samples1, perturbed_boundary.get(), 0,
                    &quadtree, domain_points, findiff1);
    perturbed_boundary->perturbation_parameters[0] = current_ang1;

    double findiff2[4];
    double samples2[4] = {current_ang2 - 2 * h, current_ang2 - h,
                          current_ang2 + h, current_ang2 + 2 * h
                         };
    get_sample_vals(config, samples2, perturbed_boundary.get(), 1,
                    &quadtree, domain_points, findiff2);

    double grad1 = (findiff1[0] - 8 * findiff1[1] + 8 * findiff1[2]
                    - findiff1[3]) / (12 * h);
    double grad2 = (findiff2[0] - 8 * findiff2[1] + 8 * findiff2[2]
                    - findiff2[3]) / (12 * h);

    // Now, perform line search
    double alpha = start_alpha;
    while (alpha > 0.01) {
      double trial_ang1 = current_ang1 + alpha * grad1;
      double trial_ang2 = current_ang2 + alpha * grad2;
      enforce_separation(&trial_ang1, &trial_ang2);

      // Calculate new obj val, check wolfe cond satisfaction,
      // else update param, repeat.
      perturbed_boundary->perturbation_parameters[0] = trial_ang1;
      perturbed_boundary->perturbation_parameters[1] = trial_ang2;
      perturbed_boundary->initialize(config.num_boundary_points,
                                     config.boundary_condition);
      quadtree.perturb(*perturbed_boundary);
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);
      double gradient = (solution.get(1, 0) - solution.get(0, 0))
                        / 0.0002;

      if (prev_gradient < gradient) {
        prev_gradient = gradient;
        current_ang1 = trial_ang1;
        current_ang2 = trial_ang2;
        std::cout << gradient << std::endl;
        perturbed_boundary->perturbation_parameters[0] = current_ang1;
        perturbed_boundary->perturbation_parameters[1] = current_ang2;

        io::write_solution_to_file("output/bake/sol/" + std::to_string(
                                     step)  + ".txt", solution, domain_points,
                                   config.solution_dimension);
        io::write_boundary_to_file("output/bake/boundary/" + std::to_string(
                                     step) + ".txt",
                                   perturbed_boundary->points);
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

  // ie_Mat solution = boundary_integral_solve(config, &quadtree,
  //                   domain_points);
  // io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
  //                            domain_points,
  //                            config.solution_dimension);
  // io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
  //                            boundary->points);

}

}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment3();
  return 0;
}

