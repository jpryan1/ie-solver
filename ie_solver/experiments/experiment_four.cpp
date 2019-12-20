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
#include "ie_solver/boundaries/ex4boundary.h"

namespace ie_solver {


ie_solver_config get_experiment_four_config() {
  ie_solver_config config;
  config.id_tol = 1e-6;
  config.pde = ie_solver_config::Pde::STOKES;
  config.num_boundary_points = pow(2, 12);
  config.domain_size = 400;
  config.domain_dimension = 2;
  config.solution_dimension = 2;
  config.boundary_condition = BoundaryCondition::DEFAULT;
  config.boundary_shape = Boundary::BoundaryShape::EX4;
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
  double start = omp_get_wtime();
  for (int i = 0; i < 4; i++) {
    quadtree.copy_into(&(trees[i]));
  }

  for (int i = 0; i < 4; i++) {
    boundary->perturbation_parameters[perturbed_param] = samples[i];
    boundary->initialize(config.num_boundary_points,
                         config.boundary_condition);
    trees[i].perturb(*boundary);
  }
  double end = omp_get_wtime();
  std::cout << "copy took " << (end - start) << std::endl;

  // #pragma omp parallel for num_threads(2)
  for (int i = 0; i < 4; i++) {
    ie_Mat solution = boundary_integral_solve(config, &trees[i],
                      domain_points);

    double gradient = solution.get(0,0);
    findiff[i] = gradient;
  }
}

void run_experiment4() {
  // double start = omp_get_wtime();
  ie_solver_config config = get_experiment_four_config();

  // We'll iteratively reinitialized another Boundary and use that
  // to update the quadtree's Boundary.
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex4Boundary());

  boundary->initialize(config.num_boundary_points,
                       config.boundary_condition);

  double current_ang1 = 0.5;

  boundary->perturbation_parameters[0] = current_ang1;
  boundary->initialize(config.num_boundary_points,
                       config.boundary_condition);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), std::vector<double>(),
                           config.solution_dimension, config.domain_dimension);

  std::vector<double> domain_points;

  std::unique_ptr<Boundary> perturbed_boundary =
    std::unique_ptr<Boundary>(new Ex4Boundary());

  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);
  perturbed_boundary->perturbation_parameters[0] = current_ang1;
  perturbed_boundary->initialize(config.num_boundary_points,
                                 config.boundary_condition);

  // get_domain_points(config.domain_size, &domain_points, quadtree.min,
  //                   quadtree.max, quadtree.min, quadtree.max);
  // TODO(John) the fact that round numbers screw things up is a problem -
  // out of domain should salt these maybe?
  domain_points.push_back(0.7);
  domain_points.push_back(0.5001);

  ie_Mat solution = boundary_integral_solve(config, &quadtree,
                    domain_points);

  double prev_obj_fun = solution.get(0, 0);
  double start_alpha = 1;
  double alpha_decay = 0.8;
  double h = 1e-4;
  int FRAME_CAP = 15;

  double prev_step_start = omp_get_wtime();
  for (int step = 0; step < FRAME_CAP; step++) {

    double next_step_start = omp_get_wtime();
    if (step > 0) {
      std::cout << "Step took " << (next_step_start - prev_step_start) << std::endl;
    }
    prev_step_start = next_step_start;
    // First, find gradient.
    double findiff1[4];
    double samples1[4] = {current_ang1 - 2 * h, current_ang1 - h,
                          current_ang1 + h, current_ang1 + 2 * h
                         };
      get_sample_vals(config, samples1, perturbed_boundary.get(), 0,
                    quadtree, domain_points, findiff1);
    perturbed_boundary->perturbation_parameters[0] = current_ang1;

    double grad1 = (findiff1[0] - 8 * findiff1[1] + 8 * findiff1[2]
                    - findiff1[3]) / (12 * h);
   
    // Now, perform line search
    double alpha = start_alpha;
    while (alpha > 0.01) {
      double trial_ang1 = current_ang1 - alpha * grad1;

      // Calculate new obj val, check wolfe cond satisfaction,
      // else update param, repeat.
      perturbed_boundary->perturbation_parameters[0] = trial_ang1;
      perturbed_boundary->initialize(config.num_boundary_points,
                                     config.boundary_condition);
      quadtree.perturb(*perturbed_boundary);
      ie_Mat solution = boundary_integral_solve(config, &quadtree,
                        domain_points);
      double cur_obj_fun = solution.get(0, 0);
std::cout << "line srch obj function val "
                  << cur_obj_fun <<" vs " <<prev_obj_fun<<std::endl;
      
      if (prev_obj_fun > cur_obj_fun) {
        prev_obj_fun = cur_obj_fun;
        current_ang1 = trial_ang1;
        std::cout << "Current obj function val "
                  << cur_obj_fun << std::endl;
        perturbed_boundary->perturbation_parameters[0] = current_ang1;

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

  io::write_solution_to_file("output/data/ie_solver_solution.txt", solution,
                             domain_points,
                             config.solution_dimension);
  io::write_boundary_to_file("output/data/ie_solver_boundary.txt",
                             boundary->points);

}


}  // namespace ie_solver


int main(int argc, char** argv) {
  srand(0);  // omp_get_wtime());
  ie_solver::run_experiment4();
  return 0;
}

