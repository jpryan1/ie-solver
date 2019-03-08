// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <string.h>
#include <fstream>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include "ie-solver/ie_mat.h"
#include "ie-solver/initialization.h"
#include "ie-solver/tools/ie_solver_tools.h"
#include "ie-solver/quadtree.h"
#include "ie-solver/kernel.h"
#include "ie-solver/log.h"
#include "ie-solver/helpers.h"

#define TIMING_ITERATIONS 5

namespace ie_solver {

void boundary_integral_solve(const ie_solver_config& config,
                             std::vector<double>* skel_times = nullptr) {
  bool is_stokes = false;
  if (config.pde == ie_solver_config::STOKES) {
    is_stokes = true;
  }

  bool is_time_trial = (skel_times != nullptr);
  double id_tol = config.id_tol;

  // config.boundary->perturbation_size = 50;
  config.boundary->initialize(config.N, config.boundary_condition);

  if (!is_time_trial) {
    write_boundary_to_file(config.boundary->points);
  }
  int dofs = config.boundary->points.size() / 2;
  if (is_stokes) dofs *= 2;

  QuadTree quadtree;
  quadtree.initialize_tree(config.boundary.get(), is_stokes);

  if (!is_time_trial) {
    quadtree.write_quadtree_to_file();
  }

  // Consider making init instead of constructor for readability
  bool strong_admissibility =
    (config.admissibility == ie_solver_config::STRONG);
  IeSolverTools ie_solver_tools(id_tol, strong_admissibility, is_stokes);

  Kernel kernel;
  kernel.load(config.boundary.get(), is_stokes);
  if (is_time_trial) {
    double elapsed = 0;
    for (int i = 0; i < TIMING_ITERATIONS; i++) {
      double start = omp_get_wtime();
      ie_solver_tools.skeletonize(kernel, &quadtree);
      double end = omp_get_wtime();
      elapsed += (end - start);
      quadtree.reset();
    }
    skel_times->push_back(elapsed / TIMING_ITERATIONS);
    return;
  }

  ie_Mat f = config.boundary->boundary_values;
  ie_Mat phi(dofs, 1);
  // std::vector<unsigned int> all;
  // for (int i = 0; i < config.boundary->points.size(); i++) {
  //   if (!is_stokes && i == config.boundary->points.size() / 2) {
  //     break;
  //   }
  //   all.push_back(i);
  // }

  // ie_Mat kern = kernel(all, all);

  // kern.left_multiply_inverse(f, &phi);
  ie_solver_tools.skeletonize(kernel, &quadtree);
  ie_solver_tools.solve(kernel, quadtree, &phi, f);

  // This will be done as a sparse mat vec in the future, for now we do
  // dense matvec

  ie_Mat K_domain, domain;
  if (is_stokes) {
    K_domain = ie_Mat(2 * TEST_SIZE * TEST_SIZE, dofs);
    domain = ie_Mat(2 * TEST_SIZE * TEST_SIZE, 1);
  } else {
    K_domain = ie_Mat(TEST_SIZE * TEST_SIZE, dofs);
    domain = ie_Mat(TEST_SIZE * TEST_SIZE, 1);
  }

  Initialization init;
  std::vector<double> domain_points;
  get_domain_points(&domain_points, quadtree.min, quadtree.max);

  init.InitializeDomainKernel(&K_domain,
                              domain_points, TEST_SIZE, &kernel,
                              is_stokes);

  ie_Mat::gemv(NORMAL, 1., K_domain, phi, 0., &domain);

  write_solution_to_file("output/data/ie_solver_solution.txt", domain,
                         domain_points, is_stokes);
  if (!is_stokes) {
    check_laplace_solution(domain, config.id_tol, domain_points,
                           config.boundary.get());
  } else {
    check_stokes_solution(domain, config.id_tol, domain_points,
                           config.boundary.get());
  }
}

}  // namespace ie_solver

int main(int argc, char** argv) {
  // TODO(John) allow for command line args for setting parameters
  srand(0);  // omp_get_wtime());

  ie_solver::ie_solver_config config;
  if (!ie_solver::parse_input_into_config(argc, argv, &config)) {
    return 1;
  }

  if (config.boundary_condition
      ==  ie_solver::Boundary::BoundaryCondition::STOKES
      || config.pde == ie_solver::ie_solver_config::STOKES) {
    config.boundary_condition =  ie_solver::Boundary::BoundaryCondition::STOKES;
    config.pde = ie_solver::ie_solver_config::STOKES;
  }


  if (!config.scaling) {
    ie_solver::boundary_integral_solve(config);
  }

  if (config.scaling) {
    std::vector<double> n_times, eps_times;
    int scale_n[] = {5000, 6000, 7000, 8000, 9000};
    double scale_eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};

    for (int n : scale_n) {
      config.N = n;
      config.id_tol = DEFAULT_ID_TOL;
      ie_solver::boundary_integral_solve(config, &n_times);
    }

    for (double eps : scale_eps) {
      config.N = DEFAULT_NUM_DISCRETIZATION_POINTS * 10;
      config.id_tol = eps;
      ie_solver::boundary_integral_solve(config, &eps_times);
    }

    ie_solver::write_times_to_files(scale_n, n_times, scale_eps, eps_times);
  }

  return 0;
}

