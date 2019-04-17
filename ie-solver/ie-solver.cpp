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

// First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();

  if (!config.scaling) {
    double error = ie_solver::boundary_integral_solve(
                     config);
    std::cout << "Error: " << error << std::endl;;
  }

  if (config.scaling) {
    std::vector<double> n_times, eps_times;
    int scale_n[] = {5000, 6000, 7000, 8000, 9000};
    double scale_eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};

    for (int n : scale_n) {
      config.num_boundary_points = n;
      config.id_tol = DEFAULT_ID_TOL;
      ie_solver::boundary_integral_solve(config, &n_times);
    }

    for (double eps : scale_eps) {
      config.num_boundary_points = DEFAULT_NUM_DISCRETIZATION_POINTS * 10;
      config.id_tol = eps;
      ie_solver::boundary_integral_solve(config, &eps_times);
    }

    ie_solver::write_times_to_files(scale_n, n_times, scale_eps, eps_times);
  }

  return 0;
}

