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
  int domain_dimension = 2;
  int solution_dimension = 1;
  if (config.pde == ie_solver_config::STOKES) {
    solution_dimension = 2;
  }

  bool is_time_trial = (skel_times != nullptr);
  double id_tol = config.id_tol;

  if (!is_time_trial) {
    write_boundary_to_file(config.boundary->points);
  }


  QuadTree quadtree;
  quadtree.initialize_tree(config.boundary.get(), std::vector<double>(),
                           solution_dimension,
                           domain_dimension);

  if (!is_time_trial) {
    quadtree.write_quadtree_to_file();
  }

  // Consider making init instead of constructor for readability
  bool strong_admissibility =
    (config.admissibility == ie_solver_config::STRONG);
  IeSolverTools ie_solver_tools(id_tol, strong_admissibility,
                                solution_dimension, domain_dimension);

  Kernel kernel;
  kernel.load(config.boundary.get(), config.pde, solution_dimension,
              domain_dimension);

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
  ie_Mat phi(config.num_boundary_points * solution_dimension, 1);
  ie_solver_tools.skeletonize(kernel, &quadtree);
  ie_solver_tools.solve(kernel, quadtree, &phi, f);

  // This will be done as a sparse mat vec in the future, for now we do
  // dense matvec
  std::vector<double> domain_points;
  get_domain_points(&domain_points, quadtree.min, quadtree.max);
  QuadTree boundary_to_domain;
  boundary_to_domain.initialize_tree(config.boundary.get(),
                                     domain_points, solution_dimension,
                                     domain_dimension);
  // ie_solver_tools.b2dskeletonize(kernel, &boundary_to_domain);
  ie_Mat domain(TEST_SIZE * TEST_SIZE * solution_dimension, 1);

  // ie_solver_tools.b2dsparse_matvec(kernel, boundary_to_domain, phi, &domain);

  ie_Mat K_domain(TEST_SIZE * TEST_SIZE * solution_dimension,
                  config.num_boundary_points * solution_dimension);
  Initialization init;
  init.InitializeDomainKernel(&K_domain,
                              domain_points, TEST_SIZE, &kernel,
                              solution_dimension);
  ie_Mat::gemv(NORMAL, 1., K_domain, phi, 0., &domain);

  write_solution_to_file("output/data/ie_solver_solution.txt", domain,
                         domain_points, solution_dimension);
  switch (config.pde) {
    case ie_solver_config::LAPLACE:
      check_laplace_solution(domain, config.id_tol, domain_points,
                             config.boundary.get());
      break;
    case ie_solver_config::STOKES:
      check_stokes_solution(domain, config.id_tol, domain_points,
                            config.boundary.get());
      break;
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

// First we init the boundary so we can correct the num_boundary_points
  config.boundary->initialize(config.num_boundary_points,
                              config.boundary_condition);
  config.num_boundary_points = config.boundary->weights.size();

  if (!config.scaling) {
    ie_solver::boundary_integral_solve(config);
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

