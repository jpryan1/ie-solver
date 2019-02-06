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
  bool is_stokes = (config.pde == ie_solver_config::STOKES);
  bool is_time_trial = (skel_times != nullptr);
  double id_tol = config.id_tol;

  // TODO(John) insert comment here explaining why this is necessary
  // if(!timing && N > 10000){
  //  printf("Turn down N or disable accuracy checking please\n");
  //  return;
  // }
  config.boundary->perturbation_size = 0;
  config.boundary->initialize(config.N, config.boundary_condition);

  if (!is_time_trial) {
    write_boundary_to_file(config.boundary->points);
  }
  // TODO(John) why not just N here?
  int dofs = config.boundary->points.size() / 2;
  QuadTree quadtree;
  quadtree.initialize_tree(config.boundary.get(), is_stokes);
  if (!is_time_trial) {
    quadtree.write_quadtree_to_file();
  }

  // Consider making init instead of constructor for readability
  // TODO(John) why do the tools AND the tree need the points and normals and
  // weights again?
  bool strong_admissibility =
    (config.admissibility == ie_solver_config::STRONG);
  IeSolverTools ie_solver_tools(id_tol, strong_admissibility, is_stokes);

  Kernel kernel;

  kernel.load(config.boundary.get(), is_stokes);

  // Now we calculate the solution to the PDE inside the domain by setting
  // up the relevant linear system.

  int dim = is_stokes ? 2 : 1;

  std::vector<double> domain_points;
  get_domain_points(&domain_points, quadtree.min, quadtree.max);

  Initialization init;

  ie_Mat f(dim * dofs, 1);
  // TODO(John) get rid of these damn if(is_stokes) statements, push them to the
  // functions
  if (is_stokes) {
    init.Stokes_InitializeBoundary(&f, config.boundary->normals);
    // notice here we are passing the normals since the flow will just be
    // unit tangent to the boundary.
  } else {
    f = config.boundary->boundary_values;
  }

  ie_Mat phi(dim * dofs, 1);

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

  ie_solver_tools.skeletonize(kernel, &quadtree);
  ie_solver_tools.check_factorization_against_kernel(kernel, &quadtree);

  std::vector<double> old_points = config.boundary->points;
  config.boundary->perturbation_size = 100;
  config.boundary->initialize(config.N, config.boundary_condition);
  f = config.boundary->boundary_values;
  std::cout << "\nPerturbed" << std::endl;
  std::vector<double> new_points = config.boundary->points;

  // Quadtree needs editing now
  // quadtree.perturb(old_points, new_points);

  quadtree.reset();
  ie_solver_tools.skeletonize(kernel, &quadtree);
  ie_solver_tools.check_factorization_against_kernel(kernel, &quadtree);

  ie_Mat phi1(dim * kernel.boundary->points.size(), 1);
  ie_solver_tools.sparse_matvec(kernel, quadtree, f, &phi1);

  quadtree.reset();
  std::cout << "\nReset" << std::endl;

  ie_solver_tools.skeletonize(kernel, &quadtree);
  ie_solver_tools.check_factorization_against_kernel(kernel, &quadtree);
  ie_Mat phi2(dim * kernel.boundary->points.size(), 1);
  ie_solver_tools.sparse_matvec(kernel, quadtree, f, &phi2);

  phi2 -= phi1;
  std::cout << "\nDiff norm: " << phi2.norm2() << std::endl;
}

}  // namespace ie_solver

int main(int argc, char** argv) {
  srand(omp_get_wtime());

  ie_solver::ie_solver_config config;
  if (!ie_solver::parse_input_into_config(argc, argv, &config)) {
    return 1;
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

