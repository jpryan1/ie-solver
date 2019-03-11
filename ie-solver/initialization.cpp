// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/initialization.h"

namespace ie_solver {


// TODO(John) now points vec might need to be boundary_points instead
void Initialization::InitializeDomainKernel(ie_Mat* K,
    const std::vector<double>& domain_points, int test_size,
    Kernel* kernel, int solution_dimension) {
  // is stokes TODO
  std::vector<double> points = kernel->boundary->points;
  std::vector<double> normals = kernel->boundary->normals;
  std::vector<double> weights = kernel->boundary->weights;

  if (solution_dimension == 2) {
    Stokes_InitializeDomainKernel(K, points, normals, weights, domain_points,
                                  test_size, kernel);
    return;
  }
  // columns for phi (aka dofs), rows for spatial domain
  int dofs = points.size() / 2;
  double scale = 1.0 / (2 * M_PI);
  // omp_set_num_threads(4);
  // #pragma omp parallel for
  for (int i = 0; i < test_size * test_size; i++) {
    Dof domain_point;
    domain_point.point = Vec2(domain_points[2 * i], domain_points[2 * i + 1]);
    bool in_domain = kernel->boundary->is_in_domain(domain_point.point);
    for (int j = 0; j < dofs; j++) {
      if (!in_domain) {
        K->set(i, j, 0);
        continue;
      }
      Dof boundary_point;
      boundary_point.point = Vec2(points[2 * j], points[2 * j + 1]);
      boundary_point.normal = Vec2(normals[2 * j], normals[2 * j + 1]);
      boundary_point.weight = weights[j];
      double potential = kernel->laplace_kernel(domain_point, boundary_point).get(0,0);
      K->set(i, j, potential);
    }
  }
}


void Initialization::Electric_InitializeBoundary(ie_Mat* f,
    const std::vector<double>& points) {
  for (unsigned int i = 0; i < f->height(); i++) {
    f->set(i, 0, 1);
  }
}

}  // namespace ie_solver
