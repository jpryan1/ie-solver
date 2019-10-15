// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include "ie_solver/initialization.h"

namespace ie_solver {


// TODO(John) now points vec might need to be boundary_points instead
void Initialization::InitializeDomainKernel(ie_Mat* K,
    const std::vector<double>& domain_points, int test_size,
    const Kernel& kernel, int solution_dimension) {
  double domain_init_start = omp_get_wtime();
  double domain_init_end;
  // is stokes TODO
  std::vector<double> points = kernel.boundary->points;
  std::vector<double> normals = kernel.boundary->normals;
  std::vector<double> weights = kernel.boundary->weights;
  if (kernel.pde == ie_solver_config::Pde::STOKES) {
    Stokes_InitializeDomainKernel(K, points, normals, weights, domain_points,
                                  test_size, kernel);
    domain_init_end = omp_get_wtime();
    std::cout << "timing: domain_init " << (domain_init_end - domain_init_start) <<
              std::endl;
    return;
  } else if (kernel.pde == ie_solver_config::Pde::LAPLACE_NEUMANN) {
    LaplaceNeumann_InitializeDomainKernel(K, points, normals, weights,
                                          domain_points,
                                          test_size, kernel);
    domain_init_end = omp_get_wtime();
    std::cout << "timing: domain_init " << (domain_init_end - domain_init_start) <<
              std::endl;
    return;
  }

  // columns for phi (aka dofs), rows for spatial domain
  int dofs = points.size() / 2;
  for (int i = 0; i < test_size * test_size; i++) {
    Dof domain_point;
    domain_point.is_boundary = false;
    domain_point.point = Vec2(domain_points[2 * i], domain_points[2 * i + 1]);
    bool in_domain = kernel.boundary->is_in_domain(domain_point.point);

    for (int j = 0; j < dofs; j++) {
      if (!in_domain) {
        K->set(i, j, 0);
        continue;
      }
      Dof boundary_point;
      boundary_point.is_boundary = true;
      boundary_point.point = Vec2(points[2 * j], points[2 * j + 1]);
      boundary_point.normal = Vec2(normals[2 * j], normals[2 * j + 1]);
      boundary_point.weight = weights[j];
      double potential =
        kernel.laplace_kernel(domain_point, boundary_point).get(0, 0);
      K->set(i, j, potential);
    }
  }
  domain_init_end = omp_get_wtime();
  std::cout << "timing: domain_init " << (domain_init_end - domain_init_start) <<
            std::endl;
}


void Initialization::LaplaceNeumann_InitializeDomainKernel(ie_Mat* K,
    const std::vector<double>& points,
    const std::vector<double>& normals,
    const std::vector<double>& weights,
    const std::vector<double>& domain_points, int test_size,
    const Kernel& kernel) {

  // columns for phi (aka dofs), rows for spatial domain
  int dofs = points.size() / 2;
  for (int i = 0; i < test_size * test_size; i++) {
    Dof domain_point;
    domain_point.is_boundary = false;
    domain_point.point = Vec2(domain_points[2 * i], domain_points[2 * i + 1]);
    bool in_domain = kernel.boundary->is_in_domain(domain_point.point);

    for (int j = 0; j < dofs; j++) {
      if (!in_domain) {
        K->set(i, j, 0);
        continue;
      }
      Dof boundary_point;
      boundary_point.is_boundary = true;
      boundary_point.point = Vec2(points[2 * j], points[2 * j + 1]);
      boundary_point.normal = Vec2(normals[2 * j], normals[2 * j + 1]);
      boundary_point.weight = weights[j];
      double potential =
        kernel.laplace_neumann_kernel_forward(domain_point, boundary_point).get(0, 0);
      K->set(i, j, potential);
    }
  }
}


void Initialization::Stokes_InitializeDomainKernel(ie_Mat* K,
    const std::vector<double>& points,
    const std::vector<double>& normals,
    const std::vector<double>& weights,
    const std::vector<double>& domain_points, int test_size,
    const Kernel& kernel) {
  // columns for phi (aka dofs), rows for spatial domain

  assert(points.size() == normals.size() &&
         "Points and normals must have same size in Stokes domain init.");
  assert(points.size() == 2 * weights.size() &&
         "In dim 2, pts must be 2*size of weights in stokes domain init.");
  #pragma omp parallel for num_threads(8)
  for (int i = 0; i < test_size * test_size; i++) {
    Vec2 x(domain_points[2 * i], domain_points[2 * i + 1]);

    bool in_domain = kernel.boundary->is_in_domain(x);
    for (unsigned int j = 0; j < points.size(); j += 2) {
      if (!in_domain) {
        K->set(2 * i  , j  , 0);
        K->set(2 * i + 1, j  , 0);
        K->set(2 * i  , j + 1, 0);
        K->set(2 * i + 1, j + 1, 0);
        continue;
      }

      Vec2 y(points[j], points[j + 1]);

      Dof a, b;
      a.is_boundary = false;
      a.point = x;
      b.is_boundary = true;
      b.point = y;
      b.normal = Vec2(normals[j], normals[j + 1]);
      b.weight = weights[j / 2];
      ie_Mat tensor = kernel.stokes_kernel(a, b);

      K->set(2 * i  , j  , tensor.get(0, 0));
      K->set(2 * i + 1, j  , tensor.get(1, 0));
      K->set(2 * i  , j + 1, tensor.get(0, 1));
      K->set(2 * i + 1, j + 1, tensor.get(1, 1));
    }
  }
}


}  // namespace ie_solver
