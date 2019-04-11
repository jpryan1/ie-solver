// Copyright 2019 John Paul Ryan
#include <cassert>
#include <cmath>
#include "ie-solver/initialization.h"

namespace ie_solver {
// TODO(John) - subclasses for different PDES

void Initialization::Stokes_InitializeDomainKernel(ie_Mat* K,
    const std::vector<double>& points,
    const std::vector<double>& normals,
    const std::vector<double>& weights,
    const std::vector<double>& domain_points, int test_size,
    Kernel* kernel) {
  // columns for phi (aka dofs), rows for spatial domain

  assert(points.size() == normals.size() &&
         "Points and normals must have same size in Stokes domain init.");
  assert(points.size() == 2 * weights.size() &&
         "In dim 2, pts must be 2*size of weights in stokes domain init.");

  for (int i = 0; i < test_size * test_size; i++) {
    Vec2 x(domain_points[2 * i], domain_points[2 * i + 1]);
    bool in_domain = kernel->boundary->is_in_domain(x);
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
      ie_Mat tensor = kernel->stokes_kernel(a, b);

      K->set(2 * i  , j  , tensor.get(0, 0));
      K->set(2 * i + 1, j  , tensor.get(1, 0));
      K->set(2 * i  , j + 1, tensor.get(0, 1));
      K->set(2 * i + 1, j + 1, tensor.get(1, 1));
    }
  }
}


}  // namespace ie_solver
