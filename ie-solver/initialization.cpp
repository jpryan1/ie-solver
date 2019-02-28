// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/initialization.h"

namespace ie_solver {


// TODO(John) now points vec might need to be boundary_points instead
void Initialization::InitializeDomainKernel(ie_Mat* K,
    const std::vector<double>& domain_points, int test_size,
    Boundary* boundary, bool is_stokes) {
  // is stokes TODO
  std::vector<double> points = boundary->points;
  std::vector<double> normals = boundary->normals;
  std::vector<double> weights = boundary->weights;

  if (is_stokes) {
    Stokes_InitializeDomainKernel(K, points, normals, weights, domain_points,
                                  test_size, boundary);
    return;
  }

  // columns for phi (aka dofs), rows for spatial domain
  int dofs = points.size() / 2;
  double scale = 1.0 / (2 * M_PI);
  // omp_set_num_threads(4);
  // #pragma omp parallel for
  for (int i = 0; i < test_size * test_size; i++) {
    Vec2 x(domain_points[2 * i], domain_points[2 * i + 1]);
    bool in_domain = boundary->is_in_domain(x);
    for (int j = 0; j < dofs; j++) {
      Vec2 y(points[2 * j], points[2 * j + 1]);

      if (!in_domain) {
        K->set(i, j, 0);
        continue;
      }
      Vec2 r = x - y;
      Vec2 n(normals[2 * j], normals[2 * j + 1]);
      double potential = -weights[j] * scale * (r.dot(n)) / (r.dot(r));
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
