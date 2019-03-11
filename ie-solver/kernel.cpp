// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie-solver/kernel.h"

namespace ie_solver {

ie_Mat Kernel::get(const Dof& a, const Dof& b) const {
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      return laplace_kernel(a, b);
      break;
    case ie_solver_config::Pde::STOKES:
      return stokes_kernel(a, b);
      break;
  }
}


double Kernel::get(unsigned int i, unsigned int j) const {
  Dof a, b;
  unsigned int i_point_index = i / solution_dimension;
  unsigned int i_points_vec_index = i_point_index * domain_dimension;
  unsigned int j_point_index = j / solution_dimension;
  unsigned int j_points_vec_index = j_point_index * domain_dimension;

  a.point = Vec2(boundary->points[i_points_vec_index],
                 boundary->points[i_points_vec_index + 1]);
  a.normal = Vec2(boundary->normals[i_points_vec_index],
                  boundary->normals[i_points_vec_index + 1]);
  a.curvature = boundary->curvatures[i_point_index];
  a.weight = boundary->weights[i_point_index];

  b.point = Vec2(boundary->points[j_points_vec_index],
                 boundary->points[j_points_vec_index + 1]);
  b.normal = Vec2(boundary->normals[j_points_vec_index],
                  boundary->normals[j_points_vec_index + 1]);
  b.curvature = boundary->curvatures[j_point_index];
  b.weight = boundary->weights[j_point_index];

  ie_Mat tensor = get(a, b);
  return tensor.get(i % solution_dimension, j % solution_dimension);
}


ie_Mat Kernel::stokes_kernel(const Dof& a, const Dof& b) const {
  // double layer
  double scale = 1.0 / (M_PI);

  ie_Mat tensor(2, 2);

  if (a.point.a[0] == b.point.a[0] && a.point.a[1] == b.point.a[1]) {
    // tangent
    double t0 = -a.normal.a[1];
    double t1 =  a.normal.a[0];
    double potential = - 0.5 * a.curvature * a.weight * scale;
    tensor.set(0, 0, -0.5 + potential * t0 * t0);
    tensor.set(1, 1, -0.5 + potential * t1 * t1);
    tensor.set(0, 1, potential * t0 * t1);
    tensor.set(1, 0, potential * t1 * t0);
    return tensor;
  }

  Vec2 r = a.point - b.point;

  double r0 = r.a[0];
  double r1 = r.a[1];
  double potential = b.weight * scale * (r.dot(b.normal)) / (pow(r.dot(
                       r), 2));

  tensor.set(0, 0, potential * r0 * r0);
  tensor.set(1, 1, potential * r1 * r1);
  tensor.set(0, 1, potential * r0 * r1);
  tensor.set(1, 0, potential * r1 * r0);
  return tensor;
}

ie_Mat Kernel::laplace_kernel(const Dof& a, const Dof& b) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat tensor(1, 1);
  if (a.point.a[0] == b.point.a[0] && a.point.a[1] == b.point.a[1]) {
    tensor.set(0, 0, 0.5 + 0.5 * a.curvature * a.weight * scale);
    return tensor;
  }
  Vec2 r = a.point - b.point;
  tensor.set(0, 0, -b.weight * scale * (r.dot(b.normal)) / (r.dot(r)));
  return tensor;
}

// This function stores the DoF data,  and calculates the diagonals of the mat
void Kernel::load(Boundary* boundary_, ie_solver_config::Pde pde_,
                  int solution_dimension_, int domain_dimension_) {
  solution_dimension = solution_dimension_;
  domain_dimension = domain_dimension_;
  boundary = boundary_;
  pde = pde_;
}


// TODO(John) shouldn't this->I have the underscore after it, not this arg?
ie_Mat Kernel::operator()(const std::vector<unsigned int>& I_,
                          const std::vector<unsigned int>& J_) const {
  ie_Mat ret(I_.size(), J_.size());

  int olda_ = I_.size();
  for (unsigned int i = 0; i < I_.size(); i++) {
    for (unsigned int j = 0; j < J_.size(); j++) {
      ret.mat[i + olda_ * j] = get(I_[i], J_[j]);
    }
  }

  return ret;
}

}  // namespace ie_solver
