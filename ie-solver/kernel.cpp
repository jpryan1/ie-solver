// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie-solver/kernel.h"

namespace ie_solver {

ie_Mat Kernel::get(const Dof& tgt, const Dof& src) const {
  ie_Mat tensor;
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      tensor = laplace_kernel(tgt, src);
      break;
    case ie_solver_config::Pde::STOKES:
      tensor = stokes_kernel(tgt, src);
      break;
  }
  return tensor;
}


double Kernel::forward_get(unsigned int tgt_ind, unsigned int src_ind) const {
  Dof tgt, src;
  unsigned int i_point_index = tgt_ind / solution_dimension;
  unsigned int i_points_vec_index = i_point_index * domain_dimension;
  unsigned int j_point_index = src_ind / solution_dimension;
  unsigned int j_points_vec_index = j_point_index * domain_dimension;

  tgt.point = Vec2(domain_points[i_points_vec_index],
                   domain_points[i_points_vec_index + 1]);
  tgt.is_boundary = false;

  src.point = Vec2(boundary->points[j_points_vec_index],
                   boundary->points[j_points_vec_index + 1]);
  src.normal = Vec2(boundary->normals[j_points_vec_index],
                    boundary->normals[j_points_vec_index + 1]);
  src.curvature = boundary->curvatures[j_point_index];
  src.weight = boundary->weights[j_point_index];
  src.is_boundary = true;
  ie_Mat tensor = get(tgt, src);
  return tensor.get(tgt_ind % solution_dimension, src_ind % solution_dimension);
}


double Kernel::get(unsigned int tgt_ind, unsigned int src_ind) const {
  Dof tgt, src;
  unsigned int i_point_index = tgt_ind / solution_dimension;
  unsigned int i_points_vec_index = i_point_index * domain_dimension;
  unsigned int j_point_index = src_ind / solution_dimension;
  unsigned int j_points_vec_index = j_point_index * domain_dimension;

  tgt.point = Vec2(boundary->points[i_points_vec_index],
                   boundary->points[i_points_vec_index + 1]);
  tgt.normal = Vec2(boundary->normals[i_points_vec_index],
                    boundary->normals[i_points_vec_index + 1]);
  tgt.curvature = boundary->curvatures[i_point_index];
  tgt.weight = boundary->weights[i_point_index];
  tgt.is_boundary = true;

  src.point = Vec2(boundary->points[j_points_vec_index],
                   boundary->points[j_points_vec_index + 1]);
  src.normal = Vec2(boundary->normals[j_points_vec_index],
                    boundary->normals[j_points_vec_index + 1]);
  src.curvature = boundary->curvatures[j_point_index];
  src.weight = boundary->weights[j_point_index];
  src.is_boundary = true;

  ie_Mat tensor = get(tgt, src);
  return tensor.get(tgt_ind % solution_dimension, src_ind % solution_dimension);
}


ie_Mat Kernel::stokes_kernel(const Dof& tgt, const Dof& src) const {
  // double layer
  double scale = 1.0 / (M_PI);

  ie_Mat tensor(2, 2);
  ie_Mat normop(2, 2);
  if (tgt.is_boundary) {
    normop.set(0, 0, src.weight * tgt.normal.a[0] * src.normal.a[0]);
    normop.set(0, 1, src.weight * tgt.normal.a[0] * src.normal.a[1]);
    normop.set(1, 0, src.weight * tgt.normal.a[1] * src.normal.a[0]);
    normop.set(1, 1, src.weight * tgt.normal.a[1] * src.normal.a[1]);
  }
  if (tgt.point.a[0] == src.point.a[0] && tgt.point.a[1] == src.point.a[1]) {
    // tangent
    double t0 = -src.normal.a[1];
    double t1 =  src.normal.a[0];
    double potential = - 0.5 * src.curvature * src.weight * scale;
    tensor.set(0, 0, -0.5 + potential * t0 * t0);
    tensor.set(1, 1, -0.5 + potential * t1 * t1);
    tensor.set(0, 1, potential * t0 * t1);
    tensor.set(1, 0, potential * t1 * t0);

    tensor += normop;
    return tensor;
  }
  Vec2 r = tgt.point - src.point;
  double r0 = r.a[0];
  double r1 = r.a[1];
  double potential = src.weight * scale * (r.dot(src.normal)) / (pow(r.dot(
                       r), 2));
  tensor.set(0, 0, potential * r0 * r0);
  tensor.set(1, 1, potential * r1 * r1);
  tensor.set(0, 1, potential * r0 * r1);
  tensor.set(1, 0, potential * r1 * r0);
  tensor += normop;
  return tensor;
}


ie_Mat Kernel::laplace_kernel(const Dof& tgt, const Dof& src) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat tensor(1, 1);
  if (tgt.point.a[0] == src.point.a[0] && tgt.point.a[1] == src.point.a[1]) {
    tensor.set(0, 0, 0.5 + 0.5 * src.curvature * src.weight * scale);
    return tensor;
  }
  Vec2 r = tgt.point - src.point;
  tensor.set(0, 0, -src.weight * scale * (r.dot(src.normal)) / (r.dot(r)));
  return tensor;
}

// This function stores the DoF data,  and calculates the diagonals of the mat
void Kernel::load(Boundary* boundary_,
                  const std::vector<double>& domain_points_,
                  ie_solver_config::Pde pde_, int solution_dimension_,
                  int domain_dimension_) {
  solution_dimension = solution_dimension_;
  domain_points = domain_points_;
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


ie_Mat Kernel::forward_get(const std::vector<unsigned int>& I_,
                           const std::vector<unsigned int>& J_) const {
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (unsigned int i = 0; i < I_.size(); i++) {
    for (unsigned int j = 0; j < J_.size(); j++) {
      ret.mat[i + olda_ * j] = forward_get(I_[i], J_[j]);
    }
  }
  return ret;
}


ie_Mat Kernel::operator()(const std::vector<Dof>& tgts,
                          const std::vector<Dof>& srcs) const {
  ie_Mat ret(solution_dimension * tgts.size(),
             solution_dimension * srcs.size());
  for (unsigned int i = 0; i < tgts.size(); i++) {
    for (unsigned int j = 0; j < srcs.size(); j++) {
      ie_Mat tensor = get(tgts[i], srcs[j]);
      for (int k = 0; k < solution_dimension; k++) {
        for (int l = 0; l < solution_dimension; l++) {
          ret.set(solution_dimension * i + k,
                  solution_dimension * j + l, tensor.get(k, l));
        }
      }
    }
  }
  return ret;
}


}  // namespace ie_solver
