// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie-solver/kernel.h"

namespace ie_solver {
  
ie_Mat Kernel::get(const Dof& a, const Dof& b) const{
   switch(pde){
    case ie_solver_config::Pde::LAPLACE:
      return laplace_kernel(a, b);
      break;
    case ie_solver_config::Pde::STOKES:
      return stokes_kernel(a,b);
      break;
   }
}


double Kernel::get(unsigned int i, unsigned int j) const {
  Dof a, b;
  // TODO(John) the switch should be as small as the above, a and b should
  // be defined with solution_dimension and domain_dimension
  switch(pde){
    
    case ie_solver_config::Pde::LAPLACE:
      a.point = Vec2(boundary->points[2 * i], boundary->points[2 * i + 1]);
      a.normal = Vec2(boundary->normals[2 * i], boundary->normals[2 * i + 1]);
      a.curvature = boundary->curvatures[i];
      a.weight = boundary->weights[i];
  
      b.point = Vec2(boundary->points[2 * j], boundary->points[2 * j + 1]);
      b.normal = Vec2(boundary->normals[2 * j], boundary->normals[2 * j + 1]);
      b.curvature = boundary->curvatures[j];
      b.weight = boundary->weights[j];
      return laplace_kernel(a, b).get(0,0);
      break;
    case ie_solver_config::Pde::STOKES:
      unsigned int dof_i = i / 2;
      unsigned int dof_j = j / 2;
      a.point = Vec2(boundary->points[2 * dof_i],
                     boundary->points[2 * dof_i + 1]);
      a.normal = Vec2(boundary->normals[2 * dof_i],
                      boundary->normals[2 * dof_i + 1]);
      a.curvature = boundary->curvatures[dof_i];
      a.weight = boundary->weights[dof_i];
  
      b.point = Vec2(boundary->points[2 * dof_j],
                     boundary->points[2 * dof_j + 1]);
      b.normal = Vec2(boundary->normals[2 * dof_j],
                      boundary->normals[2 * dof_j + 1]);
      b.curvature = boundary->curvatures[dof_j];
      b.weight = boundary->weights[dof_j];
      ie_Mat tensor = stokes_kernel(a, b);
      if (i % 2 == 0) {
        if (j % 2 == 0) {
          return tensor.get(0, 0);
        } else {
          return tensor.get(0, 1);
        }
      } else {
        if (j % 2 == 0) {
          return tensor.get(1, 0);
        } else {
          return tensor.get(1, 1);
        }
      }
      break;
  }
}


ie_Mat Kernel::stokes_kernel(const Dof& a, const Dof& b) const {
  // double layer
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
  ie_Mat tensor(1,1);
  if (a.point.a[0] == b.point.a[0] && a.point.a[1] == b.point.a[1]) {
    tensor.set(0,0, 0.5 + 0.5 * a.curvature * a.weight * scale);
    return tensor;
  }
  Vec2 r = a.point - b.point;
  tensor.set(0, 0, -b.weight * scale * (r.dot(b.normal)) / (r.dot(r)));
  return tensor;
}

// This function stores the DoF data,  and calculates the diagonals of the mat
void Kernel::load(Boundary* boundary_, ie_solver_config::Pde pde_) {
  boundary = boundary_;
  pde = pde_;
  if (pde == ie_solver_config::Pde::STOKES) { //MOVE THIS TO KERNEL FUNCTIONS
    scale = 1 / (M_PI);
  } else {
    scale = 1 / (2 * M_PI);
  }
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
