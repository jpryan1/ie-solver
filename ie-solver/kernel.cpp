// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie-solver/kernel.h"

namespace ie_solver {

double Kernel::get(unsigned int i, unsigned int j) const {
  Dof a, b;
  if (is_stokes) {
    unsigned int dof_i = i/2;
    unsigned int dof_j = j/2;
    a.point = Vec2(boundary->points[2 * dof_i], boundary->points[2 * dof_i + 1]);
    a.normal = Vec2(boundary->normals[2 * dof_i], boundary->normals[2 * dof_i + 1]);
    a.curvature = boundary->curvatures[dof_i];
    a.weight = boundary->weights[dof_i];

    b.point = Vec2(boundary->points[2 * dof_j], boundary->points[2 * dof_j + 1]);
    b.normal = Vec2(boundary->normals[2 * dof_j], boundary->normals[2 * dof_j + 1]);
    b.curvature = boundary->curvatures[dof_j];
    b.weight = boundary->weights[dof_j];
    ie_Mat tensor = stokes_kernel(a,b);
    if(i%2 == 0){
      if(j%2==0){
        return tensor.get(0,0);
      }else{
        return tensor.get(0,1);
      }
    }else{
      if(j%2==0){
        return tensor.get(1,0);
      }else{
        return tensor.get(0,1);
      }
    }
  } else {
    a.point = Vec2(boundary->points[2 * i], boundary->points[2 * i + 1]);
    a.normal = Vec2(boundary->normals[2 * i], boundary->normals[2 * i + 1]);
    a.curvature = boundary->curvatures[i];
    a.weight = boundary->weights[i];

    b.point = Vec2(boundary->points[2 * j], boundary->points[2 * j + 1]);
    b.normal = Vec2(boundary->normals[2 * j], boundary->normals[2 * j + 1]);
    b.curvature = boundary->curvatures[j];
    b.weight = boundary->weights[j];

    return laplace_kernel(a, b);
    // return electric_kernel(i,j);
  }
}


double Kernel::electric_kernel(unsigned int i, unsigned int j) const {
  if (i == j) {
    return 2;
  }

  Vec2 x(boundary->points[2 * i], boundary->points[2 * i + 1]);
  Vec2 y(boundary->points[2 * j], boundary->points[2 * j + 1]);
  Vec2 r = x - y;


  double potential = log(r.norm());

  return potential;
}


ie_Mat Kernel::stokes_kernel(const Dof& a, const Dof& b) const {
 
  // double layer
  ie_Mat tensor(2,2);

  if (a.point.a[0] == b.point.a[0] && a.point.a[1] == b.point.a[1]) {
    
    // tangent
    double t0 = -a.normal.a[1];
    double t1 =  a.normal.a[0];
    double potential = -0.5 - 0.5 * a.curvature * a.weight * scale;
    tensor.set(0,0, potential * t0 * t0);
    tensor.set(1,1, potential * t1 * t1);
    tensor.set(0,1, potential * t0 * t1);
    tensor.set(1,0, potential * t1 * t0);
    return tensor;
  }

  Vec2 r = a.point - b.point;

  double r0 = r.a[0];
  double r1 = r.a[1];
  double potential = b.weight * scale * (r.dot(b.normal)) / (pow(r.dot(
                       r), 2));
  
  tensor.set(0,0, potential * r0 * r0);
  tensor.set(1,1, potential * r1 * r1);
  tensor.set(0,1, potential * r0 * r1);
  tensor.set(1,0, potential * r1 * r0);
  return tensor;
}

double Kernel::laplace_kernel(const Dof& a, const Dof& b) const {
  if (a.point.a[0] == b.point.a[0] && a.point.a[1] == b.point.a[1]) {
    return 0.5 + 0.5 * a.curvature * a.weight * scale;
  }

  Vec2 r = a.point - b.point;

  double potential = -b.weight * scale * (r.dot(b.normal)) / (r.dot(r));
  return potential;
}

// This function stores the DoF data,  and calculates the diagonals of the mat
void Kernel::load(Boundary* boundary_, bool is_stokes_) {
  boundary = boundary_;
  is_stokes = is_stokes_;
  if (is_stokes) {
    scale = 1 / (M_PI);
  } else {
    scale = 1 / (2 * M_PI);
  }

  if (is_stokes) {
    double avg = 0;
    for (double weight : boundary->weights) {
      avg += weight;
    }
    avg /= boundary->weights.size();

    double alpha = avg / 2.0;
    double beta  = alpha * alpha;
    diag_00 = 8 * beta + 2 * beta * log(1 / (2 * beta)) - beta * M_PI;
    diag_01 = 0;
    diag_11 = diag_00;

    // printf("Diagonals are %f\n", diag_00);
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
