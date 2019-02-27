// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie-solver/kernel.h"

namespace ie_solver {

double Kernel::get(unsigned int i, unsigned int j) const {
  if (is_stokes) {
    return stokes_kernel(i, j);
  } else {
    return laplace_kernel(i, j);
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


double Kernel::stokes_kernel(unsigned int i, unsigned int j) const {
  // So this is much more awkwardly written than the function in
  // stokes_init.cpp that just writes the entire matrix at once, but the cost
  // of writing the whole matrix is just stupid. So we power through this
  // function
  // below commented is single layer
  // // We first need to ascertain which DoFs are relevant here.
  // int dof_i = i/2;
  // int dof_j = j/2;
  // // If the DoFs are the same, follow the singular procedure

  // if(dof_i==dof_j){
  //   if(i!=j) return scale*diag_01;
  //   return scale*diag_00;
  // }
  // // If they are different, figure out which of the three tensor components
  // // is desired, then return it.

  // Vec2 x(boundary->points[2*dof_i], boundary->points[2*dof_i+1]);
  // Vec2 y(boundary->points[2*dof_j], boundary->points[2*dof_j+1]);
  // Vec2 r = x-y;

  // double r0 = r.a[0];
  // double r1 = r.a[1];

  // if(i%2==0 && j%2==0){
  //   return boundary->weights[dof_j]*scale*(log(1.0/r.norm()) +
  //   (r0*r0/r.dot(r)));
  // }
  // else if(i%2==1 && j%2==1){
  //   return boundary->weights[dof_j]*scale*(log(1.0/r.norm())
  //   + (r1*r1/r.dot(r)));
  // }
  // else{
  //   return boundary->weights[dof_j]*scale*((r1*r0/r.dot(r)));
  // }

  int dof_i = i / 2;
  int dof_j = j / 2;

  if (dof_i == dof_j) {
    double t0 = -boundary->normals[2 * dof_i + 1];
    double t1 =  boundary->normals[2 * dof_i];

    if (i % 2 == 0 && j % 2 == 0) {
      return -0.5 -
             0.5 * boundary->curvatures[dof_i] * boundary->weights[dof_i]
             * scale * t0 * t0;
    } else if (i % 2 == 1 && j % 2 == 1) {
      return -0.5 -
             0.5 * boundary->curvatures[dof_i] * boundary->weights[dof_i]
             * scale * t1 * t1;
    } else {
      return -0.5 * boundary->curvatures[dof_i] * boundary->weights[dof_i]
             * scale * t0 * t1;
    }
  }

  Vec2 x(boundary->points[2 * dof_i], boundary->points[2 * dof_i + 1]);
  Vec2 y(boundary->points[2 * dof_j], boundary->points[2 * dof_j + 1]);
  Vec2 r = x - y;

  double r0 = r.a[0];
  double r1 = r.a[1];
  Vec2 n(boundary->normals[2 * dof_j], boundary->normals[2 * dof_j + 1]);
  double potential = boundary->weights[dof_j] * scale * (r.dot(n)) / (pow(r.dot(
                       r), 2));

  if (i % 2 == 0 && j % 2 == 0) {
    return potential * r0 * r0;
  } else if (i % 2 == 1 && j % 2 == 1) {
    return  potential * r1 * r1;
  } else {
    return potential * r0 * r1;
  }
}

double Kernel::laplace_kernel(unsigned int i, unsigned int j) const {
  if (i == j) {
    return 0.5 + 0.5 * boundary->curvatures[i] * boundary->weights[i] * scale;
  }

  Vec2 x(boundary->points[2 * i], boundary->points[2 * i + 1]);
  Vec2 y(boundary->points[2 * j], boundary->points[2 * j + 1]);
  Vec2 r = x - y;

  Vec2 n(boundary->normals[2 * j], boundary->normals[2 * j + 1]);
  double potential = -boundary->weights[j] * scale * (r.dot(n)) / (r.dot(r));
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
