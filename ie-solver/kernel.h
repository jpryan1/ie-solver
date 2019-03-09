// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_KERNEL_H_
#define IE_SOLVER_KERNEL_H_

#include <vector>
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/ie_mat.h"
#include "ie-solver/helpers.h"

namespace ie_solver {


struct Dof {
  Vec2 point;
  Vec2 normal;
  double curvature;
  double weight;
};  // struct


struct Kernel {
  double scale, diag_00, diag_01, diag_11;
  ie_solver_config::Pde pde;
  // TODO(John) don't have the kernel store the boundary
  Boundary* boundary;
  double electric_kernel(unsigned int i, unsigned int j) const;

  double get(unsigned int i, unsigned int j) const;

  ie_Mat stokes_kernel(const Dof& a, const Dof& b) const;
  double laplace_kernel(const Dof& a, const Dof& b) const;

  // This function stores the DoF data,  and calculates the diagonals of the mat
  void load(Boundary* boundary, ie_solver_config::Pde pde);
// TODO(John) shouldn't this->I have the underscore after it, not this arg?
  ie_Mat operator()(const std::vector<unsigned int>& I_,
                    const std::vector<unsigned int>& J_) const;
};  // struct

}  // namespace ie_solver

#endif  // IE_SOLVER_KERNEL_H_
