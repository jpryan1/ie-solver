// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_KERNEL_H_
#define IE_SOLVER_KERNEL_H_

#include <vector>
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/ie_mat.h"

namespace ie_solver {

struct Kernel {
  double scale, diag_00, diag_01, diag_11;
  bool is_stokes;
  // TODO(John) don't have the kernel store the boundary
  Boundary* boundary;
  double electric_kernel(unsigned int i, unsigned int j) const;

  double get(unsigned int i, unsigned int j) const;

  double stokes_kernel(unsigned int i, unsigned int j) const;
  double laplace_kernel(unsigned int i, unsigned int j) const;

  // This function stores the DoF data,  and calculates the diagonals of the mat
  void load(Boundary* boundary, bool is_stokes_);
// TODO(John) shouldn't this->I have the underscore after it, not this arg?
  ie_Mat operator()(const std::vector<unsigned int>& I_,
                    const std::vector<unsigned int>& J_) const;
};  // struct

}  // namespace ie_solver

#endif  // IE_SOLVER_KERNEL_H_
