// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_KERNEL_H_
#define IE_SOLVER_KERNEL_H_

#include <vector>
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/ie_mat.h"
#include "ie-solver/ie_solver_config.h"

namespace ie_solver {


struct Dof {
  Vec2 point;
  Vec2 normal;
  double curvature;
  double weight;
  bool is_boundary;
};  // struct


struct Kernel {
  double diag_00, diag_01, diag_11;
  int solution_dimension, domain_dimension;

  ie_solver_config::Pde pde;
  // TODO(John) don't have the kernel store the boundary
  Boundary* boundary;
  std::vector<double> domain_points;
  double get(unsigned int i, unsigned int j) const;
  double forward_get(unsigned int i, unsigned int j) const;
  ie_Mat get(const Dof& a, const Dof& b) const;

  ie_Mat stokes_kernel(const Dof& a, const Dof& b) const;
  ie_Mat laplace_kernel(const Dof& a, const Dof& b) const;

  void load(Boundary* boundary, const std::vector<double>& domain_points, ie_solver_config::Pde pde,
    int solution_dimension, int domain_dimension);

// TODO(John) shouldn't this->I have the underscore after it, not this arg?
  ie_Mat operator()(const std::vector<unsigned int>& I_,
                    const std::vector<unsigned int>& J_) const;
                    
  ie_Mat forward_get(const std::vector<unsigned int>& I_,
                          const std::vector<unsigned int>& J_) const;
                    
  ie_Mat operator()(const std::vector<Dof>& tgts,
                          const std::vector<Dof>& srcs) const;
};  // struct

}  // namespace ie_solver

#endif  // IE_SOLVER_KERNEL_H_
