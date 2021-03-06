// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_INITIALIZATION_H_
#define IE_SOLVER_INITIALIZATION_H_

#include <vector>
#include "ie_solver/ie_mat.h"
#include "ie_solver/boundaries/boundary.h"
#include "ie_solver/kernel/kernel.h"

namespace ie_solver {

struct Initialization {
  Initialization() {}
  ~Initialization() {}

  static void InitializeDomainKernel(ie_Mat* K,
                                     const std::vector<double>& domain_points,
                                     const Kernel& kernel,
                                     int solution_dimension);
  static void LaplaceNeumann_InitializeDomainKernel(ie_Mat* K,
      const std::vector<double>& points,
      const std::vector<double>& normals,
      const std::vector<double>& weights,
      const std::vector<double>& domain_points,
      const Kernel& kernel);
  static void Stokes_InitializeDomainKernel(ie_Mat* K,
      const std::vector<double>& points,
      const std::vector<double>& normals,
      const std::vector<double>& weights,
      const std::vector<double>& domain_points,
      const Kernel& kernel);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_INITIALIZATION_H_
