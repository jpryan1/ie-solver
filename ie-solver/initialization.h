// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_INITIALIZATION_H_
#define IE_SOLVER_INITIALIZATION_H_

#include <vector>
#include "ie-solver/ie_mat.h"
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/kernel.h"

namespace ie_solver {

struct Initialization {
  Initialization() {}
  ~Initialization() {}

  void InitializeDomainKernel(ie_Mat* K,
                              const std::vector<double>& domain_points,
                              int test_size, Kernel* kernel,
                              bool is_stokes);

  void Stokes_InitializeDomainKernel(ie_Mat* K,
                                     const std::vector<double>& points,
                                     const std::vector<double>& normals,
                                     const std::vector<double>& weights,
                                     const std::vector<double>& domain_points,
                                     int test_size, Kernel* kernel);

  void Electric_InitializeBoundary(ie_Mat* f,
                                   const std::vector<double>& points);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_INITIALIZATION_H_
