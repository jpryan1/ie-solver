// Copyright John Paul Ryan 2019
#ifndef IE_SOLVER_IE_SOLVER_CONFIG_H_
#define IE_SOLVER_IE_SOLVER_CONFIG_H_

#include <fstream>
#include <memory>

#define DEFAULT_NUM_DISCRETIZATION_POINTS 1000
#define DEFAULT_ID_TOL 1e-6

namespace ie_solver {
struct ie_solver_config {
  enum Pde {
    LAPLACE,
    STOKES
  };
  enum Admissibility {
    WEAK,
    STRONG
  };
  int num_boundary_points = DEFAULT_NUM_DISCRETIZATION_POINTS;
  double id_tol = DEFAULT_ID_TOL;
  Pde pde = LAPLACE;
  Admissibility admissibility = WEAK;
  std::unique_ptr<Boundary> boundary;
  Boundary::BoundaryCondition boundary_condition =
    Boundary::BoundaryCondition::SINGLE_ELECTRON;
  bool scaling = false;
  bool testing = false;
  std::ofstream n_scaling_output, error_scaling_output;
};

}  // namespace ie_solver
#endif  // IE_SOLVER_IE_SOLVER_CONFIG_H_
