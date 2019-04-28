// Copyright John Paul Ryan 2019
#ifndef IE_SOLVER_IE_SOLVER_CONFIG_H_
#define IE_SOLVER_IE_SOLVER_CONFIG_H_

#include <fstream>
#include <memory>

#define DEFAULT_NUM_DISCRETIZATION_POINTS 1000
#define DEFAULT_ID_TOL 1e-6
#define DEFAULT_DOMAIN_SIZE 20

namespace ie_solver {
struct ie_solver_config {
  enum Pde {
    LAPLACE,
    STOKES
  };
  int num_boundary_points = DEFAULT_NUM_DISCRETIZATION_POINTS;
  int domain_size = DEFAULT_DOMAIN_SIZE;
  int domain_dimension = 2;
  int solution_dimension = 1;
  double id_tol = DEFAULT_ID_TOL;
  Pde pde = LAPLACE;
  bool is_strong_admissibility = false;
  Boundary::BoundaryCondition boundary_condition =
    Boundary::BoundaryCondition::SINGLE_ELECTRON;
  Boundary::BoundaryShape boundary_shape =
    Boundary::BoundaryShape::CIRCLE;
  bool scaling = false;
  bool animation = false;
  bool testing = false;
  std::ofstream n_scaling_output, error_scaling_output;
};

}  // namespace ie_solver
#endif  // IE_SOLVER_IE_SOLVER_CONFIG_H_
