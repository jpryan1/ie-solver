// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_IO_IO_H_
#define IE_SOLVER_IO_IO_H_

#include <string>
#include <vector>
#include <fstream>
#include "ie_solver/ie_mat.h"
#include "ie_solver/quadtree/quadtree.h"

#define DEFAULT_NUM_DISCRETIZATION_POINTS 1000
#define DEFAULT_ID_TOL 1e-6
#define DEFAULT_DOMAIN_SIZE 20

namespace ie_solver {

struct ie_solver_config {
  enum Pde {
    LAPLACE,
    LAPLACE_NEUMANN,
    STOKES
  };
  int num_boundary_points = DEFAULT_NUM_DISCRETIZATION_POINTS;
  int domain_size = DEFAULT_DOMAIN_SIZE;
  int domain_dimension = 2;
  int solution_dimension = 1;
  double id_tol = DEFAULT_ID_TOL;
  Pde pde = LAPLACE;
  bool is_strong_admissibility = false;
  BoundaryCondition boundary_condition =
    BoundaryCondition::SINGLE_ELECTRON;
  Boundary::BoundaryShape boundary_shape =
    Boundary::BoundaryShape::CIRCLE;
  bool scaling = false;
  bool animation = false;
  bool testing = false;
  std::ofstream n_scaling_output, error_scaling_output;
};

struct io {
  static void write_boundary_to_file(const std::string& filename,
                                     const std::vector<double>& points);
  static void write_potential_to_file();
  static void write_times_to_files(int* scale_n,
                                   const std::vector<double>& n_times,
                                   double* scale_eps,
                                   const std::vector<double>& eps_times);
  static void write_solution_to_file(const std::string& filename,
                                     const ie_Mat& domain,
                                     const std::vector<double>& domain_points,
                                     int solution_dimension);
  static void write_quadtree_to_file(const std::string& filename,
                                     const QuadTree& quadtree);

  static void write_ex2_gradients_to_file(const std::string& filename,
                                          std::vector<double> angs_and_grads);
  static void write_ex3_flows_to_file(const std::string& filename,
                                      std::vector<double> ang_and_flow);
  static int parse_input_into_config(int argc, char** argv,
                                     ie_solver_config* config);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_IO_IO_H_
