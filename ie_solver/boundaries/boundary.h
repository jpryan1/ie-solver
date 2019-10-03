// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_BOUNDARIES_BOUNDARY_H_
#define IE_SOLVER_BOUNDARIES_BOUNDARY_H_

#include <vector>
#include "ie_solver/vec2.h"
#include "ie_solver/ie_mat.h"

namespace ie_solver {

struct Hole {
  Vec2 center;
  double radius;
  int num_nodes;
};

enum BoundaryCondition {
  SINGLE_ELECTRON,
  ALL_ONES,
  ALL_HALFS,
  ALL_ZEROS,
  BUMP_FUNCTION,
  STOKES,
  TANGENT_VEC,
  REVERSE_TANGENT_VEC,
  NORMAL_VEC,
  REVERSE_NORMAL_VEC,
  LEFT_TO_RIGHT_FLOW,
  NO_SLIP
};

class Boundary {
 public:
  std::vector<double> perturbation_parameters;
  std::vector<double> points, normals, curvatures, weights;
  std::vector<Hole> holes;
  ie_Mat boundary_values;
  enum BoundaryShape {
    CIRCLE,
    ROUNDED_SQUARE,
    ROUNDED_SQUARE_WITH_BUMP,
    SQUIGGLY,
    ANNULUS,
    CUBIC_SPLINE,
    EX1,
    EX2,
    EX3
  };

  BoundaryShape boundary_shape;
  BoundaryCondition boundary_condition;
  virtual void initialize(int n, BoundaryCondition bc) = 0;
  virtual bool is_in_domain(const Vec2& a) = 0;
};

class CubicBoundary : public Boundary {
 public:
  // Note, the outer nodes must appear first in the point data vectors.
  int num_outer_nodes;

  virtual void get_spline_points(std::vector<double>* outer_x0_spline_points,
                                 std::vector<double>* outer_x1_spline_points) = 0;

  void get_cubics(const std::vector<double>& x0_points,
                  const std::vector<double>& x1_points,
                  std::vector<std::vector<double>>* x0_cubics,
                  std::vector<std::vector<double>>* x1_cubics);

  void interpolate(int bc_index, bool is_interior,
                   int nodes_per_spline, BoundaryCondition boundary_condition,
                   const std::vector<std::vector<double>>& x0_cubics,
                   const std::vector<std::vector<double>>& x1_cubics);

  void find_real_roots_of_cubic(const std::vector<double>& y_cubic,
                                std::vector<double>* t_vals);
  int num_right_intersections(double x, double y, int index);
  bool is_in_domain(const Vec2& a);

};

}  // namespace ie_solver

#endif  // IE_SOLVER_BOUNDARIES_BOUNDARY_H_
