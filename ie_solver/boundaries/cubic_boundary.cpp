// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {
  
     
void CubicBoundary::get_cubics(const std::vector<double>& x0_points,
                const std::vector<double>& x1_points,
                std::vector<std::vector<double>>* x0_cubics,
                std::vector<std::vector<double>>* x1_cubics){
  int num_spline_points = x0_points.size();
  ie_Mat x0_d_vec(num_spline_points, 1);
  ie_Mat x0_y_vec(num_spline_points, 1);

  ie_Mat x1_d_vec(num_spline_points, 1);
  ie_Mat x1_y_vec(num_spline_points, 1);
  for (int i = 0; i < num_spline_points; i++) {
    double after = x0_points[(i + 1) % num_spline_points];
    double before = x0_points[(i + num_spline_points - 1) %
                              num_spline_points];
    x0_y_vec.set(i, 0,  3 * (after - before));

    after = x1_points[(i + 1) % num_spline_points];
    before = x1_points[(i + num_spline_points - 1) % num_spline_points];
    x1_y_vec.set(i, 0,  3 * (after - before));
  }

  ie_Mat A(num_spline_points, num_spline_points);

  A.set(0, num_spline_points - 1, 1.0);
  A.set(num_spline_points - 1, 0, 1.0);

  for (int i = 0; i < num_spline_points; i++) {
    A.set(i, i, 4.0);
    if (i != 0) {
      A.set(i, i - 1, 1.0);
    }
    if (i != num_spline_points - 1) {
      A.set(i, i + 1, 1.0);
    }
  }

  A.left_multiply_inverse(x0_y_vec, &x0_d_vec);
  A.left_multiply_inverse(x1_y_vec, &x1_d_vec);


  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> cubic;
    cubic.push_back(x0_points[i]);
    cubic.push_back(x0_d_vec.get(i, 0));
    int i_plus_one = (i + 1) % num_spline_points;
    cubic.push_back(3 * (x0_points[i_plus_one] - x0_points[i])
                    - 2 * x0_d_vec.get(i, 0)
                    - x0_d_vec.get(i_plus_one, 0));
    cubic.push_back(2 * (x0_points[i] - x0_points[i_plus_one])
                    + x0_d_vec.get(i, 0)
                    + x0_d_vec.get(i_plus_one, 0));

    x0_cubics->push_back(cubic);
  }
  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> cubic;
    cubic.push_back(x1_points[i]);
    cubic.push_back(x1_d_vec.get(i, 0));
    int i_plus_one = (i + 1) % num_spline_points;
    cubic.push_back(3 * (x1_points[i_plus_one] - x1_points[i])
                    - 2 * x1_d_vec.get(i, 0)
                    - x1_d_vec.get(i_plus_one, 0));
    cubic.push_back(2 * (x1_points[i] - x1_points[i_plus_one])
                    + x1_d_vec.get(i, 0)
                    + x1_d_vec.get(i_plus_one, 0));
    x1_cubics->push_back(cubic);
  }
}


void CubicBoundary::find_real_roots_of_cubic(const std::vector<double>& y_cubic,
                              std::vector<double>* t_vals){
  ie_Mat companion(3, 3);
  companion.set(1, 0, 1.0);
  companion.set(2, 1, 1.0);
  companion.set(0, 2, -y_cubic[0] / y_cubic[3]);
  companion.set(1, 2, -y_cubic[1] / y_cubic[3]);
  companion.set(2, 2, -y_cubic[2] / y_cubic[3]);
  *t_vals = companion.real_eigenvalues();
}


int CubicBoundary::num_right_intersections(double x, double y, int index){
  // Find intersection on spline, if no t in [0,1] return 0
  // For each t in that range, find corresponding x val, check if
  // greater than x.

  std::vector<double> x_cubic = all_cubics_x0[index];
  std::vector<double> y_cubic = all_cubics_x1[index];
  y_cubic[0] -= y;
  std::vector<double> t_vals;

  find_real_roots_of_cubic(y_cubic, &t_vals);
  int intersections = 0;
  for (double t : t_vals) {
    if (t > 0 && t < 1) {
      double dif = x_cubic[0] + t * x_cubic[1] + pow(t, 2) * x_cubic[2]
                   + pow(t, 3) * x_cubic[3] - x;
      if (fabs(dif) < 0.05) return -1;
      if (dif > 0) {
        intersections++;
      }
    }
  }
  return intersections;
}


void CubicBoundary::interpolate(int bc_index, bool is_interior, int nodes_per_spline,
                   BoundaryCondition boundary_condition,
                   const std::vector<std::vector<double>>& x0_cubics,
                   const std::vector<std::vector<double>>& x1_cubics)
 {
  // Must fill points, normals, curvatures, weights.
  // Points = duh
  // Normals = tangent of points rotated 90 deg clockwise
  // Curvatures = (x'y'' - x''y') / (x'^2 + y'^2)^1.5
  // assert positive curvature when testing pl0x
  // Weights = uggggh just estimate it
  int start = bc_index;
  int num_spline_points = x0_cubics.size();
  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> x_cubic = x0_cubics[i];
    std::vector<double> y_cubic = x1_cubics[i];
    for (int j = 0; j < nodes_per_spline; j++) {
      double t = j / (nodes_per_spline+0.0);
      double x = x_cubic[0] + t * x_cubic[1] + pow(t, 2) * x_cubic[2]
                 + pow(t, 3) * x_cubic[3];
      double y = y_cubic[0] + t * y_cubic[1] + pow(t, 2) * y_cubic[2]
                 + pow(t, 3) * y_cubic[3];

      points.push_back(x);
      points.push_back(y);

      double x_prime = x_cubic[1] + 2 * t * x_cubic[2]
                       + 3 * pow(t, 2) * x_cubic[3];
      double y_prime = y_cubic[1] + 2 * t * y_cubic[2]
                       + 3 * pow(t, 2) * y_cubic[3];

      double x_prime_prime = 2 * x_cubic[2] + 6 * t * x_cubic[3];
      double y_prime_prime = 2 * y_cubic[2] + 6 * t * y_cubic[3];

      double curvature = (x_prime * y_prime_prime - x_prime_prime * y_prime)
                         / pow(sqrt(pow(x_prime, 2) + pow(y_prime, 2)), 3);
      double norm = sqrt(pow(x_prime, 2) + pow(y_prime, 2));
      x_prime /= norm;
      y_prime /= norm;
      //BUG have observed x_prime = 0 ?
      if (is_interior) {
        curvatures.push_back(-curvature);
        normals.push_back(-y_prime);
        normals.push_back(x_prime);
      } else {
        curvatures.push_back(curvature);
        normals.push_back(y_prime);
        normals.push_back(-x_prime);
      }

      switch (boundary_condition) {
        case BoundaryCondition::SINGLE_ELECTRON:
          boundary_values.set(bc_index, 0, log(sqrt(pow(x + 2, 2) + pow(y + 2, 2))) / (2 * M_PI));
          break;
        case BoundaryCondition::ALL_ONES:
          boundary_values.set(bc_index, 0, 1.0);
          break;
        case BoundaryCondition::BUMP_FUNCTION: {
          double N = boundary_values.height();
          // This x_val is -1 at i=0 and +1 at i=N-1
          double x_val = ((2.0 * i + 1.0 - N) / (N - 1.0));
          boundary_values.set(bc_index, 0, exp(-1.0 / (1.0 - pow(x_val, 2))));
          break;
        }
        case BoundaryCondition::STOKES:
          boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index + 1]);
          boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index]);
          break;
        case TANGENT_VEC:
          boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index + 1]);
          boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index]);
          break;
        case REVERSE_TANGENT_VEC:
          boundary_values.set(2 * bc_index, 0, normals[2 * bc_index + 1]);
          boundary_values.set(2 * bc_index + 1, 0, -normals[2 * bc_index]);
          break;
        case NORMAL_VEC:
          boundary_values.set(2 * bc_index, 0, normals[2 * bc_index ]);
          boundary_values.set(2 * bc_index + 1, 0, normals[2 * bc_index + 1]);
          break;
        case REVERSE_NORMAL_VEC:
          boundary_values.set(2 * bc_index, 0, -normals[2 * bc_index ]);
          boundary_values.set(2 * bc_index + 1, 0, -normals[2 * bc_index + 1]);
          break;
        case LEFT_TO_RIGHT_FLOW:
          if (x < -1.01 || x > 1.99) {
            boundary_values.set(2 * bc_index, 0, 1);
            boundary_values.set(2 * bc_index + 1, 0, 0);
          } else {
            boundary_values.set(2 * bc_index, 0, 0);
            boundary_values.set(2 * bc_index + 1, 0, 0);
          }
          break;
        case NO_SLIP:
          boundary_values.set(2 * bc_index, 0, 0.);
          boundary_values.set(2 * bc_index + 1, 0, 0.);
          break;
      }
      bc_index++;
    }
  }
  int end = bc_index;

 double dist1 = sqrt(pow(points[2 * start]
                          - points[2 * (end - 1)], 2)
                      + pow(points[2 * start + 1]
                            - points[2 * (end - 1) + 1], 2));

  double dist2 = sqrt(pow(points[2 * start]
                          - points[2 * (start + 1)], 2)
                      + pow(points[2 * start + 1]
                            - points[2 * (start + 1) + 1], 2));
  weights.push_back((dist1 + dist2) / 2.0);
  for (unsigned int i = start + 1; i < end - 1; i++) {
    dist1 = sqrt(pow(points[2 * i] - points[2 * (i - 1)], 2)
                 + pow(points[2 * i + 1] - points[2 * (i - 1) + 1], 2));
    dist2 = sqrt(pow(points[2 * i] - points[2 * (i + 1)], 2)
                 + pow(points[2 * i + 1] - points[2 * (i + 1) + 1], 2));
    weights.push_back((dist1 + dist2) / 2.0);
  }

  dist1 = sqrt(pow(points[2 * (end - 1)]
                   - points[2 * ((end - 1) - 1)], 2)
               + pow(points[2 * (end - 1) + 1]
                     - points[2 * ((end - 1) - 1) + 1], 2));
  dist2 = sqrt(pow(points[2 * (end - 1)]
                   - points[2 * (start)], 2)
               + pow(points[2 * (end - 1) + 1]
                     - points[2 * (start) + 1], 2));
  weights.push_back((dist1 + dist2) / 2.0);
}



}  // namespace ie_solver