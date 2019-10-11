// Copyright 2019 John Paul Ryan
#include <cmath>
#include <cassert>
#include <iostream>
#include "ie_solver/boundaries/boundary.h"

namespace ie_solver {


void CubicBoundary::get_cubics(const std::vector<double>& x0_points,
                               const std::vector<double>& x1_points,
                               std::vector<std::vector<double>>* x0_cubics,
                               std::vector<std::vector<double>>* x1_cubics) {
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
    std::vector<double>* t_vals) {
  ie_Mat companion(3, 3);
  companion.set(1, 0, 1.0);
  companion.set(2, 1, 1.0);
  companion.set(0, 2, -y_cubic[0] / y_cubic[3]);
  companion.set(1, 2, -y_cubic[1] / y_cubic[3]);
  companion.set(2, 2, -y_cubic[2] / y_cubic[3]);
  *t_vals = companion.real_eigenvalues();
}


void CubicBoundary::interpolate(int bc_index, bool is_interior,
                                int nodes_per_spline,
                                BoundaryCondition boundary_condition,
                                const std::vector<std::vector<double>>& x0_cubics,
                                const std::vector<std::vector<double>>& x1_cubics) {
  // Must fill points, normals, curvatures, weights.
  // Points = duh
  // Normals = tangent of points rotated 90 deg clockwise
  // Curvatures = (x'y'' - x''y') / (x'^2 + y'^2)^1.5
  // assert positive curvature when testing pl0x
  // Weights = estimate integral with uniform quadrature (quad_size)
  int start = bc_index;
  int num_spline_points = x0_cubics.size();
  int quad_size = 10;
  std::vector<double> quad(2 * quad_size * num_spline_points * nodes_per_spline);
  int quad_idx = 0;
  for (int i = 0; i < num_spline_points; i++) {
    std::vector<double> x_cubic = x0_cubics[i];
    std::vector<double> y_cubic = x1_cubics[i];
    for (int j = 0; j < nodes_per_spline; j++) {
      double t = j / (nodes_per_spline + 0.0);

      double x = x_cubic[0] + t * x_cubic[1] + pow(t, 2) * x_cubic[2]
                 + pow(t, 3) * x_cubic[3];
      double y = y_cubic[0] + t * y_cubic[1] + pow(t, 2) * y_cubic[2]
                 + pow(t, 3) * y_cubic[3];

      points.push_back(x);
      points.push_back(y);

      for (int qt = 0; qt < quad_size; qt++) {
        double t_quad = ((j * quad_size) + qt) / (0.0 + quad_size * nodes_per_spline);
        double x_quad = x_cubic[0] + t_ * x_cubic[1] + pow(t_, 2) * x_cubic[2]
                    + pow(t_, 3) * x_cubic[3];
        double y_quad = y_cubic[0] + t_ * y_cubic[1] + pow(t_, 2) * y_cubic[2]
                    + pow(t_, 3) * y_cubic[3];
        quad[quad_idx++] = x_quad;
        quad[quad_idx++] = y_quad;
      }

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
          boundary_values.set(bc_index, 0, log(sqrt(pow(x + 3, 2) + pow(y + 2,
                                               2))) / (2 * M_PI));
          break;
        case BoundaryCondition::ALL_ONES:
          boundary_values.set(bc_index, 0, 1.0);
          break;
        case BoundaryCondition::ALL_NEG_ONES:
          boundary_values.set(bc_index, 0, -1.);
          break;
        case BoundaryCondition::ALL_ZEROS:
          boundary_values.set(bc_index, 0, 0.);
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
          if (x < -0.99 || x > 1.99) {
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

  std::vector<double> distances(end - start);
  int d_idx = 0;
  for (unsigned int i = 0; i < distances.size(); i++) {
    double dist_sum = 0;
    for (int j = 0; j < quad_size; j++) {
      int idx = 2 * (quad_size * i + j);

      dist_sum += sqrt(pow(quad[idx] - quad[(idx + 2) % quad.size()], 2)
                       + pow(quad[idx + 1] - quad[(idx + 3) % quad.size()], 2));

    }
    distances[d_idx++] = dist_sum;
  }
  weights.push_back((distances[distances.size() - 1] + distances[0]) / 2.);
  for (int i = 1; i < distances.size(); i++) {
    weights.push_back((distances[i - 1] + distances[i]) / 2.);
  }
}


bool CubicBoundary::is_in_domain(const Vec2& a) {
  const double* v = a.a;

  int winding_number = 0;
  for (int i = 0; i < 2 * num_outer_nodes; i += 2) {
    double dist = sqrt(pow(v[0] - points[i], 2) + pow(v[1] - points[i + 1], 2));
    if (dist < 1e-2) {
      return false;
    }
    int next_i = (i + 2) % (2 * num_outer_nodes);
    if (points[i] > v[0]) {
      if (points[i + 1] < v[1] && points[next_i + 1] > v[1]) {
        winding_number++;
      } else if (points[i + 1] > v[1] && points[next_i + 1] < v[1]) {
        winding_number--;
      }
    }
  }
  if (winding_number % 2 == 0) {
    return false;
  }

  int node_idx = num_outer_nodes;

  for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
    Hole hole = holes[hole_idx];

    winding_number = 0;
    for (int i = 2 * node_idx; i < 2 * node_idx + 2 * hole.num_nodes; i += 2) {
      double dist = sqrt(pow(v[0] - points[i], 2) +
                         pow(v[1] - points[i + 1], 2));
      if (dist < 1e-2) {
        return false;
      }
      int next_i = 2 * node_idx +
                   ((i + 2 - (2 * node_idx)) % (2 * hole.num_nodes));
      if (points[i] > v[0]) {
        if (points[i + 1] < v[1] && points[next_i + 1] > v[1]) {
          winding_number++;
        } else if (points[i + 1] > v[1] && points[next_i + 1] < v[1]) {
          winding_number--;
        }
      }
    }
    if (winding_number % 2 == 1) {
      return false;
    }
    node_idx += hole.num_nodes;
  }


  return true;
}

}  // namespace ie_solver
