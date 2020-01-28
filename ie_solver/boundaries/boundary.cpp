// Copyright 2019 John Paul Ryan
#include <cmath>
#include <cassert>
#include <iostream>
#include "ie_solver/boundaries/boundary.h"


namespace ie_solver {


void Boundary::set_boundary_values_size(BoundaryCondition bc) {
  int num_points = weights.size();
  switch (bc) {
    case SINGLE_ELECTRON:
    case EX3A:
    case ALL_ONES:
    case ALL_NEG_ONES:
    case ALL_ZEROS:
      boundary_values = ie_Mat(num_points, 1);
      break;
    // Everything above is for 1D solutions.
    // Everything below is for 2D solutions.
    case TANGENT_VEC:
    case REVERSE_TANGENT_VEC:
    case NORMAL_VEC:
    case REVERSE_NORMAL_VEC:
    case LEFT_TO_RIGHT_FLOW:
    case EX3B:
    case HORIZONTAL_VEC:
    case NO_SLIP:
      boundary_values = ie_Mat(2 * num_points, 1);
      break;
    case DEFAULT: {
      std::cout << "DEFAULT BoundaryCondition enum not to be given " <<
                "directly to set_boundary_value_size function." << std::endl;
      exit(-1);
      break;
    }
  }
}


void Boundary::apply_boundary_condition(int start_point_idx, int end_point_idx,
                                        BoundaryCondition bc) {
  for (int point_idx = start_point_idx;
       point_idx < end_point_idx;
       point_idx++) {
    switch (bc) {
      case SINGLE_ELECTRON: {
        boundary_values.set(point_idx, 0,
                            log(sqrt(pow(points[2 * point_idx] + 3, 2)
                                     + pow(points[2 * point_idx + 1] + 2, 2)))
                            / (2 * M_PI));
        break;
      }
      case ALL_ONES: {
        boundary_values.set(point_idx, 0, 1.);
        break;
      }
      case ALL_NEG_ONES: {
        boundary_values.set(point_idx, 0, -1.);
        break;
      }
      case ALL_ZEROS: {
        boundary_values.set(point_idx, 0, 0.);
        break;
      }
      // Everything above is for 1D solutions.
      // Everything below is for 2D solutions.
      case TANGENT_VEC: {
        boundary_values.set(2 * point_idx, 0, -normals[2 * point_idx + 1]);
        boundary_values.set(2 * point_idx + 1, 0, normals[2 * point_idx]);
        break;
      }
      case REVERSE_TANGENT_VEC: {
        boundary_values.set(2 * point_idx, 0, normals[2 * point_idx + 1]);
        boundary_values.set(2 * point_idx + 1, 0, -normals[2 * point_idx]);
        break;
      }
      case NORMAL_VEC: {
        boundary_values.set(2 * point_idx, 0, normals[2 * point_idx ]);
        boundary_values.set(2 * point_idx + 1, 0, normals[2 * point_idx + 1]);
        break;
      }
      case REVERSE_NORMAL_VEC: {
        boundary_values.set(2 * point_idx, 0, -normals[2 * point_idx ]);
        boundary_values.set(2 * point_idx + 1, 0, -normals[2 * point_idx + 1]);
        break;
      } case HORIZONTAL_VEC: {
        boundary_values.set(2 * point_idx, 0, 1.);
        boundary_values.set(2 * point_idx + 1, 0, 0.);
        break;
      }
      case LEFT_TO_RIGHT_FLOW: {
        double x = points[2 * point_idx];
        if (x < -0.99 || x > 1.99) {
          boundary_values.set(2 * point_idx, 0, 1);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        } else {
          boundary_values.set(2 * point_idx, 0, 0);
          boundary_values.set(2 * point_idx + 1, 0, 0);
        }
        break;
      }
      case NO_SLIP: {
        boundary_values.set(2 * point_idx, 0, 0.);
        boundary_values.set(2 * point_idx + 1, 0, 0.);
        break;
      }
      case EX3A:
      case EX3B:
      case DEFAULT: {
        std::cout << "BoundaryCondition enum not to be given " <<
                  "directly to apply_boundary_condition function." << std::endl;
        exit(-1);
        break;
      }
    }
  }
}

}  // namespace ie_solver
