// Copyright 2019 John Paul Ryan
#include "ie_solver/boundaries/squiggly.h"
#include <vector>
#include <boost/math/special_functions/ellint_2.hpp>

namespace ie_solver {

void Squiggly::draw_squiggle(int bc_index, int num_points, double start_x,
                             double start_y, double end_x, double end_y) {
  // By default, the normal vector is on the left, and the amplitude is 1.
  // The strategy here is going to be to draw the squiggle on [0,3pi], then
  // appropriately scale, rotate, and translate

  double x_diff = end_x - start_x;
  double y_diff = end_y - start_y;
  double angle = atan2(y_diff, x_diff);
  double scale = sqrt(pow(x_diff, 2) + pow(y_diff, 2)) / (3.0 * M_PI);
  for (int i = 0; i < num_points; i++) {
    double t = (3 * M_PI * i) / num_points;
    double x = t;
    double y = sin(t);

    double x_prime = 1;
    double y_prime = cos(t);

    x *= scale;
    x_prime *= scale;
    double prime_norm = sqrt(pow(x_prime, 2) + pow(y_prime, 2));
    double norm_x = -y_prime / prime_norm;
    double norm_y = x_prime / prime_norm;

    // curvature = (x'y''-y'x'') / |(x',y')|^3
    double curvature = x_prime * (-sin(t)) / pow(x_prime +
                       pow(cos(t), 2) , 1.5);
    curvatures.push_back(curvature);

    double transformed_x = cos(angle) * x - sin(angle) * y + start_x;
    double transformed_y = sin(angle) * x + cos(angle) * y + start_y;
    points.push_back(transformed_x);
    points.push_back(transformed_y);

    double transformed_x_norm = cos(angle) * norm_x - sin(angle) * norm_y;
    double transformed_y_norm = sin(angle) * norm_x + cos(angle) * norm_y;
    normals.push_back(transformed_x_norm);
    normals.push_back(transformed_y_norm);
  }


  // now we need to figure out some arclengths
  double k = 1.0 / sqrt(2.0);
  int side_scale = num_points / 6;
  std::vector<double> E(side_scale + 1);
  E[0] = 0.0;  // this isn't used, just to make the indexing a bit easier
  for (int i = 1; i < side_scale + 1; i++) {
    double ang = 0.5 * M_PI * ((i + 0.0) / side_scale);
    E[i] = scale * sqrt(2.0) * boost::math::ellint_2(k, ang);
  }

  for (int k = 0; k < 3; k++) {
    weights.push_back(E[1]);

    for (int i = 1; i < side_scale; i++) {
      weights.push_back((E[i + 1] - E[i - 1]) / 2.0);
    }


    weights.push_back(E[side_scale] - E[side_scale - 1]);
    for (int i = side_scale - 1; i >= 1; i--) {
      weights.push_back((E[i + 1] - E[i - 1]) / 2.0);
    }
  }
}


void Squiggly::initialize(int N, BoundaryCondition bc) {
  boundary_shape = BoundaryShape::SQUIGGLY;
  boundary_condition = bc;
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  //            SHAPE PROTOCOL
  //
  // Except in trivial cases, we need to be careful if we want to evenly
  // distribute discretization points on the boundary, while compartmentalizing
  // drawing functions. For that purpose, we establish a SCALE_UNIT and require
  // that all draw functions take some multiple of this unit. For example, if
  // the ratio of shape A's points to shape B's points should be x/y, then shape
  // A can have x SCALE_UNIT's and shape B can have y.
  //
  //    Squiggly
  //
  //    Shape       Num of SCALE_UNIT's     Num of Shape in Boundary
  //    side        6                       4
  // __________________________________________
  //
  //    Num of SCALE_UNIT's = 6*4 = 24
  N = 24 * (N / 24);
  int SCALE_UNIT = N / 24;
  int side_points = 6 * SCALE_UNIT;
  int bc_index = 0;
  if (bc == BoundaryCondition::STOKES) {
    boundary_values = ie_Mat(2 * N, 1);
  } else {
    boundary_values = ie_Mat(N, 1);
  }
  // Currently, we set boundary conditions in initialize, because we scale the
  // points after drawing. TODO(John) scale before, so bc init can go in draw fn
  draw_squiggle(bc_index, side_points, 0, 0, 0, 3 * M_PI);
  bc_index += side_points;
  draw_squiggle(bc_index, side_points, 0, 3 * M_PI, 3 * M_PI, 3 * M_PI);
  bc_index += side_points;
  draw_squiggle(bc_index, side_points, 3 * M_PI, 3 * M_PI, 3 * M_PI, 0);
  bc_index += side_points;
  draw_squiggle(bc_index, side_points, 3 * M_PI, 0, 0, 0);
  bc_index += side_points;

  // Brave post-processing - lets put everyone in [0,1]x[0,1]

  // // points
  for (unsigned int i = 0; i < points.size(); i += 2) {
    points[i] = (1.05 + points[i]) * (0.9 / (3 * M_PI + 2.05));
    points[i + 1] = (1.05 + points[i + 1]) * (0.9 / (3 * M_PI + 2.05));
    double potential;
    switch (boundary_condition) {
      case BoundaryCondition::SINGLE_ELECTRON:
        potential = log(sqrt(pow(points[i] + 2, 2) +
                             pow(points[i + 1] + 2, 2))) / (2 * M_PI);
        boundary_values.set(i / 2, 0, potential);
        break;
      case BoundaryCondition::ALL_ONES:
        boundary_values.set(i / 2, 0, 1.0);
        break;
      case BoundaryCondition::BUMP_FUNCTION: {
        double N = boundary_values.height();
        // This x_val is -1 at i=0 and +1 at i=N-1
        double x_val = ((2.0 * i + 1.0 - N) / (N - 1.0));
        potential = exp(-1.0 / (1.0 - pow(x_val, 2)));
        boundary_values.set(i / 2, 0, potential);
        break;
      }
      case BoundaryCondition::STOKES:
        boundary_values.set(i, 0, -normals[i + 1]);
        boundary_values.set(i + 1, 0, normals[i]);
        break;
    }
  }

  // normals stay the same
  // weights
  for (unsigned int i = 0; i < weights.size(); i++) {
    weights[i] = weights[i] * (0.9 / (3 * M_PI + 2.05));
  }

  // curvatures
  for (unsigned int i = 0; i < curvatures.size(); i++) {
    curvatures[i] = curvatures[i] / (0.9 / (3 * M_PI + 2.05));
  }

}


bool Squiggly::is_in_domain(const Vec2 & a) {
  // Hacky, TODO(John) don't be hacky
  // Unscale the point, compare to original representation
  double x = -1.05 + (a.a[0] / (0.9 / (3 * M_PI + 2.05)));
  double y = -1.05 + (a.a[1] / (0.9 / (3 * M_PI + 2.05)));
  // double x = a.a[0];
  // double y = a.a[1];

  double eps = 1e-1;

  if (y - eps > -sin(x) && y + eps <  sin(x) + 3 * M_PI &&
      x - eps > -sin(y) && x + eps <  sin(y) + 3 * M_PI) {
    return true;
  }
  return false;
}


}  // namespace ie_solver
