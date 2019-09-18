// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include <omp.h>
#include "ie_solver/kernel/kernel.h"

namespace ie_solver {

ie_Mat Kernel::get(const Dof& tgt, const Dof& src) const {
  ie_Mat tensor;
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      tensor = laplace_kernel(tgt, src);
      break;
    case ie_solver_config::Pde::STOKES:
      tensor = stokes_kernel(tgt, src);
      break;
  }
  return tensor;
}


double Kernel::forward_get(unsigned int tgt_ind, unsigned int src_ind) const {
  Dof tgt, src;
  unsigned int i_point_index = tgt_ind / solution_dimension;
  unsigned int i_points_vec_index = i_point_index * domain_dimension;
  unsigned int j_point_index = src_ind / solution_dimension;
  unsigned int j_points_vec_index = j_point_index * domain_dimension;

  tgt.point = Vec2(domain_points[i_points_vec_index],
                   domain_points[i_points_vec_index + 1]);
  tgt.is_boundary = false;

  src.point = Vec2(boundary->points[j_points_vec_index],
                   boundary->points[j_points_vec_index + 1]);
  src.normal = Vec2(boundary->normals[j_points_vec_index],
                    boundary->normals[j_points_vec_index + 1]);
  src.curvature = boundary->curvatures[j_point_index];
  src.weight = boundary->weights[j_point_index];
  src.is_boundary = true;
  ie_Mat tensor = get(tgt, src);
  return tensor.get(tgt_ind % solution_dimension, src_ind % solution_dimension);
}


double Kernel::get(unsigned int tgt_ind, unsigned int src_ind) const {
  Dof tgt, src;
  unsigned int i_point_index = tgt_ind / solution_dimension;
  unsigned int i_points_vec_index = i_point_index * domain_dimension;
  unsigned int j_point_index = src_ind / solution_dimension;
  unsigned int j_points_vec_index = j_point_index * domain_dimension;

  tgt.point = Vec2(boundary->points[i_points_vec_index],
                   boundary->points[i_points_vec_index + 1]);
  tgt.normal = Vec2(boundary->normals[i_points_vec_index],
                    boundary->normals[i_points_vec_index + 1]);
  tgt.curvature = boundary->curvatures[i_point_index];
  tgt.weight = boundary->weights[i_point_index];
  tgt.is_boundary = true;

  src.point = Vec2(boundary->points[j_points_vec_index],
                   boundary->points[j_points_vec_index + 1]);
  src.normal = Vec2(boundary->normals[j_points_vec_index],
                    boundary->normals[j_points_vec_index + 1]);
  src.curvature = boundary->curvatures[j_point_index];
  src.weight = boundary->weights[j_point_index];
  src.is_boundary = true;

  ie_Mat tensor = get(tgt, src);
  return tensor.get(tgt_ind % solution_dimension, src_ind % solution_dimension);
}


ie_Mat Kernel::stokes_kernel(const Dof& tgt, const Dof& src) const {
  // double layer
  double scale = 1.0 / (M_PI);

  ie_Mat tensor(2, 2);
  ie_Mat normop(2, 2);
  if (tgt.is_boundary) {
    normop.set(0, 0, src.weight * tgt.normal.a[0] * src.normal.a[0]);
    normop.set(0, 1, src.weight * tgt.normal.a[0] * src.normal.a[1]);
    normop.set(1, 0, src.weight * tgt.normal.a[1] * src.normal.a[0]);
    normop.set(1, 1, src.weight * tgt.normal.a[1] * src.normal.a[1]);
  }
  if (tgt.point.a[0] == src.point.a[0] && tgt.point.a[1] == src.point.a[1]) {
    // tangent
    double t0 = -src.normal.a[1];
    double t1 =  src.normal.a[0];
    double potential = - 0.5 * src.curvature * src.weight * scale;
    tensor.set(0, 0, -0.5 + potential * t0 * t0);
    tensor.set(1, 1, -0.5 + potential * t1 * t1);
    tensor.set(0, 1, potential * t0 * t1);
    tensor.set(1, 0, potential * t1 * t0);

    return tensor + normop;
  }
  Vec2 r = tgt.point - src.point;
  double r0 = r.a[0];
  double r1 = r.a[1];
  double potential = src.weight * scale * (r.dot(src.normal)) /
                     (pow(r.dot(r), 2));
  tensor.set(0, 0, potential * r0 * r0);
  tensor.set(1, 1, potential * r1 * r1);
  tensor.set(0, 1, potential * r0 * r1);
  tensor.set(1, 0, potential * r1 * r0);
  return tensor + normop;
}


ie_Mat Kernel::laplace_kernel(const Dof& tgt, const Dof& src) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat tensor(1, 1);
  if (tgt.point.a[0] == src.point.a[0] && tgt.point.a[1] == src.point.a[1]) {
    tensor.set(0, 0, 0.5 + 0.5 * src.curvature * src.weight * scale);
    return tensor;
  }
  Vec2 r = tgt.point - src.point;
  tensor.set(0, 0, -src.weight * scale * (r.dot(src.normal)) / (r.dot(r)));
  return tensor;
}

// This function stores the DoF data,  and calculates the diagonals of the mat
void Kernel::load(Boundary* boundary_,
                  const std::vector<double>& domain_points_,
                  ie_solver_config::Pde pde_, int solution_dimension_,
                  int domain_dimension_) {
  solution_dimension = solution_dimension_;
  domain_points = domain_points_;
  domain_dimension = domain_dimension_;
  boundary = boundary_;
  pde = pde_;
}


// // TODO(John) shouldn't this->I have the underscore after it, not this arg?
// ie_Mat Kernel::operator()(const std::vector<unsigned int>& I_,
//                           const std::vector<unsigned int>& J_, double* timing) const {
//   double start, end;
//   if (timing != nullptr) {
//     start = omp_get_wtime();
//   }
//   ie_Mat ret(I_.size(), J_.size());
//   int olda_ = I_.size();
//   for (unsigned int j = 0; j < J_.size(); j++) {
//     for (unsigned int i = 0; i < I_.size(); i++) {
//       ret.mat[i + olda_ * j] = get(I_[i], J_[j]);
//     }
//   }
//   if (timing != nullptr) {
//     end = omp_get_wtime();
//     *timing += end - start;
//   }
//   return ret;
// }


// TODO(John) shouldn't this->I have the underscore after it, not this arg?
ie_Mat Kernel::operator()(const std::vector<unsigned int>& I_,
                          const std::vector<unsigned int>& J_, double* timing) const {
  double start, end;
  double scale = 1.0 / (M_PI);

  if (timing != nullptr) {
    start = omp_get_wtime();
  }
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (unsigned int j = 0; j < J_.size(); j++) {
    unsigned int src_ind = J_[j];
    unsigned int j_point_index = src_ind / 2;
    unsigned int j_points_vec_index = j_point_index * 2;

    double sp1 = boundary->points[j_points_vec_index];
    double sp2 =  boundary->points[j_points_vec_index + 1];
    double sn1 =  boundary->normals[j_points_vec_index];
    double sn2 = boundary->normals[j_points_vec_index + 1];
    double sw =  boundary->weights[j_point_index];
    double sc = boundary->curvatures[j_point_index];

    for (unsigned int i = 0; i < I_.size(); i++) {
      unsigned int tgt_ind = I_[i];
      unsigned int i_point_index = tgt_ind / 2;
      unsigned int i_points_vec_index = i_point_index * 2;

      double tp1 = boundary->points[i_points_vec_index];
      double tp2 = boundary->points[i_points_vec_index + 1];
      double tn1 = boundary->normals[i_points_vec_index];
      double tn2 = boundary->normals[i_points_vec_index + 1];


      if (tp1 == sp1 && sp2 == tp2) {
        double potential = - 0.5 * sc * sw * scale;

        if (tgt_ind % 2 == 0) {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] = -0.5 + potential * sn2 * sn2 +  sw *
                                     tn1 * sn1;
          } else {
            ret.mat[i + olda_ * j] = -potential * sn1 * sn2 + sw * tn2 *
                                     sn1;
          }
        } else {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] =  -potential * sn1 * sn2 + sw * tn2 *
                                      sn1;
          } else {
            ret.mat[i + olda_ * j] = -0.5 + potential * sn1 * sn1 + sw *
                                     tn2 * sn2;
          }
        }
      } else {
        double r0 = tp1 - sp1;
        double r1 = tp2 - sp2;
        double potential = sw * scale * (r0 * sn1 + r1 * sn2) /
                           (pow(r0 * r0 + r1 * r1, 2));
        if (tgt_ind % 2 == 0) {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] = potential * r0 * r0 + sw * tn1 *
                                     sn1;
          } else {
            ret.mat[i + olda_ * j] = potential * r0 * r1 + sw * tn1 *
                                     sn2;
          }
        } else {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] = potential * r1 * r0 + sw * tn2 *
                                     sn1;
          } else {
            ret.mat[i + olda_ * j] = potential * r1 * r1 + sw * tn2 *
                                     sn2;
          }
        }
      }

    }
  }
  if (timing != nullptr) {
    end = omp_get_wtime();
    *timing += end - start;
  }
  return ret;
}



ie_Mat Kernel::forward_get(const std::vector<unsigned int>& I_,
                           const std::vector<unsigned int>& J_) const {
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (unsigned int i = 0; i < I_.size(); i++) {
    for (unsigned int j = 0; j < J_.size(); j++) {
      ret.mat[i + olda_ * j] = forward_get(I_[i], J_[j]);
    }
  }
  return ret;
}


ie_Mat Kernel::operator()(const std::vector<Dof>& tgts,
                          const std::vector<Dof>& srcs) const {
  ie_Mat ret(solution_dimension * tgts.size(),
             solution_dimension * srcs.size());
  for (unsigned int i = 0; i < tgts.size(); i++) {
    for (unsigned int j = 0; j < srcs.size(); j++) {
      ie_Mat tensor = get(tgts[i], srcs[j]);
      for (int k = 0; k < solution_dimension; k++) {
        for (int l = 0; l < solution_dimension; l++) {
          ret.set(solution_dimension * i + k,
                  solution_dimension * j + l, tensor.get(k, l));
        }
      }
    }
  }
  return ret;
}


// double laplace_error(const ie_Mat & domain, double id_tol,
//                      const std::vector<double>& domain_points,
//                      Boundary * boundary) {
//   if (boundary->holes.size() > 0) {
//     std::cout <<
//               "Error: laplace error not currently calculated for multiply"
//               << " connected domain." << std::endl;
//     return -1;
//   }
//   double max = 0;
//   double diff_norm = 0;
//   double avg = 0;
//   double norm_of_true = 0;
//   for (unsigned int i = 0; i < domain_points.size(); i += 2) {
//     double x0 = domain_points[i];
//     double x1 = domain_points[i + 1];
//     Vec2 x(x0, x1);
//     if (!boundary->is_in_domain(x)) {
//       continue;
//     }
//     double potential;
//     switch (boundary->boundary_condition) {
//       case BoundaryCondition::SINGLE_ELECTRON:
//         potential = log(sqrt(pow(x0 + 2, 2) + pow(x1 + 2, 2))) / (2 * M_PI);
//         break;
//       case BoundaryCondition::ALL_ONES:
//         potential = 1.0;
//         break;
//       case BoundaryCondition::BUMP_FUNCTION: {
//         std::cout << "Error: check Laplace called on Bump BC;"
//                   << " no analytic solution known to check against."
//                   << std::endl;
//         break;
//       }
//       case BoundaryCondition::STOKES:
//         std::cout << "Error: check Laplace called on Stokes BC." << std::endl;
//         break;
//     }
//     if (std::isnan(domain.get(i / 2, 0))) {
//       continue;
//     }
//     double diff = std::abs(potential - domain.get(i / 2, 0));
//     avg += diff / potential;
//     diff_norm += pow(diff, 2);
//     norm_of_true += pow(potential, 2);
//     max = std::max(max, diff / potential);
//   }
//   avg /= domain_points.size();
//   diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
//   return diff_norm;
// }


// double stokes_error(const ie_Mat & domain_solution, double id_tol,
//                     const std::vector<double>& domain_points,
//                     Boundary * boundary) {
//   if (boundary->boundary_shape != Boundary::ANNULUS) {
//     std::cout << "Error: cannot currently check stokes error on non-annulus" <<
//               std::endl;
//     return -1;
//   }
//   if (boundary->holes.size() != 1) {
//     std::cout << "Error: can only check error on boundary with one hole" <<
//               std::endl;
//     return -1;
//   } else if (boundary->holes[0].center.a[0] != 0.5
//              || boundary->holes[0].center.a[1] != 0.5
//              || boundary->holes[0].radius != 0.05) {
//     std::cout << "Error: can only check error on boundary with hole at center "
//               << "and radius 0.05" << std::endl;
//     return -1;
//   }
//   double truth_size = 0;
//   double total_diff = 0;
//   for (unsigned int i = 0; i < domain_points.size(); i += 2) {
//     double x0 = domain_points[i];
//     double x1 = domain_points[i + 1];
//     Vec2 x(x0, x1);
//     if (!boundary->is_in_domain(x)) {
//       continue;
//     }
//     Vec2 center(0.5, 0.5);



//     Vec2 r = x - center;
//     Vec2 sol = Vec2(domain_solution.get(i, 0), domain_solution.get(i + 1, 0));

//     Vec2 truth = Vec2(-r.a[1], r.a[0]);
//     truth = truth * (1 / truth.norm());
//     double om1 = -30;
//     double om2 = 4;
//     double r1 = 0.05;
//     double r2 = 0.25;
//     double c1 = (om2 * pow(r2, 2) - om1 * pow(r1, 2))
//                 / (pow(r2, 2) - pow(r1, 2));
//     double c2 = ((om1 - om2) * pow(r2, 2) * pow(r1, 2))
//                 / (pow(r2, 2) - pow(r1, 2));

//     double truth_length = c1 * r.norm() + (c2 / r.norm());

//     truth = truth * truth_length;
//     // domain_solution.set(i, 0, truth.a[0]);
//     // domain_solution.set(i + 1, 0, truth.a[1]);

//     double diff = (truth - sol).norm();
//     truth_size += fabs(truth_length);
//     total_diff += diff;
//   }

//   return total_diff / truth_size;
// }



}  // namespace ie_solver
