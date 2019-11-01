// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie_solver/kernel/kernel.h"

namespace ie_solver {

ie_Mat Kernel::get(const Dof& tgt, const Dof& src) const {
  ie_Mat tensor;
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      tensor = laplace_kernel(tgt, src);
      break;
    case ie_solver_config::Pde::LAPLACE_NEUMANN:
      tensor = laplace_neumann_kernel(tgt, src);
      break;
    case ie_solver_config::Pde::STOKES:
      tensor = stokes_kernel(tgt, src);
      break;
  }
  return tensor;
}


ie_Mat Kernel::stokes_kernel(const Dof & tgt, const Dof & src) const {
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


ie_Mat Kernel::laplace_kernel(const Dof & tgt, const Dof & src) const {
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


ie_Mat Kernel::laplace_neumann_kernel(const Dof & tgt, const Dof & src) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat tensor(1, 1);
  if (tgt.point.a[0] == src.point.a[0] && tgt.point.a[1] == src.point.a[1]) {
    tensor.set(0, 0, -0.5 + 0.5 * src.curvature * src.weight * scale
               + src.weight);
    return tensor;
  }
  Vec2 r = tgt.point - src.point;
  tensor.set(0, 0, src.weight * scale * (r.dot(tgt.normal)) / (r.dot(
               r)) + src.weight);
  return tensor;
}


ie_Mat Kernel::laplace_neumann_kernel_forward(const Dof & tgt,
    const Dof & src) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat tensor(1, 1);

  Vec2 r = tgt.point - src.point;
  tensor.set(0, 0, src.weight * scale * log(r.norm()));
  return tensor;
}


// This function stores the DoF data,  and calculates the diagonals of the mat
void Kernel::load(Boundary * boundary_,
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
ie_Mat Kernel::operator()(const std::vector<unsigned int>& I_,
                          const std::vector<unsigned int>& J_,
                          double * timing) const {
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      return fast_laplace_get(I_, J_, timing);
      break;
    case ie_solver_config::Pde::LAPLACE_NEUMANN:
      return fast_laplace_neumann_get(I_, J_, timing);
      break;
    case ie_solver_config::Pde::STOKES:
      return fast_stokes_get(I_, J_, timing);
      break;
  }
}


ie_Mat Kernel::fast_laplace_get(const std::vector<unsigned int>& I_,
                                const std::vector<unsigned int>& J_,
                                double * timing) const {
  double start, end;
  if (timing != nullptr) {
    start = omp_get_wtime();
  }
  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (unsigned int j = 0; j < J_.size(); j++) {
    unsigned int src_ind = J_[j];

    double sp1 = boundary->points[2 * src_ind];
    double sp2 =  boundary->points[2 * src_ind + 1];
    double sn1 =  boundary->normals[2 * src_ind];
    double sn2 = boundary->normals[2 * src_ind + 1];
    double sw =  boundary->weights[src_ind];
    double sc = boundary->curvatures[src_ind];

    for (unsigned int i = 0; i < I_.size(); i++) {
      unsigned int tgt_ind = I_[i];

      double tp1 = boundary->points[2 * tgt_ind];
      double tp2 = boundary->points[2 * tgt_ind + 1];

      if (tp1 == sp1 && sp2 == tp2) {
        ret.mat[i + olda_ * j] =  0.5 + 0.5 * sc * sw * scale;
      } else {
        double r0 = tp1 - sp1;
        double r1 = tp2 - sp2;
        ret.mat[i + olda_ * j] = -sw * scale * (r0 * sn1 + r1 * sn2) /
                                 (r0 * r0 + r1 * r1);
      }
    }
  }
  if (timing != nullptr) {
    end = omp_get_wtime();
    *timing += end - start;
  }
  return ret;
}


ie_Mat Kernel::fast_laplace_neumann_get(const std::vector<unsigned int>& I_,
                                        const std::vector<unsigned int>& J_,
                                        double * timing) const {
  double start, end;
  if (timing != nullptr) {
    start = omp_get_wtime();
  }
  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (unsigned int j = 0; j < J_.size(); j++) {
    //   for (unsigned int i = 0; i < I_.size(); i++) {
    //     ret.mat[i + olda_ * j] = get(I_[i], J_[j]);
    //   }
    // }
    unsigned int src_ind = J_[j];

    double sp1 = boundary->points[2 * src_ind];
    double sp2 =  boundary->points[2 * src_ind + 1];
    double sw =  boundary->weights[src_ind];
    double sc = boundary->curvatures[src_ind];

    for (unsigned int i = 0; i < I_.size(); i++) {
      unsigned int tgt_ind = I_[i];

      double tp1 = boundary->points[2 * tgt_ind];
      double tp2 = boundary->points[2 * tgt_ind + 1];
      double tn1 = boundary->normals[2 * tgt_ind];
      double tn2 = boundary->normals[2 * tgt_ind + 1];
      if (tp1 == sp1 && sp2 == tp2) {
        ret.mat[i + olda_ * j] =  sw - 0.5 + 0.5 * sc * sw * scale;
      } else {
        double r0 = tp1 - sp1;
        double r1 = tp2 - sp2;
        ret.mat[i + olda_ * j] = sw + sw * scale * (r0 * tn1 + r1 * tn2) /
                                 (r0 * r0 + r1 * r1);
      }
    }
  }
  if (timing != nullptr) {
    end = omp_get_wtime();
    *timing += end - start;
  }
  return ret;
}


ie_Mat Kernel::fast_stokes_get(const std::vector<unsigned int>& I_,
                               const std::vector<unsigned int>& J_,
                               double * timing) const {
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


void Kernel::make_id_mat(ie_Mat* mat, const QuadTree* tree,
                         const QuadTreeNode* node, bool strong_admissibility) const {
  double cntr_x = node->corners[0] + node->side_length / 2.0;
  double cntr_y = node->corners[1] + node->side_length / 2.0;

  std::vector<unsigned int> active_box = node->src_dof_lists.active_box;
  if (!strong_admissibility) {
    // Grab all points inside the proxy circle which are outside the box
    std::vector<unsigned int> inner_circle, outside_box;
    // int num_outside_circle = 0;

    // for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
    //   if (level_node->id != node->id) {
    for (unsigned int matrix_index : node->src_dof_lists.near) {
      // outside_box.push_back(matrix_index);
      unsigned int point_index = matrix_index / solution_dimension;
      unsigned int points_vec_index = point_index * domain_dimension;
      double x = tree->boundary->points[points_vec_index];
      double y = tree->boundary->points[points_vec_index + 1];
      double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
      if (dist < RADIUS_RATIO * node->side_length) {
        inner_circle.push_back(matrix_index);
      }
      // else {
      //   num_outside_circle++;
      // }
    }
    //   }
    // }

    // TODO(John) reimplement this
    // We only do a proxy circle if there are lots outside the circle
    // if (num_outside_circle < NUM_PROXY_POINTS) {
    //   double pxy_ind_str = omp_get_wtime();

    //   // Now all the matrices are gathered, put them into *mat.
    //   *mat = ie_Mat(2 * outside_box.size(), active_box.size());
    //   mat->set_submatrix(0, outside_box.size(), 0, active_box.size(),
    //                      (*this)(outside_box, active_box), false, true);
    //   mat->set_submatrix(outside_box.size(), 2 * outside_box.size(),
    //                      0, active_box.size(), (*this)(active_box, outside_box),
    //                      true, true);
    // } else {

    // Now all the matrices are gathered, put them into *mat.
    *mat = ie_Mat(2 * inner_circle.size() + solution_dimension * 2 *
                  NUM_PROXY_POINTS, active_box.size());

    double knl_str = omp_get_wtime();
    mat->set_submatrix(0, inner_circle.size(),
                       0, active_box.size(), (*this)(inner_circle, active_box),
                       false, true);
    mat->set_submatrix(inner_circle.size(), 2 * inner_circle.size(),
                       0, active_box.size(), (*this)(active_box, inner_circle),
                       true, true);
    double knl_end = omp_get_wtime();
    ie_Mat::kernel_time += (knl_end - knl_str);

    double pxy_start = omp_get_wtime();
    mat->set_submatrix(2 * inner_circle.size(),  solution_dimension * 2 *
                       NUM_PROXY_POINTS + 2 * inner_circle.size(), 0,
                       active_box.size(),  make_proxy_mat(cntr_x, cntr_y, node->side_length
                           * RADIUS_RATIO, tree, active_box), false, true);
    double pxy_end = omp_get_wtime();
    ie_Mat::proxy_time += (pxy_end - pxy_start);
    // }
  } else {
    *mat = make_proxy_mat(cntr_x, cntr_y, node->side_length * 1.5, tree,
                          active_box);
  }
}


ie_Mat Kernel::make_proxy_mat(double cntr_x, double cntr_y,
                              double r, const QuadTree * tree,
                              const std::vector<unsigned int>& box_inds) const {
  // each row is a pxy point, cols are box dofs
  double proxy_weight = 2.0 * M_PI * r / NUM_PROXY_POINTS;
  double proxy_curvature = 1.0 / r;
  std::vector<double>  pxy_p, pxy_n;
  for (int i = 0; i < NUM_PROXY_POINTS; i++) {
    double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
    pxy_p.push_back(cntr_x + r * cos(ang));
    pxy_p.push_back(cntr_y + r * sin(ang));
    pxy_n.push_back(cos(ang));
    pxy_n.push_back(sin(ang));
  }

  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      return fast_laplace_proxy_get(pxy_p, pxy_n, proxy_curvature ,
                                    proxy_weight, box_inds); TODO(JOHN) test this
    case ie_solver_config::Pde::LAPLACE_NEUMANN:
      TODO(JOHN)
    case ie_solver_config::Pde::STOKES:
      return fast_stokes_proxy_get(pxy_p, pxy_n, proxy_curvature , proxy_weight,
                                   box_inds);
  }
}


ie_Mat Kernel::fast_stokes_proxy_get(const std::vector<double> & pxy_p,
                                     const std::vector<double> & pxy_n,
                                     double pxy_c, double pxy_w,
                                     const std::vector<unsigned int> & box_inds) const {
  double scale = 1.0 / (M_PI);
  ie_Mat ret(2 * pxy_p.size(), box_inds.size());
  int lda = 2 * pxy_p.size();

  // Then Active to Proxy
  for (unsigned int j = 0; j < box_inds.size(); j++) {
    unsigned int src_ind = box_inds[j];
    unsigned int j_point_index = src_ind / 2;
    unsigned int j_points_vec_index = j_point_index * 2;

    double sp1 = boundary->points[j_points_vec_index];
    double sp2 =  boundary->points[j_points_vec_index + 1];
    double sn1 =  boundary->normals[j_points_vec_index];
    double sn2 = boundary->normals[j_points_vec_index + 1];
    double sw =  boundary->weights[j_point_index];
    double sc = boundary->curvatures[j_point_index];

    for (unsigned int i = 0; i < pxy_p.size(); i += 2) {
      double tp1 = pxy_p[i];
      double tp2 = pxy_p[i + 1];
      double tn1 = pxy_n[i];
      double tn2 = pxy_n[i + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;
      double potential = sw * scale * (r0 * sn1 + r1 * sn2) /
                         (pow(r0 * r0 + r1 * r1, 2));
      if (src_ind % 2 == 0) {
        ret.mat[i + lda * j] =
          potential * r0 * r0 + sw * tn1 * sn1;
        ret.mat[ 1 + i + lda * j] =
          potential * r1 * r0 + sw * tn2 * sn1;
      } else {
        ret.mat[ i + lda * j] =
          potential * r0 * r1 + sw * tn1 * sn2;
        ret.mat[1 + i + lda * j] =
          potential * r1 * r1 + sw * tn2 * sn2;
      }
    }
  }

  // First Proxy To Active
  for (unsigned int j = 0; j < pxy_p.size(); j += 2) {
    double sp1 = pxy_p[j];
    double sp2 =  pxy_p[j + 1];
    double sn1 =  pxy_n[j];
    double sn2 = pxy_n[j + 1];
    double sw =  pxy_w;
    double sc = pxy_c;
// TODO(John) attempt calculating normals/points on the fly?
    for (unsigned int i = 0; i < box_inds.size(); i++) {
      unsigned int tgt_ind = box_inds[i];
      unsigned int i_point_index = tgt_ind / 2;
      unsigned int i_points_vec_index = i_point_index * 2;

      double tp1 = boundary->points[i_points_vec_index];
      double tp2 = boundary->points[i_points_vec_index + 1];
      double tn1 = boundary->normals[i_points_vec_index];
      double tn2 = boundary->normals[i_points_vec_index + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;
      double potential = sw * scale * (r0 * sn1 + r1 * sn2) /
                         (pow(r0 * r0 + r1 * r1, 2));
      if (tgt_ind % 2 == 0) {
        ret.mat[pxy_p.size() + j + i * lda] =
          potential * r0 * r0 + sw * tn1 * sn1;
        ret.mat[pxy_p.size() + j + 1 + i * lda] =
          potential * r0 * r1 + sw * tn1 * sn2;
      } else {
        ret.mat[pxy_p.size() + j + i * lda] =
          potential * r1 * r0 + sw * tn2 * sn1;
        ret.mat[pxy_p.size() + j + 1 + i * lda] =
          potential * r1 * r1 + sw * tn2 * sn2;
      }
    }
  }
  return ret;
}



ie_Mat Kernel::fast_laplace_proxy_get(const std::vector<double> & pxy_p,
                                      const std::vector<double> & pxy_n,
                                      double pxy_c, double pxy_w,
                                      const std::vector<unsigned int> & box_inds) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(2 * pxy_p.size(), box_inds.size());
  int lda = 2 * pxy_p.size();

  // Then Active to Proxy
  for (unsigned int j = 0; j < box_inds.size(); j++) {
    unsigned int src_ind = box_inds[j];
    unsigned int j_point_index = src_ind / 2;
    unsigned int j_points_vec_index = j_point_index * 2;

    double sp1 = boundary->points[j_points_vec_index];
    double sp2 =  boundary->points[j_points_vec_index + 1];
    double sn1 =  boundary->normals[j_points_vec_index];
    double sn2 = boundary->normals[j_points_vec_index + 1];
    double sw =  boundary->weights[j_point_index];
    double sc = boundary->curvatures[j_point_index];

    for (unsigned int i = 0; i < pxy_p.size(); i += 2) {
      double tp1 = pxy_p[i];
      double tp2 = pxy_p[i + 1];
      double tn1 = pxy_n[i];
      double tn2 = pxy_n[i + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;
      ret.mat[i + lda * j] = -sw * scale * (r0 * sn1 + r1 * sn2) /
                             (r0 * r0 + r1 * r1);
    }
  }

  // First Proxy To Active
  for (unsigned int j = 0; j < pxy_p.size(); j += 2) {
    double sp1 = pxy_p[j];
    double sp2 =  pxy_p[j + 1];
    double sn1 =  pxy_n[j];
    double sn2 = pxy_n[j + 1];
    double sw =  pxy_w;
    double sc = pxy_c;
// TODO(John) attempt calculating normals/points on the fly?
    for (unsigned int i = 0; i < box_inds.size(); i++) {
      unsigned int tgt_ind = box_inds[i];
      unsigned int i_point_index = tgt_ind / 2;
      unsigned int i_points_vec_index = i_point_index * 2;

      double tp1 = boundary->points[i_points_vec_index];
      double tp2 = boundary->points[i_points_vec_index + 1];
      double tn1 = boundary->normals[i_points_vec_index];
      double tn2 = boundary->normals[i_points_vec_index + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;

      ret.mat[pxy_p.size() + i + lda * j] = -sw * scale *
                                            (r0 * sn1 + r1 * sn2) /
                                            (r0 * r0 + r1 * r1);

    }
  }
  return ret;
}



















// ie_Mat Kernel::fast_laplace_get(const std::vector<Dof>& tgts,
//                                 const std::vector<Dof>& srcs,
//                                 double* timing) const {
//   double start, end;
//   if (timing != nullptr) {
//     start = omp_get_wtime();
//   }
//   double scale = 1.0 / (2 * M_PI);
//   ie_Mat ret(tgts.size(), srcs.size());
//   int olda_ = tgts.size();
//   for (unsigned int j = 0; j < srcs.size(); j++) {
//     Dof src = srcs[j];

//     double sp1 = src.point.a[0];
//     double sp2 =  src.point.a[1];
//     double sn1 =  src.normal.a[0];
//     double sn2 = src.normal.a[1];
//     double sw =  src.weight;
//     double sc = src.curvature;

//     for (unsigned int i = 0; i < tgts.size(); i++) {
//       Dof tgt = tgts[i];

//       double tp1 = tgt.point.a[0];
//       double tp2 =  tgt.point.a[1];
//       double tn1 =  tgt.normal.a[0];
//       double tn2 = tgt.normal.a[1];

//       if (tp1 == sp1 && sp2 == tp2) {
//         ret.mat[i + olda_ * j] =  0.5 + 0.5 * sc * sw * scale;
//       } else {
//         double r0 = tp1 - sp1;
//         double r1 = tp2 - sp2;
//         ret.mat[i + olda_ * j] = -sw * scale * (r0 * sn1 + r1 * sn2) /
//                                  (r0 * r0 + r1 * r1);
//       }
//     }
//   }
//   if (timing != nullptr) {
//     end = omp_get_wtime();
//     *timing += end - start;
//   }
//   return ret;
// }


// ie_Mat Kernel::fast_stokes_get(const std::vector<Dof>& tgts,
//                                const std::vector<Dof>& srcs,
//                                double* timing) const {
//   double start, end;
//   double scale = 1.0 / (M_PI);

//   if (timing != nullptr) {
//     start = omp_get_wtime();
//   }
//   ie_Mat ret(2 * tgts.size(), 2 * srcs.size());
//   int olda_ = 2 * tgts.size();

//   for (unsigned int j = 0; j < srcs.size(); j++) {
//     Dof src = srcs[j];

//     double sp1 = src.point.a[0];
//     double sp2 =  src.point.a[1];
//     double sn1 =  src.normal.a[0];
//     double sn2 = src.normal.a[1];
//     double sw =  src.weight;
//     double sc = src.curvature;

//     for (unsigned int i = 0; i < tgts.size(); i++) {
//       Dof tgt = tgts[i];

//       double tp1 = tgt.point.a[0];
//       double tp2 =  tgt.point.a[1];
//       double tn1 =  tgt.normal.a[0];
//       double tn2 = tgt.normal.a[1];

//       if (tp1 == sp1 && sp2 == tp2) {
//         double potential = - 0.5 * sc * sw * scale;


//         ret.mat[2 * i + olda_ * (2 * j)] = -0.5
//                                            + potential * sn2 * sn2 +  sw *
//                                            tn1 * sn1;

//         ret.mat[2 * i + olda_ * (2 * j + 1)] = -potential * sn1 * sn2
//                                                 + sw * tn2 * sn1;

//         ret.mat[2 * i + 1 + olda_ * (2 * j)] =  -potential * sn1 * sn2
//                                                 + sw * tn2 * sn1;

//         ret.mat[2 * i + 1 + olda_ * (2 * j + 1)] = -0.5 + potential * sn1
//                                                    * sn1 + sw * tn2 * sn2;

//       } else {
//         double r0 = tp1 - sp1;
//         double r1 = tp2 - sp2;
//         double potential = sw * scale * (r0 * sn1 + r1 * sn2) /
//                            (pow(r0 * r0 + r1 * r1, 2));
//         ret.mat[2 * i + olda_ * (2 * j)] = potential * r0 * r0 + sw * tn1 *
//                                            sn1;
//         ret.mat[2 * i + olda_ * (2 * j + 1)] = potential * r0 * r1
//                                                + sw * tn1 * sn2;
//         ret.mat[2 * i + 1 + olda_ * (2 * j)] = potential * r1 * r0
//                                                + sw * tn2 * sn1;
//         ret.mat[2 * i + 1 + olda_ * (2 * j + 1)] = potential * r1 * r1
//                                                    + sw * tn2 * sn2;
//       }
//     }
//   }
//   if (timing != nullptr) {
//     end = omp_get_wtime();
//     *timing += end - start;
//   }
//   return ret;
// }


// double Kernel::forward_get(unsigned int tgt_ind,
//                            unsigned int src_ind) const {
//   Dof tgt, src;
//   unsigned int i_point_index = tgt_ind / solution_dimension;
//   unsigned int i_points_vec_index = i_point_index * domain_dimension;
//   unsigned int j_point_index = src_ind / solution_dimension;
//   unsigned int j_points_vec_index = j_point_index * domain_dimension;

//   tgt.point = Vec2(domain_points[i_points_vec_index],
//                    domain_points[i_points_vec_index + 1]);
//   tgt.is_boundary = false;

//   src.point = Vec2(boundary->points[j_points_vec_index],
//                    boundary->points[j_points_vec_index + 1]);
//   src.normal = Vec2(boundary->normals[j_points_vec_index],
//                     boundary->normals[j_points_vec_index + 1]);
//   src.curvature = boundary->curvatures[j_point_index];
//   src.weight = boundary->weights[j_point_index];
//   src.is_boundary = true;
//   ie_Mat tensor = get(tgt, src);
//   return tensor.get(tgt_ind % solution_dimension,
//                     src_ind % solution_dimension);
// }


// double Kernel::get(unsigned int tgt_ind, unsigned int src_ind) const {
//   Dof tgt, src;
//   unsigned int i_point_index = tgt_ind / solution_dimension;
//   unsigned int i_points_vec_index = i_point_index * domain_dimension;
//   unsigned int j_point_index = src_ind / solution_dimension;
//   unsigned int j_points_vec_index = j_point_index * domain_dimension;

//   tgt.point = Vec2(boundary->points[i_points_vec_index],
//                    boundary->points[i_points_vec_index + 1]);
//   tgt.normal = Vec2(boundary->normals[i_points_vec_index],
//                     boundary->normals[i_points_vec_index + 1]);
//   tgt.curvature = boundary->curvatures[i_point_index];
//   tgt.weight = boundary->weights[i_point_index];
//   tgt.is_boundary = true;

//   src.point = Vec2(boundary->points[j_points_vec_index],
//                    boundary->points[j_points_vec_index + 1]);
//   src.normal = Vec2(boundary->normals[j_points_vec_index],
//                     boundary->normals[j_points_vec_index + 1]);
//   src.curvature = boundary->curvatures[j_point_index];
//   src.weight = boundary->weights[j_point_index];
//   src.is_boundary = true;

//   ie_Mat tensor = get(tgt, src);
//   return tensor.get(tgt_ind % solution_dimension,
//                     src_ind % solution_dimension);
// }

// ie_Mat Kernel::forward_get(const std::vector<unsigned int>& I_,
//                            const std::vector<unsigned int>& J_) const {
//   ie_Mat ret(I_.size(), J_.size());
//   int olda_ = I_.size();
//   for (unsigned int i = 0; i < I_.size(); i++) {
//     for (unsigned int j = 0; j < J_.size(); j++) {
//       ret.mat[i + olda_ * j] = forward_get(I_[i], J_[j]);
//     }
//   }
//   return ret;
// }


// ie_Mat Kernel::operator()(const std::vector<Dof>& tgts,
//                           const std::vector<Dof>& srcs) const {
//   ie_Mat ret(solution_dimension * tgts.size(),
//              solution_dimension * srcs.size());
//   for (unsigned int i = 0; i < tgts.size(); i++) {
//     for (unsigned int j = 0; j < srcs.size(); j++) {
//       ie_Mat tensor = get(tgts[i], srcs[j]);
//       for (int k = 0; k < solution_dimension; k++) {
//         for (int l = 0; l < solution_dimension; l++) {
//           ret.set(solution_dimension * i + k,
//                   solution_dimension * j + l, tensor.get(k, l));
//         }
//       }
//     }
//   }
//   return ret;
// }


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
//         std::cout << "Error: check Laplace called on Stokes BC." <<std::endl;
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
//     std::cout << "Error: cannot currently check stokes error on non-annulus"
//               << std::endl;
//     return -1;
//   }
//   if (boundary->holes.size() != 1) {
//     std::cout << "Error: can only check error on boundary with one hole" <<
//               std::endl;
//     return -1;
//   } else if (boundary->holes[0].center.a[0] != 0.5
//              || boundary->holes[0].center.a[1] != 0.5
//              || boundary->holes[0].radius != 0.05) {
//     std::cout << "Error: can only check error on boundary with hole at "
//               << "center and radius 0.05" << std::endl;
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
//     Vec2 sol = Vec2(domain_solution.get(i, 0),
//                     domain_solution.get(i + 1, 0));

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
