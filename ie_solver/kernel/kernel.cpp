// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include "ie_solver/kernel/kernel.h"

namespace ie_solver {

int Kernel::IMPROVE_CONDITION = 0;

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
    tensor.set(0, 0,  - 0.5 + potential * t0 * t0);
    tensor.set(1, 1,  - 0.5 + potential * t1 * t1);
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
    tensor.set(0, 0, 0.5 + 0.5 * src.curvature * src.weight *
               scale);
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
    tensor.set(0, 0,  src.weight - 0.5 + 0.5 * src.curvature * src.weight * scale);
    return tensor;
  }
  Vec2 r = tgt.point - src.point;
  tensor.set(0, 0, src.weight + src.weight * scale * (r.dot(tgt.normal))
             / (r.dot(r)));
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


ie_Mat Kernel::operator()(const std::vector<int>& I_,
                          const std::vector<int>& J_,
                          double * timing) const {
  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      return laplace_get(I_, J_, timing);
      break;
    case ie_solver_config::Pde::LAPLACE_NEUMANN:
      return laplace_neumann_get(I_, J_, timing);
      break;
    default:  // to suppress compiler warning
      // case ie_solver_config::Pde::STOKES:
      return stokes_get(I_, J_, timing);
      break;
  }
}


ie_Mat Kernel::laplace_get(const std::vector<int>& I_,
                           const std::vector<int>& J_,
                           double * timing) const {
  double start, end;
  if (timing != nullptr) {
    start = omp_get_wtime();
  }
  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (int j = 0; j < J_.size(); j++) {
    int src_ind = J_[j];

    double sp1 = boundary->points[2 * src_ind];
    double sp2 =  boundary->points[2 * src_ind + 1];
    double sn1 =  boundary->normals[2 * src_ind];
    double sn2 = boundary->normals[2 * src_ind + 1];
    double sw =  boundary->weights[src_ind];
    double sc = boundary->curvatures[src_ind];

    for (int i = 0; i < I_.size(); i++) {
      int tgt_ind = I_[i];

      double tp1 = boundary->points[2 * tgt_ind];
      double tp2 = boundary->points[2 * tgt_ind + 1];

      if (tp1 == sp1 && sp2 == tp2) {
        ret.mat[i + olda_ * j] = IMPROVE_CONDITION +  0.5
                                 + 0.5 * sc * sw * scale;
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


ie_Mat Kernel::laplace_neumann_get(const std::vector<int>& I_,
                                   const std::vector<int>& J_,
                                   double * timing) const {
  double start, end;
  if (timing != nullptr) {
    start = omp_get_wtime();
  }
  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (int j = 0; j < J_.size(); j++) {
    int src_ind = J_[j];

    double sp1 = boundary->points[2 * src_ind];
    double sp2 =  boundary->points[2 * src_ind + 1];
    double sw =  boundary->weights[src_ind];
    double sc = boundary->curvatures[src_ind];

    for (int i = 0; i < I_.size(); i++) {
      int tgt_ind = I_[i];

      double tp1 = boundary->points[2 * tgt_ind];
      double tp2 = boundary->points[2 * tgt_ind + 1];
      double tn1 = boundary->normals[2 * tgt_ind];
      double tn2 = boundary->normals[2 * tgt_ind + 1];
      if (tp1 == sp1 && sp2 == tp2) {
        ret.mat[i + olda_ * j] =  IMPROVE_CONDITION + sw
                                  - 0.5 + 0.5 * sc * sw * scale;
      } else {
        double r0 = tp1 - sp1;
        double r1 = tp2 - sp2;
        ret.mat[i + olda_ * j] =  sw + sw * scale * (r0 * tn1 + r1 * tn2) /
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


ie_Mat Kernel::stokes_get(const std::vector<int>& I_,
                          const std::vector<int>& J_,
                          double * timing) const {
  double start, end;
  double scale = 1.0 / (M_PI);

  if (timing != nullptr) {
    start = omp_get_wtime();
  }
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (int j = 0; j < J_.size(); j++) {
    int src_ind = J_[j];
    int j_point_index = src_ind / 2;
    int j_points_vec_index = j_point_index * 2;

    double sp1 = boundary->points[j_points_vec_index];
    double sp2 =  boundary->points[j_points_vec_index + 1];
    double sn1 =  boundary->normals[j_points_vec_index];
    double sn2 = boundary->normals[j_points_vec_index + 1];
    double sw =  boundary->weights[j_point_index];
    double sc = boundary->curvatures[j_point_index];

    for (int i = 0; i < I_.size(); i++) {
      int tgt_ind = I_[i];
      int i_point_index = tgt_ind / 2;
      int i_points_vec_index = i_point_index * 2;

      double tp1 = boundary->points[i_points_vec_index];
      double tp2 = boundary->points[i_points_vec_index + 1];
      double tn1 = boundary->normals[i_points_vec_index];
      double tn2 = boundary->normals[i_points_vec_index + 1];

      if (tp1 == sp1 && sp2 == tp2) {
        double potential = - 0.5 * sc * sw * scale;

        if (tgt_ind % 2 == 0) {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] = IMPROVE_CONDITION - 0.5 + potential
                                     * sn2 * sn2 +  sw * tn1 * sn1;
          } else {
            ret.mat[i + olda_ * j] = -potential * sn1 * sn2 + sw * tn2 *
                                     sn1;
          }
        } else {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] =  -potential * sn1 * sn2 + sw * tn2 * sn1;
          } else {
            ret.mat[i + olda_ * j] = IMPROVE_CONDITION - 0.5 + potential
                                     * sn1 * sn1 + sw * tn2 * sn2;
          }
        }
      } else {
        double r0 = tp1 - sp1;
        double r1 = tp2 - sp2;
        double potential = sw * scale * (r0 * sn1 + r1 * sn2) /
                           (pow(r0 * r0 + r1 * r1, 2));
        if (tgt_ind % 2 == 0) {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] = potential * r0 * r0 + sw * tn1 * sn1;
          } else {
            ret.mat[i + olda_ * j] = potential * r0 * r1 + sw * tn1 * sn2;
          }
        } else {
          if (src_ind % 2 == 0) {
            ret.mat[i + olda_ * j] = potential * r1 * r0 + sw * tn2 * sn1;
          } else {
            ret.mat[i + olda_ * j] = potential * r1 * r1 + sw * tn2 * sn2;
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


ie_Mat Kernel::get_id_mat(const QuadTree* tree,
                          const QuadTreeNode* node,
                          bool strong_admissibility) const {
  double cntr_x = node->corners[0] + node->side_length / 2.0;
  double cntr_y = node->corners[1] + node->side_length / 2.0;

  std::vector<int> active_box = node->src_dof_lists.active_box;
  if (!strong_admissibility) {
    // Grab all points inside the proxy circle which are outside the box
    std::vector<int> inner_circle, outside_box;


    // So if we're at level 2 or 1, we don't use the proxy trick
    // If at level 1, just grab active from neighbors

    if (node->level == 1) {
      for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
        if (level_node->id != node->id) {
          for (int matrix_index :
               level_node->src_dof_lists.active_box) {
            outside_box.push_back(matrix_index);
          }
        }
      }
      ie_Mat mat(2 * outside_box.size(), active_box.size());
      mat.set_submatrix(0, outside_box.size(), 0, active_box.size(),
                        (*this)(outside_box, active_box), false, true);
      mat.set_submatrix(outside_box.size(), 2 * outside_box.size(),
                        0, active_box.size(),
                        (*this)(active_box, outside_box), true, true);
      return mat;
    }
    // If at level 2, grab active from all on level, plus from leaves of level 1
    if (node->level == 2) {
      for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
        if (level_node->id != node->id) {
          for (int matrix_index :
               level_node->src_dof_lists.active_box) {
            outside_box.push_back(matrix_index);
          }
        }
      }
      for (QuadTreeNode* level_node : tree->levels[node->level - 1]->nodes) {
        if (level_node->is_leaf) {
          for (int matrix_index :
               level_node->src_dof_lists.original_box) {
            outside_box.push_back(matrix_index);
          }
        }
      }

      ie_Mat mat(2 * outside_box.size(), active_box.size());
      mat.set_submatrix(0, outside_box.size(), 0, active_box.size(),
                        (*this)(outside_box, active_box), false, true);
      mat.set_submatrix(outside_box.size(), 2 * outside_box.size(),
                        0, active_box.size(),
                        (*this)(active_box, outside_box), true, true);
      return mat;
    }

// DEBUGGINGGGGGGGGGGGGGGGGG

    // for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
    //   if (level_node->id != node->id) {
    //     for (int matrix_index :
    //           level_node->src_dof_lists.active_box) {
    //       outside_box.push_back(matrix_index);
    //     }
    //   }
    // }
    // for (int lvl = node->level - 1; lvl > 0; lvl--) {

    //   for (QuadTreeNode* level_node : tree->levels[lvl]->nodes) {
    //     if (level_node->is_leaf) {
    //       for (int matrix_index :
    //             level_node->src_dof_lists.original_box) {
    //         outside_box.push_back(matrix_index);
    //       }
    //     }
    //   }
    // }

    // ie_Mat matt(2 * outside_box.size(), active_box.size());
    // matt.set_submatrix(0, outside_box.size(), 0, active_box.size(),
    //                   (*this)(outside_box, active_box), false, true);
    // matt.set_submatrix(outside_box.size(), 2 * outside_box.size(),
    //                   0, active_box.size(),
    //                   (*this)(active_box, outside_box), true, true);
    // return matt;


    for (int matrix_index : node->src_dof_lists.near) {
      // outside_box.push_back(matrix_index);
      int point_index = matrix_index / solution_dimension;
      int points_vec_index = point_index * domain_dimension;
      double x = tree->boundary->points[points_vec_index];
      double y = tree->boundary->points[points_vec_index + 1];
      double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
      if (dist < RADIUS_RATIO * node->side_length) {
        inner_circle.push_back(matrix_index);
      }
    }


    ie_Mat pxy = get_proxy_mat(cntr_x, cntr_y, node->side_length
                               * RADIUS_RATIO, tree, active_box);
    // Now all the matrices are gathered, put them into mat.
    ie_Mat mat(2 * inner_circle.size() + pxy.height(), active_box.size());

    double knl_str = omp_get_wtime();
    mat.set_submatrix(0, inner_circle.size(),
                      0, active_box.size(), (*this)(inner_circle, active_box),
                      false, true);
    mat.set_submatrix(inner_circle.size(), 2 * inner_circle.size(),
                      0, active_box.size(), (*this)(active_box, inner_circle),
                      true, true);
    double knl_end = omp_get_wtime();
    ie_Mat::kernel_time += (knl_end - knl_str);

    double pxy_start = omp_get_wtime();
    mat.set_submatrix(2 * inner_circle.size(),  pxy.height()
                      + 2 * inner_circle.size(), 0,
                      active_box.size(),  pxy, false, true);
    double pxy_end = omp_get_wtime();
    ie_Mat::proxy_time += (pxy_end - pxy_start);
    return mat;
  } else {
    return get_proxy_mat(cntr_x, cntr_y, node->side_length * 1.5, tree,
                         active_box);
  }
}


ie_Mat Kernel::get_proxy_mat(double cntr_x, double cntr_y,
                             double r, const QuadTree * tree,
                             const std::vector<int>& box_inds) const {
  // each row is a pxy point, cols are box dofs
  double proxy_weight = 2.0 * M_PI * r / NUM_PROXY_POINTS;
  std::vector<double>  pxy_p, pxy_n;
  for (int i = 0; i < NUM_PROXY_POINTS; i++) {
    double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
    for (int k = 1; k < 2; k++) {   // modify this for annulus
      double eps = (k - 1) * 0.001;
      pxy_p.push_back(cntr_x + (r + eps) * cos(ang));
      pxy_p.push_back(cntr_y + (r + eps) * sin(ang));
      pxy_n.push_back(cos(ang));
      pxy_n.push_back(sin(ang));
    }
  }

  switch (pde) {
    case ie_solver_config::Pde::LAPLACE:
      return laplace_proxy_get(pxy_p, pxy_n ,
                               proxy_weight, box_inds);
    case ie_solver_config::Pde::LAPLACE_NEUMANN:
      return laplace_neumann_proxy_get(pxy_p, pxy_n ,
                                       proxy_weight, box_inds);
      break;
    default:
      // case ie_solver_config::Pde::STOKES:
      return stokes_proxy_get(pxy_p, pxy_n , proxy_weight,
                              box_inds);
  }
}


ie_Mat Kernel::stokes_proxy_get(const std::vector<double> & pxy_p,
                                const std::vector<double> & pxy_n,
                                double pxy_w,
                                const std::vector<int> & box_inds) const {
  double scale = 1.0 / (M_PI);
  ie_Mat ret(2 * pxy_p.size(), box_inds.size());
  int lda = 2 * pxy_p.size();

  // Then Active to Proxy
  for (int j = 0; j < box_inds.size(); j++) {
    int src_ind = box_inds[j];
    int j_point_index = src_ind / 2;
    int j_points_vec_index = j_point_index * 2;

    double sp1 = boundary->points[j_points_vec_index];
    double sp2 =  boundary->points[j_points_vec_index + 1];
    double sn1 =  boundary->normals[j_points_vec_index];
    double sn2 = boundary->normals[j_points_vec_index + 1];
    double sw =  boundary->weights[j_point_index];

    for (int i = 0; i < pxy_p.size(); i += 2) {
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
  for (int j = 0; j < pxy_p.size(); j += 2) {
    double sp1 = pxy_p[j];
    double sp2 =  pxy_p[j + 1];
    double sn1 =  pxy_n[j];
    double sn2 = pxy_n[j + 1];
    double sw =  pxy_w;
    for (int i = 0; i < box_inds.size(); i++) {
      int tgt_ind = box_inds[i];
      int i_point_index = tgt_ind / 2;
      int i_points_vec_index = i_point_index * 2;

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



ie_Mat Kernel::laplace_proxy_get(const std::vector<double> & pxy_p,
                                 const std::vector<double> & pxy_n,
                                 double pxy_w,
                                 const std::vector<int> & box_inds) const {
  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(pxy_p.size(), box_inds.size());
  int lda = pxy_p.size();
  // First Active to Proxy
  for (int j = 0; j < box_inds.size(); j++) {
    int src_ind = box_inds[j];

    double sp1 = boundary->points[2 * src_ind];
    double sp2 =  boundary->points[2 * src_ind + 1];
    double sn1 =  boundary->normals[2 * src_ind];
    double sn2 = boundary->normals[2 * src_ind + 1];
    double sw =  boundary->weights[src_ind];

    for (int i = 0; i < pxy_p.size(); i += 2) {
      double tp1 = pxy_p[i];
      double tp2 = pxy_p[i + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;
      ret.mat[(i / 2) + lda * j] = -sw * scale * (r0 * sn1 + r1 * sn2) /
                                   (r0 * r0 + r1 * r1);
    }
  }

  // Then Proxy To Active
  for (int j = 0; j < pxy_p.size(); j += 2) {
    double sp1 = pxy_p[j];
    double sp2 =  pxy_p[j + 1];
    double sn1 =  pxy_n[j];
    double sn2 = pxy_n[j + 1];
    double sw =  pxy_w;
    for (int i = 0; i < box_inds.size(); i++) {
      int tgt_ind = box_inds[i];
      double tp1 = boundary->points[2 * tgt_ind];
      double tp2 = boundary->points[2 * tgt_ind + 1];
      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;

      ret.mat[(pxy_p.size() / 2) + i * lda + (j / 2)] = -sw * scale *
          (r0 * sn1 + r1 * sn2) /
          (r0 * r0 + r1 * r1);
    }
  }

  return ret;
}


ie_Mat Kernel::laplace_neumann_proxy_get(const std::vector<double> & pxy_p,
    const std::vector<double> & pxy_n,
    double pxy_w,
    const std::vector<int> & box_inds) const {

  double scale = 1.0 / (2 * M_PI);
  ie_Mat ret(pxy_p.size(), box_inds.size());
  int lda = pxy_p.size();
  // First Active to Proxy
  for (int j = 0; j < box_inds.size(); j++) {
    int src_ind = box_inds[j];

    double sp1 = boundary->points[2 * src_ind];
    double sp2 =  boundary->points[2 * src_ind + 1];
    double sw =  boundary->weights[src_ind];

    for (int i = 0; i < pxy_p.size(); i += 2) {
      double tp1 = pxy_p[i];
      double tp2 = pxy_p[i + 1];
      double tn1 = pxy_n[i];
      double tn2 = pxy_n[i + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;
      ret.mat[(i / 2) + lda * j] = sw + sw * scale *
                                   (r0 * tn1 + r1 * tn2) / (r0 * r0 + r1 * r1);
    }
  }

  // Then Proxy To Active
  for (int j = 0; j < pxy_p.size(); j += 2) {
    double sp1 = pxy_p[j];
    double sp2 =  pxy_p[j + 1];
    double sw =  pxy_w;

    for (int i = 0; i < box_inds.size(); i++) {
      int tgt_ind = box_inds[i];

      double tp1 = boundary->points[2 * tgt_ind];
      double tp2 = boundary->points[2 * tgt_ind + 1];
      double tn1 = boundary->normals[2 * tgt_ind];
      double tn2 = boundary->normals[2 * tgt_ind + 1];

      double r0 = tp1 - sp1;
      double r1 = tp2 - sp2;

      ret.mat[(pxy_p.size() / 2) + i * lda + (j / 2)] =   sw + sw * scale *
          (r0 * tn1 + r1 * tn2) / (r0 * r0 + r1 * r1);
    }
  }

  return ret;
}


// double Kernel::forward_get(int tgt_ind,
//                            int src_ind) const {
//   Dof tgt, src;
//   int i_point_index = tgt_ind / solution_dimension;
//   int i_points_vec_index = i_point_index * domain_dimension;
//   int j_point_index = src_ind / solution_dimension;
//   int j_points_vec_index = j_point_index * domain_dimension;

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


// double Kernel::get(int tgt_ind, int src_ind) const {
//   Dof tgt, src;
//   int i_point_index = tgt_ind / solution_dimension;
//   int i_points_vec_index = i_point_index * domain_dimension;
//   int j_point_index = src_ind / solution_dimension;
//   int j_points_vec_index = j_point_index * domain_dimension;

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

// ie_Mat Kernel::forward_get(const std::vector<int>& I_,
//                            const std::vector<int>& J_) const {
//   ie_Mat ret(I_.size(), J_.size());
//   int olda_ = I_.size();
//   for (int i = 0; i < I_.size(); i++) {
//     for (int j = 0; j < J_.size(); j++) {
//       ret.mat[i + olda_ * j] = forward_get(I_[i], J_[j]);
//     }
//   }
//   return ret;
// }

}  // namespace ie_solver
