// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/tools/ie_solver_tools.h"
#include "ie-solver/kernel.h"
#define NUM_PROXY_POINTS 128
namespace ie_solver {

void IeSolverTools::make_id_mat(const Kernel& kernel, ie_Mat* mat,
                                const QuadTree* tree,
                                const QuadTreeNode* node) {
  if (solution_dimension == 2) {
    make_stokes_id_mat(kernel, mat, tree, node);
    return;
  }
  double cntr_x = node->corners[0] + node->side_length / 2.0;
  double cntr_y = node->corners[1] + node->side_length / 2.0;
  double radius_ratio = 1.5;

  if (!strong_admissibility) {
    // Grab all points inside the proxy circle
    std::vector<unsigned int> inner_circle;
    for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
      if (level_node->id != node->id) {
        for (unsigned int idx : level_node->interaction_lists.active_box) {
          double x, y;
          x = tree->boundary->points[2 * idx];
          y = tree->boundary->points[2 * idx + 1];
          double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
          if (dist < radius_ratio * node->side_length) {
            inner_circle.push_back(idx);
          }
        }
      }
    }

    ie_Mat near_box = kernel(inner_circle,
                             node->interaction_lists.active_box);
    ie_Mat box_near(inner_circle.size(),
                    node->interaction_lists.active_box.size());
    kernel(node->interaction_lists.active_box,
           inner_circle).transpose_into(&box_near);

    ie_Mat proxy = ie_Mat(2 * NUM_PROXY_POINTS,
                          node->interaction_lists.active_box.size());

    make_proxy_mat(kernel, &proxy, cntr_x, cntr_y, node->side_length
                   * radius_ratio, tree, node->interaction_lists.active_box);

    *mat = ie_Mat(2 * inner_circle.size() + 2 * NUM_PROXY_POINTS,
                  node->interaction_lists.active_box.size());
    mat->set_submatrix(0, inner_circle.size(),
                       0, node->interaction_lists.active_box.size(),
                       near_box);
    mat->set_submatrix(inner_circle.size(), 2 * inner_circle.size(),
                       0, node->interaction_lists.active_box.size(),
                       box_near);
    mat->set_submatrix(2 * inner_circle.size(),
                       2 * NUM_PROXY_POINTS + 2 * inner_circle.size(),
                       0, node->interaction_lists.active_box.size(), proxy);
  } else {
    *mat = ie_Mat(2 * NUM_PROXY_POINTS,
                  node->interaction_lists.active_box.size());
    make_proxy_mat(kernel, mat, cntr_x, cntr_y, node->side_length * 1.5, tree,
                   node->interaction_lists.active_box);
  }
}


void IeSolverTools::make_proxy_mat(const Kernel& kernel, ie_Mat* pxy,
                                   double cntr_x, double cntr_y,
                                   double r, const QuadTree* tree,
                                   const std::vector<unsigned int>& box_inds) {
  // each row is a pxy point, cols are box dofs
  double proxy_weight = 2.0 * M_PI * r / NUM_PROXY_POINTS;
  double proxy_curvature = 1.0 / r;

  for (int i = 0; i < NUM_PROXY_POINTS; i++) {
    double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
    Vec2 p(cntr_x + r * cos(ang), cntr_y + r * sin(ang));
    for (unsigned int j_ = 0; j_ < box_inds.size(); j_++) {
      Dof a, b;
      a.point = p;
      a.normal = Vec2(cos(ang), sin(ang));
      a.curvature = proxy_curvature;
      a.weight = proxy_weight;

      unsigned int j = box_inds[j_];
      b.point = Vec2(tree->boundary->points[2 * j],
                     tree->boundary->points[2 * j + 1]);
      b.normal = Vec2(tree->boundary->normals[2 * j],
                      tree->boundary->normals[2 * j + 1]);
      b.curvature = tree->boundary->curvatures[j];
      b.weight = tree->boundary->weights[j];

      pxy->set(i, j_, kernel.laplace_kernel(a, b));
      pxy->set(i + NUM_PROXY_POINTS, j_, kernel.laplace_kernel(b, a));
    }
  }
}


void IeSolverTools::make_stokes_id_mat(const Kernel& kernel, ie_Mat* mat,
                                       const QuadTree* tree,
                                       const QuadTreeNode* node) {
  double cntr_x = node->corners[0] + node->side_length / 2.0;
  double cntr_y = node->corners[1] + node->side_length / 2.0;
  double radius_ratio = 1.5;
  if (!strong_admissibility) {
    // Grab all points inside the proxy circle
    std::vector<unsigned int> inner_circle;
    for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
      if (level_node->id != node->id) {
        for (unsigned int idx : level_node->interaction_lists.active_box) {
          double x, y;
          idx = 2 * (idx / 2);
          x = tree->boundary->points[idx];
          y = tree->boundary->points[idx + 1];
          double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
          if (dist < radius_ratio * node->side_length) {
            inner_circle.push_back(idx);
          }
        }
      }
    }

    ie_Mat near_box = kernel(inner_circle,
                             node->interaction_lists.active_box);
    ie_Mat box_near(inner_circle.size(),
                    node->interaction_lists.active_box.size());
    kernel(node->interaction_lists.active_box,
           inner_circle).transpose_into(&box_near);

    ie_Mat proxy = ie_Mat(2 * 2 * NUM_PROXY_POINTS,
                          node->interaction_lists.active_box.size());

    make_stokes_proxy_mat(kernel, &proxy, cntr_x, cntr_y, node->side_length
                          * radius_ratio, tree,
                          node->interaction_lists.active_box);
    *mat = ie_Mat(2 * inner_circle.size() + 2 * 2 * NUM_PROXY_POINTS,
                  node->interaction_lists.active_box.size());
    mat->set_submatrix(0, inner_circle.size(),
                       0, node->interaction_lists.active_box.size(),
                       near_box);
    mat->set_submatrix(inner_circle.size(), 2 * inner_circle.size(),
                       0, node->interaction_lists.active_box.size(),
                       box_near);
    mat->set_submatrix(2 * inner_circle.size(),
                       2 * 2 * NUM_PROXY_POINTS + 2 * inner_circle.size(),
                       0, node->interaction_lists.active_box.size(), proxy);

  } else {
    *mat = ie_Mat(2 * 2 * NUM_PROXY_POINTS,
                  node->interaction_lists.active_box.size());
    make_stokes_proxy_mat(kernel, mat, cntr_x, cntr_y, node->side_length * 1.5,
                          tree, node->interaction_lists.active_box);
  }
}


void IeSolverTools::make_stokes_proxy_mat(const Kernel& kernel, ie_Mat* pxy,
    double cntr_x, double cntr_y,
    double r, const QuadTree* tree,
    const std::vector<unsigned int>& box_inds) {
  // each row is a pxy point, cols are box dofs
  double proxy_weight = 2.0 * M_PI * r / NUM_PROXY_POINTS;
  double proxy_curvature = 1.0 / r;

  for (int i = 0; i < NUM_PROXY_POINTS; i++) {
    double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
    Vec2 p(cntr_x + r * cos(ang), cntr_y + r * sin(ang));
    // TODO(John) better variable names
    for (unsigned int j_ = 0; j_ < box_inds.size(); j_++) {
      // For influence from dof onto proxy, need col of tensor
      // For influence from proxy into dof, need row, transpose

      Dof a, b;
      a.point = p;
      a.normal = Vec2(cos(ang), sin(ang));
      a.curvature = proxy_curvature;
      a.weight = proxy_weight;

      unsigned int j = 2 * (box_inds[j_] / 2);
      b.point = Vec2(tree->boundary->points[j],
                     tree->boundary->points[j + 1]);
      b.normal = Vec2(tree->boundary->normals[j],
                      tree->boundary->normals[j + 1]);
      b.curvature = tree->boundary->curvatures[j / 2];
      b.weight = tree->boundary->weights[j / 2];

      ie_Mat ab_tensor = kernel.stokes_kernel(a, b);
      ie_Mat ba_tensor = kernel.stokes_kernel(b, a);

      if (box_inds[j_] % 2 == 0) {
        pxy->set(2 * i, j_, ab_tensor.get(0, 0));
        pxy->set(2 * i + 1, j_, ab_tensor.get(1, 0));
        pxy->set(2 * i + 2 * NUM_PROXY_POINTS, j_, ba_tensor.get(0, 0));
        pxy->set(2 * i + 1 + 2 * NUM_PROXY_POINTS, j_, ab_tensor.get(0, 1));
      } else {
        pxy->set(2 * i, j_, ab_tensor.get(0, 1));
        pxy->set(2 * i + 1, j_, ab_tensor.get(1, 1));
        pxy->set(2 * i + 2 * NUM_PROXY_POINTS, j_, ba_tensor.get(1, 0));
        pxy->set(2 * i + 1 + 2 * NUM_PROXY_POINTS, j_, ab_tensor.get(1, 1));
      }
    }
  }
}

}  // namespace ie_solver
