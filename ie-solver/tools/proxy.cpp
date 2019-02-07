// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "ie-solver/tools/ie_solver_tools.h"
#include "ie-solver/kernel.h"
#define NUM_PROXY_POINTS 128
namespace ie_solver {

void IeSolverTools::make_id_mat(const Kernel& K, ie_Mat* mat,
                                const QuadTree* tree,
                                const QuadTreeNode* node) {
  double cntr_x = node->corners[0] + node->side_length / 2.0;
  double cntr_y = node->corners[1] + node->side_length / 2.0;
  double radius_ratio = 1.5;
  if (is_stokes) {
    // // TODO fix this later, match laplace please, make cleaner please
    // mat = ie_Mat(200, node->interaction_lists.active_box.size());
    // make_stokes_proxy_mat(mat, cntr_x, cntr_y, node->side_length*2,
    //   node->interaction_lists.active_box );
  } else {
    if (!strong_admissibility) {
      std::vector<unsigned int> inner_circle;

      for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
        if (level_node->id != node->id) {
          for (unsigned int idx : level_node->interaction_lists.active_box) {
            double x = tree->boundary->points[2 * idx];
            double y = tree->boundary->points[2 * idx + 1];
            double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
            if (dist < radius_ratio * node->side_length) {
              inner_circle.push_back(idx);
            }
          }
        }
      }
      // inner_circle = node->interaction_lists.near;
      ie_Mat box_near = K(inner_circle,
                          node->interaction_lists.active_box);

      ie_Mat proxy = ie_Mat(NUM_PROXY_POINTS,
                            node->interaction_lists.active_box.size());
      make_proxy_mat(&proxy, cntr_x, cntr_y, node->side_length * radius_ratio,
                     tree, node->interaction_lists.active_box);
      *mat = ie_Mat(inner_circle.size() + NUM_PROXY_POINTS,
                    node->interaction_lists.active_box.size());
      mat->set_submatrix(0, inner_circle.size(),
                         0, node->interaction_lists.active_box.size(),
                         box_near);
      mat->set_submatrix(inner_circle.size(),
                         NUM_PROXY_POINTS + inner_circle.size(),
                         0, node->interaction_lists.active_box.size(), proxy);
    } else {
      *mat = ie_Mat(NUM_PROXY_POINTS,
                    node->interaction_lists.active_box.size());
      make_proxy_mat(mat, cntr_x, cntr_y, node->side_length * 1.5, tree,
                     node->interaction_lists.active_box);
    }
  }
}


void IeSolverTools::make_proxy_mat(ie_Mat* pxy, double cntr_x, double cntr_y,
                                   double r, const QuadTree* tree,
                                   const std::vector<unsigned int>& box_inds) {
  // each row is a pxy point, cols are box dofs
  double scale = 1.0 / (2 * M_PI);

  for (int i = 0; i < NUM_PROXY_POINTS; i++) {
    double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
    Vec2 p(cntr_x + r * cos(ang), cntr_y + r * sin(ang));
    for (unsigned int j = 0; j < box_inds.size(); j++) {
      unsigned int box_index = box_inds[j];
      Vec2 q(tree->boundary->points[2 * box_index],
             tree->boundary->points[2 * box_index + 1]);

      Vec2 r = p - q;
      Vec2 n(tree->boundary->normals[2 * box_index],
             tree->boundary->normals[2 * box_index + 1]);

      double potential = -tree->boundary->weights[box_index] * scale *
                         (r.dot(n)) / (r.dot(r));
      pxy->set(i, j, potential);
    }
  }
}


void IeSolverTools::make_stokes_proxy_mat(ie_Mat* pxy, double cntr_x,
    double cntr_y, double r, const std::vector<unsigned int>& box_indices) {

  // double scale = 1.0 / (M_PI);
  // for(int i = 0; i < 200; i += 2){
  //   double ang = 2 * M_PI * i * 0.01;
  //   Vec2 p(cntr_x + r*cos(ang), cntr_y + r*sin(ang));
  //   for(unsigned int j = 0; j < box_indices.size(); j++){

  //     unsigned int box_index = box_indices[j];
  //     Vec2 q(points[2*box_index], points[2*box_index+1]);

  //     Vec2 r = p-q;
  //     Vec2 n(normals[2*box_index], normals[2*box_index+1]);
  //     double r0 = r.a[0];
  //     double r1 = r.a[1];

  //     double potential = weights[box_index] * scale *
  //       (r.dot(n)) / (pow(r.dot(r),2));

  //     if( box_index % 2 == 0 ){
  //       pxy.set(i  , j, potential*r0*r0);
  //       pxy.set(i+1, j, potential*r0*r1);
  //     }else{
  //       pxy.set(i  , j, potential*r0*r1);
  //       pxy.set(i+1, j, potential*r1*r1);
  //     }
  //   }
  // }
}

}  // namespace ie_solver
