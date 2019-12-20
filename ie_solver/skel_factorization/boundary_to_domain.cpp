// // Copyright 2019 John Paul Ryan
// #include <iostream>
// #include <cassert>
// #include <cmath>
// #include "ie-solver/tools/ie_solver_tools.h"

// namespace ie_solver {


// void IeSolverTools::b2d_apply_diag_matrix(const ie_Mat& mat,
//     const std::vector<int>& tgt,
//     const std::vector<int>& src,
//     const ie_Mat& vec_in, ie_Mat* vec_out) {

//   if (tgt.size() * src.size() == 0) {
//     return;
//   }
//   std::vector<int> ZERO_VECTOR;
//   ZERO_VECTOR.push_back(0);

//   ie_Mat product(tgt.size(), 1);
//   ie_Mat::gemv(NORMAL, 1., mat, vec_in(src, ZERO_VECTOR), 0., &product);
//   vec_out->set_submatrix(tgt, ZERO_VECTOR,
//                         product + (*vec_out)(tgt, ZERO_VECTOR));
// }


// void IeSolverTools::b2dsparse_matvec(const Kernel& K, const QuadTree& tree,
//                                     const ie_Mat& x,
//                                     ie_Mat* b) {
//   int lvls = tree.levels.size();

//   for (int level_ = lvls - 2; level_ > 0; level_--) {
//     ie_Mat add(b->height(), 1);
//     ie_Mat start = x;

//     // T _src inv
//     for (int level = lvls - 1; level > level_; level--) {
//       QuadTreeLevel* current_level = tree.levels[level];
//       for (QuadTreeNode* current_node : current_level->nodes) {
//         if (!current_node->compressed) continue;
//         apply_sweep_matrix(current_node->src_T, &start,
//                           current_node->src_dof_lists.redundant,
//                           current_node->src_dof_lists.skel,
//                           false);
//       }
//     }

//     // K_r
//     for (int level = lvls - 1; level > level_; level--) {
//       QuadTreeLevel* current_level = tree.levels[level];
//       for (QuadTreeNode* current_node : current_level->nodes) {
//         if (!current_node->compressed) {
//           continue;
//         }
//         if (current_node->X_rs.height() * current_node->X_rs.width() > 0) {
//           b2d_apply_diag_matrix(current_node->X_rs,
//                                 current_node->tgt_dof_lists.redundant,
//                                 current_node->src_dof_lists.skel,
//                                 start, &add);
//         }
//         if (current_node->X_rr.height() * current_node->X_rr.width() > 0) {
//           b2d_apply_diag_matrix(current_node->X_rr,
//                                 current_node->tgt_dof_lists.redundant,
//                                 current_node->src_dof_lists.redundant,
//                                 start, &add);
//         }
//         if (current_node->X_sr.height() * current_node->X_sr.width() > 0) {
//           b2d_apply_diag_matrix(current_node->X_sr,
//                                 current_node->tgt_dof_lists.skel,
//                                 current_node->src_dof_lists.redundant,
//                                 start, &add);
//         }
//       }
//     }

//     // T_tgt inv
//     for (int level = level_ + 1; level < lvls; level++) {
//       QuadTreeLevel* current_level = tree.levels[level];
//       for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
//         QuadTreeNode* current_node = current_level->nodes[n];
//         if (!current_node->compressed) continue;
//         apply_sweep_matrix(current_node->tgt_T, &add,
//                           current_node->tgt_dof_lists.skel,
//                           current_node->tgt_dof_lists.redundant, true);
//       }
//     }
//     (*b) += add;
//   }

//   ie_Mat final_start = x;
//   for (int level = lvls - 1; level >= 0; level--) {
//     QuadTreeLevel* current_level = tree.levels[level];
//     for (QuadTreeNode* current_node : current_level->nodes) {
//       if (!current_node->compressed) {
//         continue;
//       }
//       apply_sweep_matrix(current_node->src_T, &final_start,
//                         current_node->src_dof_lists.redundant,
//                         current_node->src_dof_lists.skel,
//                         false);
//     }
//   }
// // We need all of the skeleton indices. This is just the negation of
// // [0,b.size()] and the redundant DoFs
// // with that in mind...

//   ie_Mat final_add(b->height(), 1);

//   QuadTreeLevel* current_level = tree.levels[0];
//   for (QuadTreeNode* current_node : current_level->nodes) {
//     if (!current_node->compressed) {
//       continue;
//     }
//     if (current_node->X_rs.height() * current_node->X_rs.width() > 0) {
//       b2d_apply_diag_matrix(current_node->X_rs,
//                             current_node->tgt_dof_lists.redundant,
//                             current_node->src_dof_lists.skel,
//                             final_start, &final_add);
//     }
//     if (current_node->X_rr.height() * current_node->X_rr.width() > 0) {
//       b2d_apply_diag_matrix(current_node->X_rr,
//                             current_node->tgt_dof_lists.redundant,
//                             current_node->src_dof_lists.redundant,
//                             final_start, &final_add);
//     }
//     if (current_node->X_sr.height() * current_node->X_sr.width() > 0) {
//       b2d_apply_diag_matrix(current_node->X_sr,
//                             current_node->tgt_dof_lists.skel,
//                             current_node->src_dof_lists.redundant,
//                             final_start, &final_add);
//     }
//   }

//   std::vector<int> src_allskel =
//     tree.root->src_dof_lists.active_box;

//   std::vector<int> tgt_allskel =
//     tree.root->tgt_dof_lists.active_box;

//   if (src_allskel.size()*tgt_allskel.size() > 0) {
//     ie_Mat allskel_mat = K.forward_get(tgt_allskel, src_allskel);
//     b2d_apply_diag_matrix(allskel_mat,
//                           tgt_allskel,
//                           src_allskel,
//                           final_start, &final_add);
//   }

//   for (int level = 0; level < lvls; level++) {
//     QuadTreeLevel* current_level = tree.levels[level];
//     // TODO(John) record in notes and explore the following observation:
//     // changing the order of this for loop affects the accuracy of the
//     // sparse_mat_vec on a random vector, BUT NOT THE SOLUTION ERROR
//     for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
//       QuadTreeNode* current_node = current_level->nodes[n];
//       if (!current_node->compressed) {
//         continue;
//       }

//       // Finally we need to apply U_T inverse
//       // U_T inverse changes the redundant elements - it makes them equal
//       // to T transpose times the skeleton elements + the redundant
//       // elements
//       apply_sweep_matrix(current_node->tgt_T, &final_add,
//                         current_node->tgt_dof_lists.skel,
//                         current_node->tgt_dof_lists.redundant, true);
//     }
//   }

//   (*b) += final_add;

//   for (int i = 0; i < tree.domain_points.size(); i += 2) {
//     if (!tree.boundary->is_in_domain(
//           Vec2(tree.domain_points[i], tree.domain_points[i + 1]))) {
//       for (int j = 0; j < solution_dimension; j++) {
//         b->set((i / 2) * solution_dimension + j, 0, 0.0);
//       }
//     }
//   }
// }


// void IeSolverTools::b2dskeletonize(const Kernel& kernel, QuadTree* tree) {
//   int lvls = tree->levels.size();
//   for (int level = lvls - 1; level > 4; level--) {
//     QuadTreeLevel* current_level = tree->levels[level];
//     // First, get all active dofs from children
//     for (QuadTreeNode * node : current_level->nodes) {
//       if (node->compressed) continue;
//       remove_inactive_dofs_at_box(node);
//     }
//     // Next, get all active near dofs from neighbors
//     for (QuadTreeNode* node_a : current_level->nodes) {
//       if (node_a->compressed) continue;
//       node_a->src_dof_lists.near.clear();
//       node_a->tgt_dof_lists.near.clear();
//       for (QuadTreeNode* neighbor : node_a->neighbors) {
//         // Some neighbors are smaller boxes from higher levels, we don't
//         // care about those, their parents have the updated information.
//         if (neighbor->level > node_a->level) {
//           continue;
//         }
//         for (int idx : neighbor->src_dof_lists.active_box) {
//           node_a->src_dof_lists.near.push_back(idx);
//         }
//         for (int idx : neighbor->tgt_dof_lists.active_box) {
//           node_a->tgt_dof_lists.near.push_back(idx);
//         }
//       }
//     }
//     for (int n = 0; n < current_level->nodes.size(); n++) {
//       QuadTreeNode* current_node = current_level->nodes[n];
//       if (current_node->compressed ||
//           current_node->src_dof_lists.active_box.size() == 0 ||
//           current_node->tgt_dof_lists.active_box.size() == 0
//         ) {
//         continue;
//       }
//       if (current_node->src_dof_lists.active_box.size() < MIN_DOFS_TO_COMPRESS
//           && current_node->tgt_dof_lists.active_box.size()
//           < MIN_DOFS_TO_COMPRESS) {
//         continue;
//       }
//       b2dinterpolative_decomposition(kernel, tree, current_node);
//     }
//   }
//   remove_inactive_dofs_at_all_boxes(tree);
// }


// int IeSolverTools::b2dinterpolative_decomposition(const Kernel& kernel,
//     const QuadTree* tree, QuadTreeNode* node) {
//   assert(node != nullptr && "InterpolativeDecomposition fails on null node.");
//   assert(node->src_dof_lists.active_box.size() +
//         node->tgt_dof_lists.active_box.size() > 0 &&
//         "Num of DOFs must be positive in InterpolativeDecomposition.");

//   // SRC TO TGT
//   ie_Mat tgt_pxy;
//   make_tgt_id_mat(kernel, &tgt_pxy, tree, node);
//   if (tgt_pxy.width() * tgt_pxy.height() != 0) {
//     std::vector<int> tgt_p;
//     int tgt_numskel = tgt_pxy.id(&tgt_p, &node->tgt_T, id_tol);

//     if (tgt_numskel == 0) return 0;
//     set_rs_ranges(&node->tgt_dof_lists, tgt_p, node->tgt_T.height(),
//                   node->tgt_T.width());
//     set_skelnear_range(&node->tgt_dof_lists);
//   } else {
//     node->tgt_T = ie_Mat(0, 0);
//   }

//   // TGT TO SRC
//   ie_Mat src_pxy;
//   make_src_id_mat(kernel, &src_pxy, tree, node);
//   if (src_pxy.width() * src_pxy.height() != 0) {
//     std::vector<int> src_p;
//     int src_numskel = src_pxy.id(&src_p, &node->src_T, id_tol);
//     if (src_numskel == 0) return 0;
//     set_rs_ranges(&node->src_dof_lists, src_p, node->src_T.height(),
//                   node->src_T.width());
//     set_skelnear_range(&node->src_dof_lists);
//   } else {
//     node->src_T = ie_Mat(0, 0);
//   }

// // get K matrices and set them!
//   std::vector<int> tgt_r = node->tgt_dof_lists.redundant;
//   std::vector<int> tgt_s = node->tgt_dof_lists.skel;
//   std::vector<int> src_r = node->src_dof_lists.redundant;
//   std::vector<int> src_s = node->src_dof_lists.skel;

//   ie_Mat Xrr = kernel.forward_get(tgt_r, src_r);
//   ie_Mat Xrs = kernel.forward_get(tgt_r, src_s);
//   ie_Mat Xsr = kernel.forward_get(tgt_s, src_r);

//   ie_Mat Ksr = kernel.forward_get(tgt_s, src_r);
//   ie_Mat Krs = kernel.forward_get(tgt_r, src_s);
//   ie_Mat Kss = kernel.forward_get(tgt_s, src_s);


// // TODO(John) i think we're doing some unnecessary matmul here.
//   if (node->tgt_T.height() * node->tgt_T.width()*node->src_T.height() *
//       node->src_T.width() > 0) {
//     ie_Mat::gemm(TRANSPOSE, NORMAL, -1., node->tgt_T, Ksr, 1., &Xrr);
//     ie_Mat::gemm(TRANSPOSE, NORMAL, -1., node->tgt_T, Kss, 1., &Xrs);

//     ie_Mat::gemm(NORMAL,    NORMAL, -1., Xrs,     node->src_T,       1., &Xrr);
//     ie_Mat::gemm(NORMAL,    NORMAL, -1., Kss, node->src_T,       1., &Xsr);
//   }

//   node->X_rr = Xrr;
//   node->X_rs = Xrs;
//   node->X_sr = Xsr;
//   node->compressed = true;

//   return node->src_T.width() + node->tgt_T.width();
// }


// void IeSolverTools::make_tgt_id_mat(const Kernel& kernel, ie_Mat* mat,
//                                     const QuadTree* tree,
//                                     const QuadTreeNode* node) {
//   // Interactions where
//   //    TGT:  Box domain points
//   //    SRC:  Far Boundary Points = inner circle boundary + pxy circle

//   double cntr_x = node->corners[0] + node->side_length / 2.0;
//   double cntr_y = node->corners[1] + node->side_length / 2.0;
//   double radius_ratio = 1.5;
//   double r = node->side_length * radius_ratio;

//   if (true) {  // node->level <= tree->no_proxy_level){
//     // No proxy circle
//     std::vector<int> far;
//     for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
//       if (level_node->id != node->id) {
//         for (int matrix_index : level_node->src_dof_lists.active_box) {
//           far.push_back(matrix_index);
//         }
//       }
//     }

//     *mat = ie_Mat(far.size(),
//                   node->tgt_dof_lists.active_box.size());
//     ie_Mat box_far = kernel.forward_get(node->tgt_dof_lists.active_box, far);
//     for (int i = 0; i < box_far.height(); i++) {
//       for (int j = 0; j < box_far.width(); j++) {
//         mat->set(j, i, box_far.get(i, j));
//       }
//     }

//     return;
//   }



//   // Grab all points inside the proxy circle
//   std::vector<int> inner_circle;
//   for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
//     if (level_node->id != node->id) {
//       for (int matrix_index : level_node->src_dof_lists.active_box) {
//         int point_index = matrix_index / solution_dimension;
//         int points_vec_index = point_index * domain_dimension;

//         double x = tree->boundary->points[points_vec_index];
//         double y = tree->boundary->points[points_vec_index + 1];
//         double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
//         if (dist < radius_ratio * node->side_length) {
//           inner_circle.push_back(matrix_index);
//         }
//       }
//     }
//   }

//   *mat = ie_Mat(inner_circle.size() + solution_dimension * NUM_PROXY_POINTS,
//                 node->tgt_dof_lists.active_box.size());

//   ie_Mat near_box = kernel.forward_get(node->tgt_dof_lists.active_box,
//                                       inner_circle);

//   for (int i = 0; i < near_box.height(); i++) {
//     for (int j = 0; j < near_box.width(); j++) {
//       mat->set(j, i, near_box.get(i, j));
//     }
//   }

//   double proxy_weight = 2.0 * M_PI * r / NUM_PROXY_POINTS;
//   double proxy_curvature = 1.0 / r;

//   for (int i = 0; i < NUM_PROXY_POINTS; i++) {
//     double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
//     Vec2 p(cntr_x + r * cos(ang), cntr_y + r * sin(ang));
//     for (int j = 0; j < node->tgt_dof_lists.active_box.size(); j++) {
//       Dof a, b;
//       a.is_boundary = true;
//       a.point = p;
//       a.normal = Vec2(cos(ang), sin(ang));
//       a.curvature = proxy_curvature;
//       a.weight = proxy_weight;

//       int matrix_index = node->tgt_dof_lists.active_box[j];
//       int point_index = matrix_index / solution_dimension;
//       int points_vec_index = point_index * domain_dimension;
//       b.point = Vec2(tree->domain_points[points_vec_index],
//                     tree->domain_points[points_vec_index + 1]);
//       b.is_boundary = false;
//       ie_Mat ba_tensor = kernel.get(b, a);

//       for (int k = 0; k < solution_dimension; k++) {
//         mat->set(inner_circle.size() + solution_dimension * i + k, j,
//                 ba_tensor.get(node->tgt_dof_lists.active_box[j]
//                               % solution_dimension, k));
//       }
//     }
//   }
// }


// void IeSolverTools::make_src_id_mat(const Kernel& kernel, ie_Mat* mat,
//                                     const QuadTree* tree,
//                                     const QuadTreeNode* node) {
//   // Interactions where
//   //    TGT:  Far Domain Points = Pxy mat points + inner circle domain points
//   //    SRC:  Box Boundary Points

//   double cntr_x = node->corners[0] + node->side_length / 2.0;
//   double cntr_y = node->corners[1] + node->side_length / 2.0;
//   double radius_ratio = 1.5;
//   double r = node->side_length * radius_ratio;

//   if (true) {  // node->level <= tree->no_proxy_level){
//     // No proxy circle
//     std::vector<int> far;
//     for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
//       if (level_node->id != node->id) {
//         for (int matrix_index : level_node->tgt_dof_lists.active_box) {
//           far.push_back(matrix_index);
//         }
//       }
//     }
//     *mat = ie_Mat(far.size(), node->src_dof_lists.active_box.size());
//     ie_Mat far_box = kernel.forward_get(far, node->src_dof_lists.active_box);
//     for (int i = 0; i < far_box.height(); i++) {
//       for (int j = 0; j < far_box.width(); j++) {
//         mat->set(i, j, far_box.get(i, j));
//       }
//     }
//     return;
//   }


//   // Grab all points inside the proxy circle
//   std::vector<int> inner_circle;
//   for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
//     if (level_node->id != node->id) {
//       for (int matrix_index : level_node->tgt_dof_lists.active_box) {
//         int point_index = matrix_index / solution_dimension;
//         int points_vec_index = point_index * domain_dimension;

//         double x = tree->domain_points[points_vec_index];
//         double y = tree->domain_points[points_vec_index + 1];
//         double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
//         if (dist < radius_ratio * node->side_length) {
//           inner_circle.push_back(matrix_index);
//         }
//       }
//     }
//   }

//   *mat = ie_Mat(inner_circle.size() + solution_dimension * NUM_PROXY_POINTS,
//                 node->src_dof_lists.active_box.size());

//   ie_Mat near_box = kernel.forward_get(inner_circle,
//                                       node->src_dof_lists.active_box);

//   for (int i = 0; i < near_box.height(); i++) {
//     for (int j = 0; j < near_box.width(); j++) {
//       mat->set(i, j, near_box.get(i, j));
//     }
//   }

//   for (int i = 0; i < NUM_PROXY_POINTS; i++) {
//     double ang = 2 * M_PI * i * (1.0 / NUM_PROXY_POINTS);
//     Vec2 p(cntr_x + r * cos(ang), cntr_y + r * sin(ang));
//     for (int j = 0; j < node->src_dof_lists.active_box.size(); j++) {
//       Dof a, b;
//       a.point = p;
//       a.is_boundary = false;
//       int matrix_index = node->src_dof_lists.active_box[j];
//       int point_index = matrix_index / solution_dimension;
//       int points_vec_index = point_index * domain_dimension;
//       b.point = Vec2(tree->boundary->points[points_vec_index],
//                     tree->boundary->points[points_vec_index + 1]);
//       b.normal = Vec2(tree->boundary->normals[points_vec_index],
//                       tree->boundary->normals[points_vec_index + 1]);
//       b.curvature = tree->boundary->curvatures[point_index];
//       b.weight = tree->boundary->weights[point_index];
//       b.is_boundary = true;
//       ie_Mat ab_tensor = kernel.get(a, b);
//       for (int k = 0; k < solution_dimension; k++) {
//         mat->set(inner_circle.size() + solution_dimension * i + k, j,
//                 ab_tensor.get(k, node->src_dof_lists.active_box[j]
//                               % solution_dimension));
//       }
//     }
//   }
// }

// }  // namespace ie_solver
