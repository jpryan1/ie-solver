// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include "ie_solver/skel_factorization/skel_factorization.h"
#include "ie_solver/kernel/kernel.h"

namespace ie_solver {

// TODO(John) URGENTLY skelnear needs to be replaced as a variable name

// TODO(John) here and elsewhere check on passing by reference

SkelFactorization::SkelFactorization(double id_tol, bool strong_admissibility_,
                                     int solution_dimension_,
                                     int domain_dimension_) {
  assert(id_tol > 0 && "id_tol must be greater than one to init tools.");
  this->id_tol = id_tol;
  strong_admissibility = strong_admissibility_;
  solution_dimension = solution_dimension_;
  domain_dimension = domain_dimension_;
  id_time = 0.;
  make_id_total = 0.;
}


// TODO(John) the permutation vector is weirdly unique among the
// src_dof_lists and is weird to deal with in perturbing, can we just ditch
// it somehow?
int SkelFactorization::id_compress(const Kernel& kernel,
                                   const QuadTree* tree, QuadTreeNode* node) {
  assert(node != nullptr && "InterpolativeDecomposition fails on null node.");
  assert(node->src_dof_lists.active_box.size() > 0 &&
         "Num of DOFs must be positive in InterpolativeDecomposition.");
  // TODO(John) better variable name
  ie_Mat pxy;
  double make_id_start = omp_get_wtime();
  make_id_mat(kernel, &pxy, tree, node);
  double make_id_end = omp_get_wtime();
  make_id_total += make_id_end - make_id_start;

  if (pxy.height() == 0) {
    return 0;
  }
  std::vector<unsigned int> p;
  double id_start = omp_get_wtime();
  unsigned int numskel = pxy.id(&p, &node->T, id_tol);
  double id_end = omp_get_wtime();
  id_time += (id_end - id_start);
  if (numskel == 0) {
    node->compression_ratio = 0.;
    return 0;
  }
  node->compression_ratio = (node->T.width() /
                             (0.0 + node->T.width() + node->T.height()));
  node->src_dof_lists.set_rs_ranges(p, node->T.height(), node->T.width());
  node->src_dof_lists.set_skelnear_range(strong_admissibility);

  return node->T.width();
}


void SkelFactorization::get_x_matrices(ie_Mat* K, const ie_Mat& Z, ie_Mat* Xrr,
                                       const std::vector<unsigned int>& r,
                                       const std::vector<unsigned int>& s,
                                       const std::vector<unsigned int>& n) {
  assert(r.size()*s.size() > 0 &&
         "For get_x_matrices, there must be redundant and skel dofs.");

  unsigned int r_size = r.size();
  unsigned int s_size = s.size();
  unsigned int n_size = n.size();
  ie_Mat  Xrs(r_size, s_size), Xsr(s_size, r_size), Xrn(r_size, n_size),
          Xnr(n_size, r_size);

  *Xrr = (*K)(r, r);
  Xrs = (*K)(r, s);
  Xsr = (*K)(s, r);
  if (n_size > 0) {
    Xrn = (*K)(r, n);
    Xnr = (*K)(n, r);
  }

  ie_Mat::gemm(TRANSPOSE, NORMAL, -1., Z, (*K)(s, r), 1., Xrr);
  ie_Mat::gemm(TRANSPOSE, NORMAL, -1., Z, (*K)(s, s), 1., &Xrs);
  ie_Mat::gemm(NORMAL,    NORMAL, -1., Xrs,        Z,          1., Xrr);
  ie_Mat::gemm(NORMAL,    NORMAL, -1., (*K)(s, s), Z,          1., &Xsr);
  if (n_size > 0) {
    ie_Mat::gemm(TRANSPOSE, NORMAL, -1., Z, (*K)(s, n), 1., &Xrn);
    ie_Mat::gemm(NORMAL,    NORMAL, -1., (*K)(n, s), Z,          1., &Xnr);
  }

  K->set_submatrix(r, s, Xrs);
  K->set_submatrix(s, r, Xsr);
  if (n_size > 0) {
    K->set_submatrix(r, n, Xrn);
    K->set_submatrix(n, r, Xnr);
  }
  // TODO(John) this seems like an awful lot of stores, can we avoid this?
}


// TODO(John) These involved multiplying into a buffer, then copying the buffer
// into K
// Can we do this in place? Or use just one buffer and resize each time?
// Does it matter?
// Also, we might not need to do ALL of these matmuls, since some of the
// blocks will be eliminated anyways. This should be worked out on a whiteboard
void SkelFactorization::schur_update(const Kernel& kernel, QuadTreeNode* node) {
  assert(node != nullptr && "SchurUpdate fails on null node.");
  assert(node->T.height()*node->T.width() > 0 &&
         "Z must have positive dimensions in SchurUpdate.");
  assert(node->src_dof_lists.active_box.size() > 0 &&
         "Num of DOFs must be positive in SchurUpdate.");
  // height of Z is number of skeleton columns
  unsigned int num_redundant = node->T.width();
  unsigned int num_skel      = node->T.height();
  // GENERATE K_BN,BN
  std::vector<unsigned int> BN;
  for (unsigned int idx : node->src_dof_lists.active_box) {
    BN.push_back(idx);
  }
  if (strong_admissibility) {
    for (unsigned int idx : node->src_dof_lists.near) BN.push_back(idx);
  }
  // Note that BN has all currently deactivated DoFs removed.

  ie_Mat update(BN.size(), BN.size());
  get_all_schur_updates(&update, BN, node, strong_admissibility);

  ie_Mat K_BN = kernel(BN, BN) - update;
  // Generate various index ranges within BN
  std::vector<unsigned int> s, r, n, sn;
  for (unsigned int i = 0; i < num_skel; i++) {
    s.push_back(node->src_dof_lists.permutation[i]);
    sn.push_back(node->src_dof_lists.permutation[i]);
  }
  for (unsigned int i = 0; i < num_redundant; i++) {
    r.push_back(node->src_dof_lists.permutation[i + num_skel]);
  }

  if (strong_admissibility) {
    for (unsigned int i = 0; i < node->src_dof_lists.near.size(); i++) {
      n.push_back(i + num_redundant + num_skel);
      sn.push_back(i + num_redundant + num_skel);
    }
  }
  node->X_rr = ie_Mat(num_redundant, num_redundant);
  get_x_matrices(&K_BN, node->T, &(node->X_rr), r, s, n);

  // Generate left and right schur complement matrices
  // TODO(John) change this variable naming
  int num_skelnear = sn.size();

  node->L = ie_Mat(num_skelnear,  num_redundant);
  node->U = ie_Mat(num_redundant, num_skelnear);

  ie_Mat K_BN_sn_r = K_BN(sn, r);
  ie_Mat K_BN_r_sn = K_BN(r, sn);

  node->X_rr.right_multiply_inverse(K_BN_sn_r, &node->L);
  node->X_rr.left_multiply_inverse(K_BN_r_sn, &node->U);

  node->schur_update = ie_Mat(sn.size(), sn.size());
  ie_Mat::gemm(NORMAL, NORMAL, 1.0, node->L, K_BN_r_sn, 0.,
               &(node->schur_update));

  node->compressed = true;
}


void SkelFactorization::skeletonize(const Kernel& kernel, QuadTree* tree) {
  double skel_start = omp_get_wtime();
  int node_counter = 0;
  int already_marked = 0;
  int just_skelled = 0;
  unsigned int lvls = tree->levels.size();
  int active_dofs = tree->boundary->points.size() / 2;

  for (unsigned int level = lvls - 1; level > 0; level--) {
    if (lvls - level > LEVEL_CAP) {
      break;
    }
    // std::cout << "Skelling level " << level << std::endl;

    tree->remove_inactive_dofs_at_level(level);
    QuadTreeLevel* current_level = tree->levels[level];

    for (unsigned int n = 0; n < current_level->nodes.size(); n++) {
      tree->remove_inactive_dofs_at_level(level);

      node_counter++;
      if (node_counter > NODE_CAP) {
        break;
      }
      QuadTreeNode* current_node = current_level->nodes[n];
      if (current_node->compressed) {
        already_marked++;
        continue;
      }
      if (current_node->src_dof_lists.active_box.size()
          < MIN_DOFS_TO_COMPRESS) {
        continue;
      }

      int redundants = id_compress(kernel, tree, current_node);

      active_dofs -= redundants;
      if (redundants == 0) {
        continue;
      }
      just_skelled++;
      schur_update(kernel, current_node);
    }
  }
  // If the above breaks due to a cap, we need to manually propagate active
  // boxes up the tree.
  tree->remove_inactive_dofs_at_all_boxes();

  std::vector<unsigned int> allskel = tree->root->src_dof_lists.active_box;
  if (allskel.size() > 0) {
    ie_Mat allskel_updates = ie_Mat(allskel.size(), allskel.size());
    get_all_schur_updates(&allskel_updates, allskel, tree->root, false);
    allskel_mat = kernel(allskel, allskel) - allskel_updates;
  }
  std::cout << "skel_count: already_marked/newly_skelled: " << already_marked <<
            " " << just_skelled << std::endl;
  // check_factorization_against_kernel(kernel, tree);
  double skel_end = omp_get_wtime();
  std::cout << "timing: skeletonize " << (skel_end - skel_start) << std::endl;
  std::cout << "timing: all_interpdecomps " << id_time << std::endl;

  std::cout << "timing: all_make_id_mats " << make_id_total << std::endl;
}



void SkelFactorization::get_all_schur_updates(ie_Mat* updates,
    const std::vector<unsigned int>& BN, const QuadTreeNode* node,
    bool get_neighbors) const {
  assert(node != nullptr && "get_all_schur_updates fails on null node.");
  assert(BN.size() > 0 && "get_all_schur_updates needs positive num of DOFs");
  if (!node->is_leaf) get_descendents_updates(updates, BN, node);

  if (get_neighbors) {
    for (QuadTreeNode* neighbor : node->neighbors) {
      if (neighbor->level != node->level) continue;
      if (neighbor->compressed) get_update(updates, BN, neighbor);
      if (!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);
    }
  }
}


void SkelFactorization::get_descendents_updates(ie_Mat* updates,
    const std::vector<unsigned int>& BN, const QuadTreeNode* node)  const {
  assert(node != nullptr && "get_descendents_updates fails on null node.");
  assert(!node->is_leaf &&
         "get_descendents_updates must be called on non-leaf.");

  // by assumption, node is not a leaf
  for (QuadTreeNode* child : node->children) {
    if (child->compressed) get_update(updates, BN, child);
    if (!child->is_leaf) get_descendents_updates(updates, BN, child);
  }
}


void SkelFactorization::get_update(ie_Mat* update,
                                   const std::vector<unsigned int>& BN,
                                   const QuadTreeNode* node)  const {
  // node needs to check all its dofs against BN, enter interactions into
  // corresponding locations
  // node only updated its own BN dofs, and the redundant ones are no longer
  // relevant, so we only care about child's SN dofs
  // First create a list of Dofs that are also in node's skelnear,
  // and with each one give the index in skelnear and the index in BN
  std::vector<unsigned int> BN_;
  std::vector<unsigned int> sn_;
  for (unsigned int sn_idx = 0; sn_idx < node->src_dof_lists.skelnear.size();
       sn_idx++) {
    for (unsigned int bn_idx = 0; bn_idx < BN.size(); bn_idx++) {
      if (BN[bn_idx] == node->src_dof_lists.skelnear[sn_idx]) {
        sn_.push_back(sn_idx);
        BN_.push_back(bn_idx);
      }
    }
  }
  // For every pair of dofs shared by both, update their interaction
  int num_shared_by_both = BN_.size();
  for (unsigned int i = 0; i < num_shared_by_both; i++) {
    for (unsigned int j = 0; j < num_shared_by_both; j++) {
      update->addset(BN_[i], BN_[j], node->schur_update.get(sn_[i], sn_[j]));
    }
  }
}

void SkelFactorization::make_id_mat(const Kernel& kernel, ie_Mat* mat,
                                    const QuadTree* tree,
                                    const QuadTreeNode* node) {
  double cntr_x = node->corners[0] + node->side_length / 2.0;
  double cntr_y = node->corners[1] + node->side_length / 2.0;

  std::vector<unsigned int> active_box = node->src_dof_lists.active_box;
  if (!strong_admissibility) {
    // Grab all points inside the proxy circle which are outside the box
    std::vector<unsigned int> inner_circle, outside_box;
    int num_outside_circle = 0;
    for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
      if (level_node->id != node->id) {
        for (unsigned int matrix_index :
             level_node->src_dof_lists.active_box) {
          outside_box.push_back(matrix_index);
          unsigned int point_index = matrix_index / solution_dimension;
          unsigned int points_vec_index = point_index * domain_dimension;
          double x = tree->boundary->points[points_vec_index];
          double y = tree->boundary->points[points_vec_index + 1];
          double dist = sqrt(pow(cntr_x - x, 2) + pow(cntr_y - y, 2));
          if (dist < RADIUS_RATIO * node->side_length) {
            inner_circle.push_back(matrix_index);
          } else {
            num_outside_circle++;
          }
        }
      }
    }
    // We only do a proxy circle if there are lots outside the circle
    if (num_outside_circle < active_box.size() * 4) {
      // Now all the matrices are gathered, put them into *mat.
      *mat = ie_Mat(2 * outside_box.size(), active_box.size());
      mat->set_submatrix(0, outside_box.size(), 0, active_box.size(),
                         kernel(outside_box, active_box));
      mat->set_submatrix(outside_box.size(), 2 * outside_box.size(),
                         0, active_box.size(), kernel(active_box, outside_box), true);
    } else {
      // Construct mat of interactions with pxy circle points
      ie_Mat proxy = ie_Mat(solution_dimension * 2 * NUM_PROXY_POINTS,
                            active_box.size());

      make_proxy_mat(kernel, &proxy, cntr_x, cntr_y, node->side_length
                     * RADIUS_RATIO, tree, active_box);

      // Now all the matrices are gathered, put them into *mat.
      *mat = ie_Mat(2 * inner_circle.size() + solution_dimension * 2 *
                    NUM_PROXY_POINTS, active_box.size());
      mat->set_submatrix(0, inner_circle.size(),
                         0, active_box.size(), kernel(inner_circle, active_box));
      mat->set_submatrix(inner_circle.size(), 2 * inner_circle.size(),
                         0, active_box.size(), kernel(active_box, inner_circle), true);
      mat->set_submatrix(2 * inner_circle.size(),  solution_dimension * 2 *
                         NUM_PROXY_POINTS + 2 * inner_circle.size(), 0,
                         active_box.size(), proxy);

    }


  } else {
    *mat = ie_Mat(solution_dimension * 2 * NUM_PROXY_POINTS, active_box.size());
    make_proxy_mat(kernel, mat, cntr_x, cntr_y, node->side_length * 1.5, tree,
                   active_box);
  }
}


void SkelFactorization::make_proxy_mat(const Kernel & kernel, ie_Mat * pxy,
                                       double cntr_x, double cntr_y,
                                       double r, const QuadTree * tree,
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

      unsigned int matrix_index = box_inds[j_];
      unsigned int point_index = matrix_index / solution_dimension;
      unsigned int points_vec_index = point_index * domain_dimension;
      b.point = Vec2(tree->boundary->points[points_vec_index],
                     tree->boundary->points[points_vec_index + 1]);
      b.normal = Vec2(tree->boundary->normals[points_vec_index],
                      tree->boundary->normals[points_vec_index + 1]);
      b.curvature = tree->boundary->curvatures[point_index];
      b.weight = tree->boundary->weights[point_index];

      ie_Mat ab_tensor = kernel.get(a, b);
      ie_Mat ba_tensor = kernel.get(b, a);

      for (int k = 0; k < solution_dimension; k++) {
        pxy->set(solution_dimension * i + k, j_,
                 ab_tensor.get(k, box_inds[j_] % solution_dimension));
        pxy->set(solution_dimension * (i + NUM_PROXY_POINTS) + k, j_,
                 ba_tensor.get(box_inds[j_] % solution_dimension, k));
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////


// Sets vec(b) = vec(b) + mat*vec(a)
void SkelFactorization::apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
    const std::vector<unsigned int>& a,
    const std::vector<unsigned int>& b,
    bool transpose = false) const {
  if (a.size()*b.size() == 0) return;
  if (transpose) {
    assert(mat.height() == a.size());
  } else {
    assert(mat.width() == a.size());
  }
  ie_Mat product(b.size(),  vec->width());

  if (transpose) {
    ie_Mat::gemm(TRANSPOSE, NORMAL, 1., mat, (*vec)(a, 0, vec->width()), 0.,
                 &product);
  } else {
    ie_Mat::gemm(NORMAL,  NORMAL, 1., mat, (*vec)(a, 0, vec->width()), 0.,
                 &product);
  }
  vec->set_submatrix(b, 0, vec->width(), product + (*vec)(b, 0, vec->width()));
}


// Sets vec(range) = mat * vec(range)
void SkelFactorization::apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
    const std::vector<unsigned int>& range)
const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(), vec->width());
  ie_Mat::gemm(NORMAL, NORMAL, 1., mat, (*vec)(range, 0, vec->width()), 0.,
               &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void SkelFactorization::apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat* vec,
    const std::vector<unsigned int>& range)
const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(),  vec->width());
  mat.left_multiply_inverse((*vec)(range,  0, vec->width()), &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void SkelFactorization::apply_diag_pinv_matrix(const ie_Mat& mat, ie_Mat* vec,
    const std::vector<unsigned int>& range)
const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(),  vec->width());
  mat.left_multiply_pseudoinverse((*vec)(range,  0, vec->width()), &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void SkelFactorization::sparse_matvec(const QuadTree& quadtree, const ie_Mat& x,
                                      ie_Mat* b) const {
  *b = x;
  int lvls = quadtree.levels.size();
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->compressed) {
        continue;
      }
      // First we need to apply L_T inverse
      // L_T inverse changes the skel elements - it makes them equal to
      // T times the redundant elements + the skeleton elements.
      apply_sweep_matrix(current_node->T, b,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skel,
                         false);

      // Next we need to apply U inverse
      // U inverse changes the redundant elements - it makes them equal to
      // L transpose times the skelnear elements + the redundant elements
      apply_sweep_matrix(current_node->U, b,
                         current_node->src_dof_lists.skelnear,
                         current_node->src_dof_lists.redundant, false);
    }
  }

  // This can go through the tree in any order, is parallelizable
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->compressed) {
        continue;
      }
      apply_diag_matrix(current_node->X_rr, b,
                        current_node->src_dof_lists.redundant);
    }
  }
  // We need all of the skeleton indices. This is just the negation of
  // [0,b.size()] and the redundant DoFs
  // with that in mind...

  std::vector<unsigned int> allskel =
    quadtree.root->src_dof_lists.active_box;

  if (allskel.size() > 0) {
    apply_diag_matrix(allskel_mat, b, allskel);
  }

  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    // TODO(John) record in notes and explore the following observation:
    // changing the order of this for loop affects the accuracy of the
    // sparse_mat_vec on a random vector, BUT NOT THE SOLUTION ERROR
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      // Next we need to apply L inverse
      // L inverse changes the skelnear elements - it makes them equal to
      // L times the redundant elements + the skelnear elements
      apply_sweep_matrix(current_node->L, b,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear, false);
      // Finally we need to apply U_T inverse
      // U_T inverse changes the redundant elements - it makes them equal
      // to T transpose times the skeleton elements + the redundant
      // elements
      apply_sweep_matrix(current_node->T, b,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
    }
  }
}


void SkelFactorization::solve(const QuadTree& quadtree, ie_Mat* x,
                              const ie_Mat& b) const {
  assert(x->height() == b.height());
  int lvls = quadtree.levels.size();
  *x = b;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->compressed) {
        continue;
      }
      // Next we need to apply L inverse
      // L inverse changes the skelnear elements - it makes them equal to
      // L times the redundant elements + the skelnear elements

      // Finally we need to apply U_T inverse
      // U_T inverse changes the redundant elements - it makes them equal
      // to T transpose times the skeleton elements + the redundant
      // elements
      apply_sweep_matrix(-current_node->T, x,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, x,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear, false);
    }
  }
  // This can go through the tree in any order, is parallelizable

  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->compressed) {
        continue;
      }
      // double cond = current_node->X_rr.condition_number();
      // if (cond > 1000) {
      //   std::cout << "Node " << current_node->id << " solve X_rr inv -- ";
      //   std::cout << "Inverting w/ condition number " << cond << std::endl;
      // }
      apply_diag_inv_matrix(current_node->X_rr, x,
                            current_node->src_dof_lists.redundant);
    }
  }

  // We need all of the skeleton indices. This is just the negation of
  // [0,b.size()] and the redundant DoFs
  // with that in mind...
  std::vector<unsigned int> allskel = quadtree.root->src_dof_lists.active_box;
  // quadtree.allskel_mat.write_singular_values_to_file("allskel_sing_vals.txt");
  if (allskel.size() > 0) {
    // double cond2 = allskel_mat.condition_number();
    // if (cond2 > 1e6) {
    //   std::cout << "Allskel inv -- ";
    //   std::cout << "Inverting w/ condition number " << cond2 << std::endl;
    // }
    apply_diag_inv_matrix(allskel_mat, x, allskel);
  }
  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      // for(unsigned int n = 0; n < current_level->nodes.size(); n++){

      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      // First we need to apply L_T inverse
      // L_T inverse changes the skel elements - it makes them equal to T
      // times the redundant elements + the skeleton elements.
      // Next we need to apply U inverse
      // U inverse changes the redundant elements - it makes them equal to
      // L transpose times the skelnear elements + the redundant elements
      apply_sweep_matrix(-current_node->U, x,
                         current_node->src_dof_lists.skelnear,
                         current_node->src_dof_lists.redundant, false);
      apply_sweep_matrix(-current_node->T, x,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skel,
                         false);
    }
  }
}


std::vector<unsigned int> big_to_small(const std::vector<unsigned int>& big,
                                       const std::unordered_map<unsigned int,
                                       unsigned int>& map) {
  std::vector<unsigned int> small;
  for (unsigned int idx : big) {
    small.push_back(map.at(idx));
  }
  return small;
}

// TODO make usage of psi height and u width more readable
void SkelFactorization::multiply_connected_solve(const QuadTree& quadtree,
    ie_Mat* mu, ie_Mat* alpha,
    const ie_Mat& b) const {
  // TODO(John) LaTeX this up into a more readable form
  assert(x->height() == b.height());
  int lvls = quadtree.levels.size();
  
  std::vector<unsigned int> allskel = quadtree.root->src_dof_lists.active_box;
  std::vector<unsigned int> sorted_allskel = allskel;
  std::sort(sorted_allskel.begin(), sorted_allskel.end());
  unsigned int skel_idx = 0;
  
  std::vector<unsigned int> allredundant;
  for (unsigned int i = 0;
       i < quadtree.boundary->weights.size()*solution_dimension; i++) {
    if (skel_idx < sorted_allskel.size() && i == sorted_allskel[skel_idx]) {
      skel_idx++;
    } else {
      allredundant.push_back(i);
    }
  }
  // In our bordered linear system, the skel and redundant indices are partitioned
  // so we create a map from their original index into their partition
  std::unordered_map<unsigned int, unsigned int> skel_big2small, red_big2small;
  for (int i = 0; i < allskel.size(); i++) {
    skel_big2small[allskel[i]] = i;
  }
  for (int i = 0; i < allredundant.size(); i++) {
    red_big2small[allredundant[i]] = i;
  }

  *mu = b;

  ie_Mat modified_Psi = Psi.transpose();
  ie_Mat modified_U = U;

  // First apply the sweep matrices to x and U to modify them.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->compressed) {
        continue;
      }

      apply_sweep_matrix(-current_node->T, mu,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, mu,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear, false);

      apply_sweep_matrix(-current_node->T, &modified_U,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, &modified_U,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear, false);
    }
  }
  
  ie_Mat w = (*mu)(allredundant, 0, 1);
  ie_Mat z = (*mu)(allskel, 0, 1);
  
  // Now apply the other sweep matrices to Psi to modify it.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, &modified_Psi,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant,
                         true);
      apply_sweep_matrix(-current_node->U, &modified_Psi,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear,
                         true);
    }
  }
  modified_Psi = modified_Psi.transpose();

  ie_Mat Dinv_w = w;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->compressed) {
        continue;
      }
      std::vector<unsigned int> small_redundants = big_to_small(
            current_node->src_dof_lists.redundant, red_big2small);

      apply_diag_inv_matrix(current_node->X_rr, &Dinv_w,
                            small_redundants);
    }
  }

  // B has modified_Psi as a subblock, but is otherwise a ton of 0s, so
  // we don't actually form B.
  ie_Mat B_Dinv_w_nonzero(modified_Psi.height(), 1);

  ie_Mat::gemm(NORMAL, NORMAL, 1., modified_Psi(0, modified_Psi.height(),
               allredundant), Dinv_w, 0., &B_Dinv_w_nonzero);

  // B_Dinv_w.set_submatrix(allskel.size(), allskel.size() + Psi.height(), 0, 1,
  //                       B_Dinv_w_nonzero);
  // The rest of the matrix is automatically zeros, as should be the case.
  
  ie_Mat M(allskel.size() + Psi.height(), 1);

  M.set_submatrix(0, allskel.size(), 0, 1, z);
  M.set_submatrix(allskel.size(), M.height(), 0, 1, -B_Dinv_w_nonzero);

  // Again, C is mostly 0s, so we just apply Dinv to the nonzero block
  ie_Mat Dinv_C_nonzero = modified_U;
  for (int level = lvls - 1; level >= 0; level--) {  // level>=0; level--){
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->compressed) {
        continue;
      }
      std::vector<unsigned int> small_redundants = big_to_small(
            current_node->src_dof_lists.redundant, red_big2small);

      apply_diag_inv_matrix(current_node->X_rr, &Dinv_C_nonzero,
                            small_redundants);
    }
  }

  ie_Mat B_Dinv_C_nonzero(modified_Psi.height(), modified_U.width());
  ie_Mat::gemm(NORMAL, NORMAL, 1., modified_Psi, Dinv_C_nonzero, 0., &B_Dinv_C_nonzero);

  ie_Mat S(allskel.size() + Psi.height(),
           allskel.size() + Psi.height());
                  
  S.set_submatrix(0, allskel.size(),
                  0, allskel.size(), allskel_mat);
  S.set_submatrix(allskel.size(), S.height(),
                  0, allskel.size(),
                  modified_Psi(0, Psi.height(), allskel));
  S.set_submatrix(0, allskel.size(),
                  allskel.size(), S.width(),
                  modified_U(allskel, 0, U.width()));
  ie_Mat ident(Psi.height(), Psi.height());
  ident.eye(Psi.height());
  S.set_submatrix(allskel.size(), S.height(), allskel.size(), S.width(),
                  - ident - B_Dinv_C_nonzero);

  ie_Mat y(S.height(), 1);
  S.left_multiply_inverse(M, &y);

  *alpha = y(allskel.size(), y.height(), 0, 1);

  ie_Mat Cy(modified_U.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., modified_U, *alpha, 0., &Cy);
  
  ie_Mat N = w - Cy;
  // b_vec.set_submatrix(allskel.size(), b_vec.height(), 0, 1,
  //                     b_vec(allskel.size(), b_vec.height(), 0, 1) - Cx);

  // ie_Mat second_paren = b_vec;
  // ie_Mat y = second_paren;
  ie_Mat Dinv_N = N;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->compressed) {
        continue;
      }
      std::vector<unsigned int> small_redundants = big_to_small(
            current_node->src_dof_lists.redundant, red_big2small);

      apply_diag_inv_matrix(current_node->X_rr, &Dinv_N,
                            small_redundants);
    }
  }

  mu->set_submatrix(allredundant, 0, 1, Dinv_N);
  mu->set_submatrix(allskel, 0, 1, y(0, allskel.size(), 0, 1));

  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      // for(unsigned int n = 0; n < current_level->nodes.size(); n++){

      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      // First we need to apply L_T inverse
      // L_T inverse changes the skel elements - it makes them equal to T
      // times the redundant elements + the skeleton elements.
      // Next we need to apply U inverse
      // U inverse changes the redundant elements - it makes them equal to
      // L transpose times the skelnear elements + the redundant elements
      apply_sweep_matrix(-current_node->U, mu,
                         current_node->src_dof_lists.skelnear,
                         current_node->src_dof_lists.redundant, false);
      apply_sweep_matrix(-current_node->T, mu,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skel,
                         false);
    }
  }
}

}  // namespace ie_solver
