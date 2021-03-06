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

std::vector<int> big_to_small(const std::vector<int>& big,
                              const std::unordered_map<int,
                              int>& map) {
  std::vector<int> small;
  for (int idx : big) {
    small.push_back(map.at(idx));
  }
  return small;
}


SkelFactorization::SkelFactorization(double id_tol,
                                     bool strong_admissibility_, int nt) {
  assert(id_tol > 0 && "id_tol must be greater than one to init tools.");
  this->id_tol = id_tol;
  strong_admissibility = strong_admissibility_;
  num_threads = nt;
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

  double make_start = omp_get_wtime();
  ie_Mat pxy = kernel.get_id_mat(tree, node, strong_admissibility);
  double make_end = omp_get_wtime();
  make_mat_time += (make_end - make_start);
  if (pxy.height() == 0) {
    return 0;
  }
  std::vector<int> p;
  double id_start = omp_get_wtime();
  int numskel = pxy.id(&p, &node->T, id_tol);
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


// TODO(John) These involved multiplying into a buffer, then copying the buffer
// into K
// Can we do this in place? Or use just one buffer and resize each time?
// Does it matter?
// Also, we might not need to do ALL of these matmuls, since some of the
// blocks will be eliminated anyways. This should be worked out on a whiteboard
void SkelFactorization::decouple(const Kernel& kernel, QuadTreeNode* node) {

  // height of Z is number of skeleton columns
  int num_redundant = node->T.width();
  int num_skel      = node->T.height();
  // GENERATE K_BN,BN
  std::vector<int> BN;
  for (int idx : node->src_dof_lists.active_box) {
    BN.push_back(idx);
  }
  if (strong_admissibility) {
    for (int idx : node->src_dof_lists.near) BN.push_back(idx);
  }
  // Note that BN has all currently deactivated DoFs removed.

  ie_Mat update(BN.size(), BN.size());
  get_all_schur_updates(&update, BN, node, strong_admissibility);

  ie_Mat K_BN = kernel(BN, BN) - update;
  // Generate various index ranges within BN
  std::vector<int> s, r, n, sn;
  for (int i = 0; i < num_skel; i++) {
    s.push_back(node->src_dof_lists.permutation[i]);
    sn.push_back(node->src_dof_lists.permutation[i]);
  }
  for (int i = 0; i < num_redundant; i++) {
    r.push_back(node->src_dof_lists.permutation[i + num_skel]);
  }

  if (strong_admissibility) {
    for (int i = 0; i < node->src_dof_lists.near.size(); i++) {
      n.push_back(i + num_redundant + num_skel);
      sn.push_back(i + num_redundant + num_skel);
    }
  }

  ie_Mat K_BN_r_sn = K_BN(r, s) - node->T.transpose() * K_BN(s, s) ;
  node->X_rr = K_BN(r, r) - node->T.transpose() * K_BN(s, r)
               - K_BN_r_sn * node->T;
  ie_Mat K_BN_sn_r = K_BN(s, r) - K_BN(s, s) * node->T;

  node->X_rr.LU_factorize(&node->X_rr_lu, &node->X_rr_piv);
  node->X_rr_is_LU_factored = true;

  node->L = ie_Mat(sn.size(),  num_redundant);
  node->U = ie_Mat(num_redundant, sn.size());

  node->X_rr_lu.right_multiply_inverse(K_BN_sn_r, node->X_rr_piv, &node->L);
  node->X_rr_lu.left_multiply_inverse(K_BN_r_sn, node->X_rr_piv,  &node->U);
  node->schur_update = node->L * K_BN_r_sn;
  node->compressed = true;
}


// void report_num_threads(int level) {
//   #pragma omp single
//   {
//     printf("Level %d: number of threads in the team - %d\n",
//     level, omp_get_num_threads());
//   }
// }


void SkelFactorization::skeletonize(const Kernel& kernel, QuadTree* tree) {
  double skel_start = omp_get_wtime();
  int node_counter = 0;
  int lvls = tree->levels.size();
  // int active_dofs = tree->boundary->points.size() / 2;
  // make_mat_time = 0.0;
  // id_time = 0.0;
  // schur_time = 0;
  ie_Mat::proxy_time = 0.;
  ie_Mat::kernel_time = 0.;
  for (int level = lvls - 1; level >  0; level--) {
    tree->remove_inactive_dofs_at_level(level);
    QuadTreeLevel* current_level = tree->levels[level];
    #pragma omp parallel for num_threads(num_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      double node_start = omp_get_wtime();
      QuadTreeNode* current_node = current_level->nodes[n];
      if (current_node->compressed || current_node->src_dof_lists.active_box.size()
          < MIN_DOFS_TO_COMPRESS) {
        continue;
      }
      if (id_compress(kernel, tree, current_node) == 0) {
        continue;
      }
      decouple(kernel, current_node);
      double node_end = omp_get_wtime();
      current_node->compress_time = node_end - node_start;
      node_counter++;
    }
  }
  // If the above breaks due to a cap, we need to manually propagate active
  // boxes up the tree.
  tree->remove_inactive_dofs_at_all_boxes();

  std::vector<int> allskel = tree->root->src_dof_lists.active_box;
  if (allskel.size() > 0) {
    ie_Mat allskel_updates = ie_Mat(allskel.size(), allskel.size());
    get_all_schur_updates(&allskel_updates, allskel, tree->root, false);
    tree->allskel_mat = kernel(allskel, allskel) - allskel_updates;
  }

  if (tree->U.width() == 0) {
    tree->allskel_mat.LU_factorize(&tree->allskel_mat_lu, &tree->allskel_mat_piv);
    double skel_end = omp_get_wtime();
    // std::cout << "timing: skeletonize " << (skel_end - skel_start) << std::endl;
    return;
  }

  std::vector<QuadTreeNode*> all_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
  }

  std::vector<int> sorted_allskel = allskel;
  std::sort(sorted_allskel.begin(), sorted_allskel.end());
  int skel_idx = 0;

  std::vector<int> allredundant;
  for (int i = 0;
       i < tree->boundary->weights.size()* tree->solution_dimension; i++) {
    if (skel_idx < sorted_allskel.size() && i == sorted_allskel[skel_idx]) {
      skel_idx++;
    } else {
      allredundant.push_back(i);
    }
  }
  if (allredundant.size() == 0) {
    std::cout << "No compression possible" << std::endl;
    exit(0);
  }

  // In our bordered linear system, the skel and redundant indices are
  // partitioned so we create a map from their original index into their
  // partition
  std::unordered_map<int, int> skel_big2small, red_big2small;
  for (int i = 0; i < allskel.size(); i++) {
    skel_big2small[allskel[i]] = i;
  }
  for (int i = 0; i < allredundant.size(); i++) {
    red_big2small[allredundant[i]] = i;
  }

  ie_Mat modified_Psi = tree->Psi.transpose();
  ie_Mat modified_U = tree->U;

  // First apply the sweep matrices to x and U to modify them.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];

    #pragma omp parallel for num_threads(num_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, &modified_U,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, &modified_U,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear, false);
    }
  }

  // Now apply the other sweep matrices to Psi to modify it.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];
    #pragma omp parallel for num_threads(num_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
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

  // Again, C is mostly 0s, so we just apply Dinv to the nonzero block
  ie_Mat Dinv_C_nonzero = modified_U(allredundant, 0, modified_U.width());

  #pragma omp parallel for num_threads(num_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];

    if (current_node->src_dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    std::vector<int> small_redundants = big_to_small(
                                          current_node->src_dof_lists.redundant, red_big2small);
    assert(current_node->X_rr_is_LU_factored);
    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv,
                          &Dinv_C_nonzero,
                          small_redundants);
  }

  ie_Mat ident(tree->Psi.height(), tree->Psi.height());
  ident.eye(tree->Psi.height());
  ie_Mat S(allskel.size() + tree->Psi.height(),
           allskel.size() + tree->Psi.height());

  S.set_submatrix(0, allskel.size(),
                  0, allskel.size(), tree->allskel_mat);
  S.set_submatrix(allskel.size(), S.height(),
                  0, allskel.size(),
                  modified_Psi(0, tree->Psi.height(), allskel));
  S.set_submatrix(0, allskel.size(),
                  allskel.size(), S.width(),
                  modified_U(allskel, 0, tree->U.width()));

  S.set_submatrix(allskel.size(), S.height(), allskel.size(), S.width(),
                  - ident - (modified_Psi(0, modified_Psi.height(),
                                          allredundant) * Dinv_C_nonzero));

  S.LU_factorize(&tree->S_LU, &tree->S_piv);

  double skel_end = omp_get_wtime();
  // std::cout << "timing: skeletonize " << (skel_end - skel_start) << std::endl;

}


// void SkelFactorization::diag_block_factorizer() {
//   while (true) {
//     while (!kill_factorizer && block_to_factorize == nullptr) {}
//     if (kill_factorizer) {
//       return;
//     }
//   }
// }


void SkelFactorization::get_all_schur_updates(ie_Mat * updates,
    const std::vector<int>& BN, const QuadTreeNode * node,
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


void SkelFactorization::get_descendents_updates(ie_Mat * updates,
    const std::vector<int>& BN, const QuadTreeNode * node)  const {
  assert(node != nullptr && "get_descendents_updates fails on null node.");
  assert(!node->is_leaf &&
         "get_descendents_updates must be called on non-leaf.");

  // by assumption, node is not a leaf
  for (QuadTreeNode* child : node->children) {
    if (child->compressed) get_update(updates, BN, child);
    if (!child->is_leaf) get_descendents_updates(updates, BN, child);
  }
}


void SkelFactorization::get_update(ie_Mat * update,
                                   const std::vector<int>& BN,
                                   const QuadTreeNode * node)  const {
  // node needs to check all its dofs against BN, enter interactions into
  // corresponding locations
  // node only updated its own BN dofs, and the redundant ones are no longer
  // relevant, so we only care about child's SN dofs
  // First create a list of Dofs that are also in node's skelnear,
  // and with each one give the index in skelnear and the index in BN
  std::vector<int> BN_;
  std::vector<int> sn_;
  for (int sn_idx = 0; sn_idx < node->src_dof_lists.skelnear.size();
       sn_idx++) {
    for (int bn_idx = 0; bn_idx < BN.size(); bn_idx++) {
      if (BN[bn_idx] == node->src_dof_lists.skelnear[sn_idx]) {
        sn_.push_back(sn_idx);
        BN_.push_back(bn_idx);
      }
    }
  }
  // For every pair of dofs shared by both, update their interaction
  int num_shared_by_both = BN_.size();
  for (int i = 0; i < num_shared_by_both; i++) {
    for (int j = 0; j < num_shared_by_both; j++) {
      update->addset(BN_[i], BN_[j], node->schur_update.get(sn_[i], sn_[j]));
    }
  }
}


///////////////////////////////////////////////////////////////////////////////


// Sets vec(b) = vec(b) + mat*vec(a)
void SkelFactorization::apply_sweep_matrix(const ie_Mat & mat, ie_Mat * vec,
    const std::vector<int>& a,
    const std::vector<int>& b,
    bool transpose = false) const {
  if (a.size()*b.size() == 0) return;
  if (transpose) {
    assert(mat.height() == a.size());
  } else {
    assert(mat.width() == a.size());
  }
  ie_Mat product;

  if (transpose) {
    product = mat.transpose() * (*vec)(a, 0, vec->width());
  } else {
    product = mat * (*vec)(a, 0, vec->width());
  }
  vec->set_submatrix(b, 0, vec->width(), product + (*vec)(b, 0, vec->width()));
}


// Sets vec(range) = mat * vec(range)
void SkelFactorization::apply_diag_matrix(const ie_Mat & mat, ie_Mat * vec,
    const std::vector<int>& range)
const {
  if (range.size() == 0) return;
  vec->set_submatrix(range,  0, vec->width(),  mat * (*vec)(range, 0,
                     vec->width()));
}


void SkelFactorization::apply_diag_inv_matrix(const ie_Mat & mat,
    const std::vector<lapack_int>& piv, ie_Mat * vec,
    const std::vector<int>& range) const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(),  vec->width());
  mat.left_multiply_inverse((*vec)(range,  0, vec->width()), piv, &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void SkelFactorization::sparse_matvec(const QuadTree & quadtree,
                                      const ie_Mat & x,
                                      ie_Mat * b) const {
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

  std::vector<int> allskel =
    quadtree.root->src_dof_lists.active_box;

  if (allskel.size() > 0) {
    apply_diag_matrix(quadtree.allskel_mat, b, allskel);
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


void SkelFactorization::solve(const QuadTree & quadtree, ie_Mat * x,
                              const ie_Mat & b) const {
  assert(x->height() == b.height());
  int lvls = quadtree.levels.size();
  *x = b;
  std::vector<QuadTreeNode*> all_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
  }
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(num_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
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

  #pragma omp parallel for num_threads(num_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];
    if (current_node->src_dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    assert(current_node->X_rr_is_LU_factored);

    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv, x,
                          current_node->src_dof_lists.redundant);
  }


  // We need all of the skeleton indices. This is just the negation of
  // [0,b.size()] and the redundant DoFs
  // with that in mind...
  std::vector<int> allskel = quadtree.root->src_dof_lists.active_box;
  if (allskel.size() > 0) {
    apply_diag_inv_matrix(quadtree.allskel_mat_lu, quadtree.allskel_mat_piv, x,
                          allskel);
  }
  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(num_threads)
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      // for(int n = 0; n < current_level->nodes.size(); n++){

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


// TODO(John) make usage of psi height and u width more readable
void SkelFactorization::multiply_connected_solve(const QuadTree & quadtree,
    ie_Mat * mu, ie_Mat * alpha,
    const ie_Mat & b) const {

  // TODO(John) LaTeX this up into a more readable form
  assert(mu->height() == b.height());
  int lvls = quadtree.levels.size();
  std::vector<QuadTreeNode*> all_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
  }

  std::vector<int> allskel = quadtree.root->src_dof_lists.active_box;
  std::vector<int> sorted_allskel = allskel;
  std::sort(sorted_allskel.begin(), sorted_allskel.end());
  int skel_idx = 0;

  std::vector<int> allredundant;
  for (int i = 0;
       i < quadtree.boundary->weights.size()* quadtree.solution_dimension; i++) {
    if (skel_idx < sorted_allskel.size() && i == sorted_allskel[skel_idx]) {
      skel_idx++;
    } else {
      allredundant.push_back(i);
    }
  }
  // In our bordered linear system, the skel and redundant indices are
  // partitioned so we create a map from their original index into their
  // partition
  std::unordered_map<int, int> skel_big2small, red_big2small;
  for (int i = 0; i < allskel.size(); i++) {
    skel_big2small[allskel[i]] = i;
  }
  for (int i = 0; i < allredundant.size(); i++) {
    red_big2small[allredundant[i]] = i;
  }

  *mu = b;

  ie_Mat modified_Psi = quadtree.Psi.transpose();
  ie_Mat modified_U = quadtree.U;

  // First apply the sweep matrices to x and U to modify them.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];

    #pragma omp parallel for num_threads(num_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
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

  // After the result of the first sweep matrices, grab w and z.
  ie_Mat w = (*mu)(allredundant, 0, 1);
  ie_Mat z = (*mu)(allskel, 0, 1);

  // Now apply the other sweep matrices to Psi to modify it.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(num_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
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
  // Again, C is mostly 0s, so we just apply Dinv to the nonzero block

  #pragma omp parallel for num_threads(num_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];
    if (current_node->src_dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    std::vector<int> small_redundants = big_to_small(
                                          current_node->src_dof_lists.redundant, red_big2small);
    assert(current_node->X_rr_is_LU_factored);
    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv,
                          &Dinv_w, small_redundants);
  }


  // B has modified_Psi as a subblock, but is otherwise a ton of 0s, so
  // we don't actually form B.

  // B_Dinv_w.set_submatrix(allskel.size(), allskel.size() + Psi.height(), 0, 1,
  //                       B_Dinv_w_nonzero);
  // The rest of the matrix is automatically zeros, as should be the case.

  ie_Mat M(allskel.size() + quadtree.Psi.height(), 1);

  M.set_submatrix(0, allskel.size(), 0, 1, z);
  M.set_submatrix(allskel.size(), M.height(), 0, 1, -(modified_Psi(0,
                  modified_Psi.height(),
                  allredundant) * Dinv_w));

  ie_Mat y(quadtree.S_LU.height(), 1);
  quadtree.S_LU.left_multiply_inverse(M, quadtree.S_piv, &y);

  *alpha =  y(allskel.size(), y.height(), 0, 1);

  ie_Mat N = w - modified_U(allredundant, 0, modified_U.width()) * (*alpha);
  ie_Mat Dinv_N = N;

  #pragma omp parallel for num_threads(num_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];
    if (current_node->src_dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    std::vector<int> small_redundants = big_to_small(
                                          current_node->src_dof_lists.redundant, red_big2small);
    assert(current_node->X_rr_is_LU_factored);

    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv,
                          &Dinv_N, small_redundants);
  }


  mu->set_submatrix(allredundant, 0, 1, Dinv_N);
  mu->set_submatrix(allskel, 0, 1, y(0, allskel.size(), 0, 1));

  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(num_threads)
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      // for(int n = 0; n < current_level->nodes.size(); n++){

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
