// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "ie-solver/tools/ie_solver_tools.h"
#include "ie-solver/kernel.h"

#define NODE_CAP INFINITY
#define LEVEL_CAP INFINITY

namespace ie_solver {

int IeSolverTools::interpolative_decomposition(const Kernel& kernel,
    const QuadTree* tree, QuadTreeNode* node) {
  assert(node != nullptr && "InterpolativeDecomposition fails on null node.");
  assert(node->interaction_lists.active_box.size() > 0 &&
         "Num of DOFs must be positive in InterpolativeDecomposition.");
  // TODO(John) better variable name
  ie_Mat pxy;

  make_id_mat(kernel, &pxy, tree, node);

  std::vector<unsigned int> p;
  unsigned int numskel = pxy.id(&p, &node->T, id_tol);

  if (numskel == 0) return 0;
  set_rs_ranges(&node->interaction_lists, p, node->T.height(),
                node->T.width());

  set_skelnear_range(&node->interaction_lists);

  return node->T.width();
}


void IeSolverTools::get_x_matrices(ie_Mat* K, const ie_Mat& Z, ie_Mat* Xrr,
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

  // this is just for readability
  *Xrr = (*K)(r, r);
  Xrs = (*K)(r, s);
  Xsr = (*K)(s, r);
  if (n_size > 0) {
    Xrn = (*K)(r, n);
    Xnr = (*K)(n, r);
  }

  ie_Mat::gemm(TRANSPOSE, NORMAL, -1., Z, (*K)(s, r), 1., Xrr);
  ie_Mat::gemm(TRANSPOSE, NORMAL, -1., Z, (*K)(s, s), 1., &Xrs);
  ie_Mat::gemm(NORMAL,    NORMAL, -1., Xrs,     Z,       1., Xrr);
  ie_Mat::gemm(NORMAL,    NORMAL, -1., (*K)(s, s), Z,       1., &Xsr);
  if (n_size > 0) {
    ie_Mat::gemm(TRANSPOSE, NORMAL, -1., Z, (*K)(s, n), 1., &Xrn);
    ie_Mat::gemm(NORMAL,   NORMAL, -1., (*K)(n, s), Z,     1., &Xnr);
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
void IeSolverTools::schur_update(const Kernel& kernel, QuadTreeNode* node) {
  assert(node != nullptr && "SchurUpdate fails on null node.");
  assert(node->T.height()*node->T.width() > 0 &&
         "Z must have positive dimensions in SchurUpdate.");
  assert(node->interaction_lists.active_box.size() > 0 &&
         "Num of DOFs must be positive in SchurUpdate.");
  // height of Z is number of skeleton columns
  unsigned int num_redundant = node->T.width();
  unsigned int num_skel     = node->T.height();
  // GENERATE K_BN,BN
  std::vector<unsigned int> BN;
  for (unsigned int idx : node->interaction_lists.active_box) {
    BN.push_back(idx);
  }
  if (strong_admissibility) {
    for (unsigned int idx : node->interaction_lists.near) BN.push_back(idx);
  }
  // Note that BN has all currently deactivated DoFs removed.
  ie_Mat K_BN = kernel(BN, BN);  // que bien!

  ie_Mat update(BN.size(), BN.size());
  get_all_schur_updates(&update, BN, node, strong_admissibility);
 
  K_BN -= update;
  // Generate various index ranges within BN
  std::vector<unsigned int> s, r, n, sn;
  for (unsigned int i = 0; i < num_skel; i++) {
    s.push_back(node->interaction_lists.permutation[i]);
    sn.push_back(node->interaction_lists.permutation[i]);
  }
  for (unsigned int i = 0; i < num_redundant; i++) {
    r.push_back(node->interaction_lists.permutation[i + num_skel]);
  }

  if (strong_admissibility) {
    for (unsigned int i = 0; i < node->interaction_lists.near.size(); i++) {
      n.push_back(i + num_redundant + num_skel);
      sn.push_back(i + num_redundant + num_skel);
    }
  }
  ie_Mat  Xrr(num_redundant, num_redundant);
  get_x_matrices(&K_BN, node->T, &Xrr, r, s, n);
  node->D_r = Xrr;

  // Generate left and right schur complement matrices
  // TODO(John) change this variable naming
  int num_skelnear = sn.size();

  node->L = ie_Mat(num_skelnear, num_redundant);
  node->U = ie_Mat(num_redundant, num_skelnear);

  double cond = Xrr.condition_number();
  if (cond > 1000) {
    std::cout << "Node " << node->id << " schur update -- ";
    std::cout << "Width " << node->side_length << " num dofs " << Xrr.width() <<
              std::endl;
    std::cout << "Inverting w/ condition number " << cond << std::endl;
  }

  Xrr.right_multiply_inverse(K_BN(sn, r), &node->L);
  Xrr.left_multiply_inverse(K_BN(r, sn), &node->U);

  ie_Mat schur(sn.size(), sn.size());
  ie_Mat::gemm(NORMAL, NORMAL, 1.0, node->L, K_BN(r, sn), 0., &schur);
  // set schur update
  node->schur_update = schur;
  node->schur_updated = true;
}


void IeSolverTools::skeletonize(const Kernel& kernel, QuadTree* tree) {
  int node_counter = 0;
  unsigned int lvls = tree->levels.size();
  int active_dofs = tree->boundary->points.size() / 2;
  for (unsigned int level = lvls - 1; level > 0; level--) {
    if (lvls - level > LEVEL_CAP) {
      break;
    }
    QuadTreeLevel* current_level = tree->levels[level];
    // First, get all active dofs from children
    for (QuadTreeNode * node : current_level->nodes) {
      if (node->schur_updated) continue;
      populate_active_box(node);
    }
    // Next, get all active near dofs from neighbors
    for (QuadTreeNode* node_a : current_level->nodes) {
      if (node_a->schur_updated) continue;
      node_a->interaction_lists.near.clear();
      for (QuadTreeNode* neighbor : node_a->neighbors) {
        // Some neighbors are smaller boxes from higher levels, we don't
        // care about those, their parents have the updated information.
        if (neighbor->level > node_a->level) {
          continue;
        }
        for (unsigned int idx : neighbor->interaction_lists.active_box) {
          node_a->interaction_lists.near.push_back(idx);
        }
      }
    }

    for (unsigned int n = 0; n < current_level->nodes.size(); n++) {
      node_counter++;
      if (node_counter > NODE_CAP) {
        break;
      }
      QuadTreeNode* current_node = current_level->nodes[n];
      if (current_node->schur_updated) {
        continue;
      }
      if (current_node->interaction_lists.active_box.size()
          < MIN_DOFS_TO_COMPRESS) {
        continue;
      }

      int redundants = interpolative_decomposition(kernel, tree, current_node);

      active_dofs -= redundants;
      if (redundants == 0) {
        continue;
      }
      schur_update(kernel, current_node);
    }
  }
  // If the above breaks due to a cap, we need to manually propagate active
  // boxes up the tree.
  // std::cout << "Final count " << active_dofs << std::endl;
  populate_all_active_boxes(tree);
  // check_factorization_against_kernel(kernel, tree);
}


}  // namespace ie_solver
