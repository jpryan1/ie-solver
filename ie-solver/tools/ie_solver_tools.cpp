// Copyright 2019 John Paul Ryan
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include "ie-solver/tools/ie_solver_tools.h"

namespace ie_solver {

// TODO(John) URGENTLY skelnear needs to be replaced as a variable name

// TODO(John) here and elsewhere check on passing by reference

IeSolverTools::IeSolverTools(double id_tol, bool strong_admissibility_,
                             int solution_dimension_) {
  assert(id_tol > 0 && "id_tol must be greater than one to init tools.");
  this->id_tol = id_tol;
  strong_admissibility = strong_admissibility_;
  solution_dimension = solution_dimension_;
}


void IeSolverTools::check_factorization_against_kernel(const Kernel& kernel,
    QuadTree* tree) {
  int check_size = 100;
  // This ensures that operations know what the remaining skels are.
  populate_all_active_boxes(tree);

  int dofs = tree->boundary->points.size() / 2;
  // Take a random 100x100 submatrix of A-A^hat, estimate 2-norm by power method
  std::vector<unsigned int> rand_x_indices, rand_y_indices;
  for (int i = 0; i < check_size; i++) {
    rand_x_indices.push_back(rand() % dofs);
    rand_y_indices.push_back(rand() % dofs);
  }

  ie_Mat A = kernel(rand_x_indices, rand_y_indices);
  ie_Mat A_hat(check_size, check_size);
  ie_Mat basis(dofs, 1);
  ie_Mat b(dofs, 1);

  for (int y = 0; y < rand_y_indices.size(); y++) {
    int rand_y_idx = rand_y_indices[y];
    for (unsigned int i = 0; i < dofs; i++) {
      basis.set(i, 0, 0);
    }
    basis.set(rand_y_idx, 0, 1.0);
    sparse_matvec(kernel, *tree, basis, &b);
    for (int x = 0; x < rand_x_indices.size(); x++) {
      int rand_x_idx = rand_x_indices[x];
      A_hat.set(x, y, b.get(rand_x_idx, 0));
    }
  }
  double truenorm = A.frob_norm();
  A -= A_hat;
  std::cout << "Total error: " << 100.0 * A.frob_norm() / truenorm << "%" <<
            std::endl;
}


void IeSolverTools::populate_all_active_boxes(QuadTree* tree) {
  int lvls = tree->levels.size();
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->schur_updated) {
        continue;
      }
      populate_active_box(current_node);
    }
  }
}


void IeSolverTools::populate_active_box(QuadTreeNode* node) {
  // this function removes from the box any DoFs which have already been made
  // redundant. It involves a bunch of annoying C++ functions and probably
  // would look nicer in matlab.

  // populate active_box
  node->interaction_lists.skel.clear();
  node->interaction_lists.skelnear.clear();
  node->interaction_lists.redundant.clear();
  node->interaction_lists.active_box.clear();
  if (!node->is_leaf) {
    for (QuadTreeNode* child : node->children) {
      if (child->schur_updated) {
        for (unsigned int i : child->interaction_lists.skel) {
          node->interaction_lists.active_box.push_back(i);
        }
      } else {
        for (unsigned int i : child->interaction_lists.active_box) {
          node->interaction_lists.active_box.push_back(i);
        }
      }
    }
  } else {
    node->interaction_lists.active_box = node->interaction_lists.original_box;
  }
}

// TODO(John) the permutation vector is weirdly unique among the
// interaction_lists and is weird to deal with in perturbing, can we just ditch
// it somehow?
void IeSolverTools::set_rs_ranges(InteractionLists* interaction_lists,
                                  const std::vector<unsigned int>& prm,
                                  unsigned int sk, unsigned int rd) {
  assert(prm.size() == sk + rd);

  for (unsigned int i = 0; i < sk; i++) {
    interaction_lists->skel.push_back(interaction_lists->active_box[prm[i]]);
    interaction_lists->permutation.push_back(prm[i]);
  }
  for (unsigned int i = sk; i < sk + rd; i++) {
    interaction_lists->redundant.push_back(
      interaction_lists->active_box[prm[i]]);
    interaction_lists->permutation.push_back(prm[i]);
  }
}

void IeSolverTools::set_skelnear_range(InteractionLists* interaction_lists) {
  for (unsigned int i = 0; i < interaction_lists->skel.size(); i++) {
    interaction_lists->skelnear.push_back(interaction_lists->skel[i]);
  }
  if (strong_admissibility) {
    for (unsigned int i = 0; i < interaction_lists->near.size(); i++) {
      interaction_lists->skelnear.push_back(interaction_lists->near[i]);
    }
  }
}

}  // namespace ie_solver
