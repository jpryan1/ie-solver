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
                             bool is_stokes_) {
  assert(id_tol > 0 && "id_tol must be greater than one to init tools.");
  this->id_tol = id_tol;
  strong_admissibility = strong_admissibility_;
  is_stokes = is_stokes_;
}


void IeSolverTools::check_factorization_against_kernel(const Kernel& kernel,
    QuadTree* tree) {
  assert(tree->boundary->points.size() % 2 == 0);
  populate_all_active_boxes(tree);
  unsigned int dofs = tree->boundary->points.size() / 2;
  int rand_idx = rand() % dofs;
  ie_Mat e1(dofs, 1);
  for (unsigned int i = 0; i < dofs; i++) {
    e1.set(i, 0, 0);
  }
  e1.set(rand_idx, 0, 1.0);
  ie_Mat b(dofs, 1);
  sparse_matvec(kernel, *tree, e1, &b);
  double error = 0;
  for (int i = 0; i < dofs; i++) {
    error += pow(kernel.get(i, rand_idx) - b.get(i, 0), 2);
  }
  std::cout << "Total row error: " << sqrt(error) << std::endl;
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


void IeSolverTools::set_rs_ranges(InteractionLists* interaction_lists,
                                  const std::vector<unsigned int>& prm,
                                  unsigned int sk, unsigned int rd) {
  assert(prm.size() == sk + rd);

  for (unsigned int i = 0; i < sk; i++) {
    interaction_lists->skel.push_back(interaction_lists->active_box[prm[i]]);
    interaction_lists->permutation.push_back(prm[i]);
  }
  for (unsigned int i = sk; i < sk + rd; i++) {
    interaction_lists->redundant.push_back(interaction_lists->active_box[prm[i]]);
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
