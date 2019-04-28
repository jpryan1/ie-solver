// Copyright 2019 John Paul Ryan
#include <cassert>
#include <iostream>
#include "ie-solver/tools/ie_solver_tools.h"

namespace ie_solver {

void IeSolverTools::get_all_schur_updates(ie_Mat* updates,
    const std::vector<unsigned int>& BN, const QuadTreeNode* node,
    bool get_neighbors) const {
  assert(node != nullptr && "get_all_schur_updates fails on null node.");
  assert(BN.size() > 0 && "get_all_schur_updates needs positive num of DOFs");
  if (!node->is_leaf) get_descendents_updates(updates, BN, node);

  if (get_neighbors) {
    for (QuadTreeNode* neighbor : node->neighbors) {
      if (neighbor->level != node->level) continue;
      if (neighbor->schur_updated) get_update(updates, BN, neighbor);
      if (!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);
    }
  }
}


void IeSolverTools::get_descendents_updates(ie_Mat* updates,
    const std::vector<unsigned int>& BN, const QuadTreeNode* node)  const {
  assert(node != nullptr && "get_descendents_updates fails on null node.");
  assert(!node->is_leaf &&
         "get_descendents_updates must be called on non-leaf.");

  // by assumption, node is not a leaf
  for (QuadTreeNode* child : node->children) {
    if (child->schur_updated) get_update(updates, BN, child);
    if (!child->is_leaf) get_descendents_updates(updates, BN, child);
  }
}


void IeSolverTools::get_update(ie_Mat* update,
                               const std::vector<unsigned int>& BN,
                               const QuadTreeNode* node)  const {
  // node needs to check all its dofs against BN, enter interactions into
  // corresponding locations
  // node only updated its own BN dofs, and the redundant ones are no longer
  // relevant, so we only care about child's SN dofs
  std::vector<unsigned int> sn = node->src_dof_lists.skelnear;
  // First create a list of Dofs that are also in node's skelnear,
  // and with each one give the index in skelnear and the index in BN
  std::vector<unsigned int> BN_;
  std::vector<unsigned int> sn_;
  for (unsigned int i = 0; i < sn.size(); i++) {
    for (unsigned int j = 0; j < BN.size(); j++) {
      if (BN[j] == sn[i]) {
        sn_  .push_back(i);
        BN_  .push_back(j);
      }
    }
  }
  for (unsigned int i = 0; i < BN_  .size(); i++) {
    for (unsigned int j = 0; j < BN_  .size(); j++) {
      update->addset(BN_[i], BN_[j], node->schur_update.get(sn_[i], sn_[j]));
    }
  }
}

}  // namespace ie_solver
