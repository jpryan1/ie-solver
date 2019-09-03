// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_QUADTREE_QUADTREE_H_
#define IE_SOLVER_QUADTREE_QUADTREE_H_

#include <cassert>
#include <vector>
#include "ie_solver/ie_mat.h"
#include "ie_solver/boundaries/boundary.h"

#define MAX_LEAF_DOFS 32

// TODO(John) all headers should have variable names explicit.

namespace ie_solver {

struct InteractionLists {
  std::vector<unsigned int> original_box, active_box, redundant, skel, near,
      skelnear, permutation;

  void set_rs_ranges(const std::vector<unsigned int>& prm, unsigned int sk,
                     unsigned int rd) {
    assert(prm.size() == sk + rd);

    for (unsigned int i = 0; i < sk; i++) {
      skel.push_back(active_box[prm[i]]);
      permutation.push_back(prm[i]);
    }
    for (unsigned int i = sk; i < sk + rd; i++) {
      redundant.push_back(active_box[prm[i]]);
      permutation.push_back(prm[i]);
    }
  }

  void set_skelnear_range(bool strong_admissibility) {
    for (unsigned int i = 0; i < skel.size(); i++) {
      skelnear.push_back(skel[i]);
    }
    if (strong_admissibility) {
      for (unsigned int i = 0; i < near.size(); i++) {
        skelnear.push_back(near[i]);
      }
    }
  }
};

struct QuadTreeNode {
  static unsigned int id_count;

  unsigned int id, level;
  bool is_leaf, compressed = false;
  double side_length;

  QuadTreeNode *tl, *tr, *bl, *br;
  QuadTreeNode* parent;
  QuadTreeNode* children[4];
  std::vector<QuadTreeNode*> neighbors;

  InteractionLists src_dof_lists, tgt_dof_lists;
  // For inverse operator
  ie_Mat T, L, U, X_rr, schur_update;
  // For forward operator
  ie_Mat src_T, tgt_T, X_rs, X_sr;


  // format is {BL, TL, TR, BR}
  double corners[8];

  QuadTreeNode() {
    id = id_count++;
    is_leaf = true;
    tl = NULL;
    tr = NULL;
    bl = NULL;
    br = NULL;
    for (int i = 0; i < 4; i++) children[i] = NULL;
  }
  // ~QuadTreeNode() {
  //   for (QuadTreeNode* child : children) {
  //     if (child) {
  //       delete child;
  //     }
  //   }
  // }
};

struct QuadTreeLevel {
  std::vector<QuadTreeNode*> nodes;
};

class QuadTree {
 public:
  int solution_dimension, domain_dimension;
  unsigned int no_proxy_level = 0;
  double min, max;
  Boundary* boundary;
  std::vector<double> domain_points;
  QuadTreeNode* root;
  std::vector<QuadTreeLevel*> levels;

  //////////////////////
  // TREE CONSTRUCTION
  //////////////////////

  ~QuadTree();
  void initialize_tree(Boundary* boundary,
                       const std::vector<double>& domain_points,
                       int solution_dimension_,
                       int domain_dimension_);
  // TODO(John) different scheme than bool is_boundary
  void recursive_add(QuadTreeNode* node, double x, double y,
                     unsigned int mat_ind, bool is_boundary);
  void get_descendent_neighbors(QuadTreeNode* big, QuadTreeNode* small);
  void node_subdivide(QuadTreeNode* node);
  void reset();
  void reset(Boundary* boundary_);

  ////////////////
  // UPDATING
  ////////////////

  void mark_neighbors_and_parents(QuadTreeNode* node);
  void perturb(const Boundary& new_boundary);

  void remove_inactive_dofs_at_level(int level);
  void populate_all_active_boxes();
  void populate_active_box(QuadTreeNode* node);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_QUADTREE_QUADTREE_H_
