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
  std::vector<int> original_box,
      active_box,
      redundant,
      skel,
      near,
      skelnear,
      permutation;

  void set_rs_ranges(const std::vector<int>& prm, int sk,
                     int rd) {
    assert(prm.size() == sk + rd);

    for (int i = 0; i < sk; i++) {
      skel.push_back(active_box[prm[i]]);
      permutation.push_back(prm[i]);
    }
    for (int i = sk; i < sk + rd; i++) {
      redundant.push_back(active_box[prm[i]]);
      permutation.push_back(prm[i]);
    }
  }

  void set_skelnear_range(bool strong_admissibility) {
    for (int i = 0; i < skel.size(); i++) {
      skelnear.push_back(skel[i]);
    }
    if (strong_admissibility) {
      for (int i = 0; i < near.size(); i++) {
        skelnear.push_back(near[i]);
      }
    }
  }
};

struct QuadTreeNode {
  static int id_count;
  int id, level, dofs_below;

  bool is_leaf, X_rr_is_LU_factored = false, compressed = false;
  double side_length, compression_ratio = 0., compress_time = 0.;

  QuadTreeNode *tl, *tr, *bl, *br;
  QuadTreeNode* parent;
  QuadTreeNode* children[4];
  std::vector<QuadTreeNode*> neighbors;

  InteractionLists src_dof_lists, tgt_dof_lists;
  // For inverse operator
  ie_Mat T, L, U, X_rr, schur_update, X_rr_lu;
  std::vector<lapack_int> X_rr_piv;
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
};


struct QuadTreeLevel {
  std::vector<QuadTreeNode*> nodes;
  ~QuadTreeLevel() {
    for (QuadTreeNode* node : nodes) {
      delete node;
    }
  }
};


class QuadTree {
 public:
  int solution_dimension, domain_dimension;
  int no_proxy_level = 0;
  double min, max;
  Boundary* boundary = nullptr;
  std::vector<double> domain_points;
  QuadTreeNode* root;

  ie_Mat allskel_mat, allskel_mat_lu, U, Psi, S_LU;
  std::vector<lapack_int> allskel_mat_piv, S_piv;

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
                     int mat_ind, bool is_boundary);
  void get_descendent_neighbors(QuadTreeNode* big, QuadTreeNode* small);
  void node_subdivide(QuadTreeNode* node);
  void consolidate_node(QuadTreeNode* node);
  void reset();
  void reset(Boundary* boundary_);

  void copy_into(QuadTree* new_tree) const;

  ////////////////
  // UPDATING
  ////////////////

  void mark_neighbors_and_parents(QuadTreeNode* node);
  void perturb(const Boundary& new_boundary);

  void remove_inactive_dofs_at_level(int level);
  void remove_inactive_dofs_at_all_boxes();
  void remove_inactive_dofs_at_box(QuadTreeNode* node);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_QUADTREE_QUADTREE_H_
