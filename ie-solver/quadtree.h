// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_QUADTREE_H_
#define IE_SOLVER_QUADTREE_H_

#include <vector>
#include "ie-solver/ie_mat.h"
#include "boundaries/boundary.h"

#define MAX_LEAF_DOFS 32

// TODO(John) all headers should have variable names explicit.

namespace ie_solver {

struct InteractionLists {
  std::vector<unsigned int> original_box, active_box, redundant, skel, near,
      skelnear, permutation;
};

struct QuadTreeNode {
  static unsigned int id_count;

  unsigned int id, level;
  bool is_leaf, schur_updated = false;
  double side_length;

  QuadTreeNode *tl, *tr, *bl, *br;
  QuadTreeNode* parent;
  QuadTreeNode* children[4];
  std::vector<QuadTreeNode*> neighbors;

  InteractionLists interaction_lists;
  ie_Mat T, L, U, D_r, schur_update;

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
  int solution_dimension;
  double min, max;
  Boundary* boundary;
  QuadTreeNode* root;
  std::vector<QuadTreeLevel*> levels;

  ~QuadTree() {
    for (QuadTreeLevel* level : levels) {
      for (QuadTreeNode* node : level->nodes) {
        delete node;
      }
    }
    for (QuadTreeLevel* level : levels) {
      if (level) {
        delete level;
      }
    }
    levels.clear();
  }

  void reset();
  void reset(Boundary* boundary_);

  void initialize_tree(Boundary* boundary, int solution_dimension);
  void recursive_add(QuadTreeNode* node, double x, double y,
                     unsigned int mat_ind);
  void get_descendent_neighbors(QuadTreeNode* big, QuadTreeNode* small);

  void node_subdivide(QuadTreeNode* node);

  void write_quadtree_to_file();

  void mark_neighbors_and_parents(QuadTreeNode* node);

  void perturb(const Boundary& new_boundary);

  // void print();
  // void rec_print(QuadTreeNode*);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_QUADTREE_H_
