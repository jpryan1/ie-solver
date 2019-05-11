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
  ie_Mat allskel_mat, U, Psi;

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
  // PRINTING
  ////////////////

  void write_quadtree_to_file();
  void print_node_norms();
  void rec_print_norms(QuadTreeNode* node);

  ////////////////
  // UPDATING
  ////////////////

  void mark_neighbors_and_parents(QuadTreeNode* node);
  void perturb(const Boundary& new_boundary);

  ////////////////
  // SOLVING
  ////////////////

  void apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
                          const std::vector<unsigned int>& a,
                          const std::vector<unsigned int>& b,
                          bool transpose) const;
  void apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
                         const std::vector<unsigned int>& range) const;
  void apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat* vec,
                             const std::vector<unsigned int>& range) const;
  void apply_diag_pinv_matrix(const ie_Mat& mat, ie_Mat* vec,
                              const std::vector<unsigned int>& range) const;
  void sparse_matvec(const ie_Mat& x, ie_Mat* b) const;

  void solve(ie_Mat* x, const ie_Mat& b) const;
  void multiply_connected_solve(ie_Mat* x, ie_Mat* alpha, const ie_Mat& b) const;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_QUADTREE_H_
