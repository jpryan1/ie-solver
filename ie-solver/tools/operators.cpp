// Copyright 2019 John Paul Ryan
#include <iostream>
#include <cassert>
#include "ie-solver/tools/ie_solver_tools.h"
namespace ie_solver {

// Sets vec(b) = vec(b) + mat*vec(a)
void IeSolverTools::apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
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
void IeSolverTools::apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
                                      const std::vector<unsigned int>& range)
const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(), vec->width());
  ie_Mat::gemm(NORMAL, NORMAL, 1., mat, (*vec)(range, 0, vec->width()), 0.,
               &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void IeSolverTools::apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat* vec,
    const std::vector<unsigned int>& range) const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(),  vec->width());
  mat.left_multiply_inverse((*vec)(range,  0, vec->width()), &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void IeSolverTools::sparse_matvec(const Kernel& K, const QuadTree& tree,
                                  const ie_Mat& x,
                                  ie_Mat* b) const {
  *b = x;
  int lvls = tree.levels.size();
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->schur_updated) {
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
    QuadTreeLevel* current_level = tree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->schur_updated) {
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
    tree.root->src_dof_lists.active_box;

  if (allskel.size() > 0) {
    ie_Mat allskel_mat(allskel.size(), allskel.size());

    get_all_schur_updates(&allskel_mat, allskel, tree.root, false);

    allskel_mat = K(allskel, allskel) - allskel_mat;
    apply_diag_matrix(allskel_mat, b, allskel);
  }

  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = tree.levels[level];
    // TODO(John) record in notes and explore the following observation:
    // changing the order of this for loop affects the accuracy of the
    // sparse_mat_vec on a random vector, BUT NOT THE SOLUTION ERROR
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->schur_updated) {
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


void IeSolverTools::solve(const Kernel& K, const QuadTree& tree, ie_Mat* x,
                          const ie_Mat& b) const {
  assert(x->height() == b.height());
  int lvls = tree.levels.size();
  *x = b;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->schur_updated) {
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

  for (int level = lvls - 1; level >= 0; level--) {  // level>=0; level--){
    QuadTreeLevel* current_level = tree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->schur_updated) {
        continue;
      }
      double cond = current_node->X_rr.condition_number();
      if (cond > 1000) {
        std::cout << "Node " << current_node->id << " solve X_rr inv -- ";
        std::cout << "Inverting w/ condition number " << cond << std::endl;
      }
      apply_diag_inv_matrix(current_node->X_rr, x,
                            current_node->src_dof_lists.redundant);
    }
  }

  // We need all of the skeleton indices. This is just the negation of
  // [0,b.size()] and the redundant DoFs
  // with that in mind...
  std::vector<unsigned int> allskel = tree.root->src_dof_lists.active_box;
  if (allskel.size() > 0) {
    ie_Mat allskel_mat(allskel.size(), allskel.size());
    get_all_schur_updates(&allskel_mat, allskel, tree.root, false);
    allskel_mat = K(allskel, allskel) - allskel_mat;
    double cond2 = allskel_mat.condition_number();
    if (cond2 > 1000) {
      std::cout << "Allskel inv -- ";
      std::cout << "Inverting w/ condition number " << cond2 << std::endl;
    }
    apply_diag_inv_matrix(allskel_mat, x, allskel);
  }
  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = tree.levels[level];
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      // for(unsigned int n = 0; n < current_level->nodes.size(); n++){

      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->schur_updated) {
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

}  // namespace ie_solver
