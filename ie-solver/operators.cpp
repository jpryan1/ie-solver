// Copyright 2019 John Paul Ryan
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include "ie-solver/quadtree.h"
namespace ie_solver {

// Sets vec(b) = vec(b) + mat*vec(a)
void QuadTree::apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
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
void QuadTree::apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
                                 const std::vector<unsigned int>& range)
const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(), vec->width());
  ie_Mat::gemm(NORMAL, NORMAL, 1., mat, (*vec)(range, 0, vec->width()), 0.,
               &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void QuadTree::apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat* vec,
                                     const std::vector<unsigned int>& range) const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(),  vec->width());
  mat.left_multiply_inverse((*vec)(range,  0, vec->width()), &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void QuadTree::apply_diag_pinv_matrix(const ie_Mat& mat, ie_Mat* vec,
                                      const std::vector<unsigned int>& range) const {
  if (range.size() == 0) return;
  ie_Mat product(range.size(),  vec->width());
  mat.left_multiply_pseudoinverse((*vec)(range,  0, vec->width()), &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void QuadTree::sparse_matvec(const ie_Mat& x, ie_Mat* b) const {
  *b = x;
  int lvls = levels.size();
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = levels[level];
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
    QuadTreeLevel* current_level = levels[level];
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
    root->src_dof_lists.active_box;

  if (allskel.size() > 0) {
    apply_diag_matrix(allskel_mat, b, allskel);
  }

  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = levels[level];
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


void QuadTree::solve(ie_Mat* x, const ie_Mat& b) const {
  assert(x->height() == b.height());
  int lvls = levels.size();
  *x = b;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = levels[level];
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
    QuadTreeLevel* current_level = levels[level];
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
  std::vector<unsigned int> allskel = root->src_dof_lists.active_box;
  allskel_mat.write_singular_values_to_file("allskel_sing_vals.txt");
  if (allskel.size() > 0) {
    double cond2 = allskel_mat.condition_number();
    if (cond2 > 1000) {
      std::cout << "Allskel inv -- ";
      std::cout << "Inverting w/ condition number " << cond2 << std::endl;
    }
    apply_diag_pinv_matrix(allskel_mat, x, allskel);
  }
  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = levels[level];
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


std::vector<unsigned int> big_to_small(const std::vector<unsigned int>& big,
                                       const std::unordered_map<unsigned int,
                                       unsigned int>& map) {
  std::vector<unsigned int> small;
  for (unsigned int idx : big) {
    small.push_back(map.at(idx));
  }
  return small;
}


void QuadTree::multiply_connected_solve(ie_Mat* x, ie_Mat* alpha,
                                        const ie_Mat& b) const {
  assert(x->height() == b.height());
  int lvls = levels.size();
  *x = b;

  ie_Mat modified_Psi = Psi.transpose();
  ie_Mat modified_U = U;


  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->schur_updated) {
        continue;
      }

      apply_sweep_matrix(-current_node->T, x,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, x,
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

// modified psi transpose = T^R transpose * Psi transpose

  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->schur_updated) {
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

  ////////////////////////////////////////////////////////////////
  // THE BELOW IS ALL SETTING THE STAGE
  std::vector<unsigned int> allskel = root->src_dof_lists.active_box;
  int total_redundant = boundary->weights.size() * solution_dimension -
                        allskel.size();

  ie_Mat A(allskel.size() + Psi.height(),
           allskel.size() + Psi.height());
  ie_Mat B(allskel.size() + Psi.height(), total_redundant);
  ie_Mat C(total_redundant, allskel.size() + Psi.height());
  ie_Mat ident(Psi.height(), Psi.height());
  ident.eye(Psi.height());

  A.set_submatrix(0, allskel.size(),
                  0, allskel.size(), allskel_mat);
  A.set_submatrix(allskel.size(), A.height(),
                  allskel.size(), A.width(),
                  -ident);
  A.set_submatrix(allskel.size(), A.height(),
                  0, allskel.size(),
                  modified_Psi(0, Psi.height(), allskel));
  A.set_submatrix(0, allskel.size(),
                  allskel.size(), A.width(),
                  modified_U(allskel, 0, U.width()));

  std::vector<unsigned int> allredundant;
  std::vector<unsigned int> sorted_allskel = allskel;
  std::sort(sorted_allskel.begin(), sorted_allskel.end());
  unsigned int skel_idx = 0;
  for (unsigned int i = 0; i < boundary->weights.size()*solution_dimension; i++) {
    if (i == sorted_allskel[skel_idx]) {
      skel_idx++;
    } else {
      allredundant.push_back(i);
    }
  }
  B.set_submatrix(allskel.size(), B.height(),
                  0, B.width(),
                  modified_Psi(0, modified_Psi.height(), allredundant));
  C.set_submatrix(0, C.height(),
                  allskel.size(), C.width(),
                  modified_U(allredundant, 0, modified_U.width()));
  std::unordered_map<unsigned int, unsigned int> skel_big2small, red_big2small;
  for (int i = 0; i < allskel.size(); i++) {
    skel_big2small[allskel[i]] = i;
  }
  for (int i = 0; i < allredundant.size(); i++) {
    red_big2small[allredundant[i]] = i;
  }

  // THE ABOVE IS ALL SETTING THE STAGE
  ////////////////////////////////////////////////////////////////

  ie_Mat b_vec = (*x)(allredundant, 0, 1);
  ie_Mat Dinv_b = b_vec;
  for (int level = lvls - 1; level >= 0; level--) {  // level>=0; level--){
    QuadTreeLevel* current_level = levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->schur_updated) {
        continue;
      }
      std::vector<unsigned int> small_redundants = big_to_small(
            current_node->src_dof_lists.redundant, red_big2small);

      apply_diag_inv_matrix(current_node->X_rr, &Dinv_b,
                            small_redundants);
    }
  }
  ie_Mat B_Dinv_b(B.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., B, Dinv_b, 0., &B_Dinv_b);


  ie_Mat a(B.height(), 1);

  a.set_submatrix(0, allskel.size(), 0, 1, (*x)(allskel, 0, 1));
  ie_Mat first_paren = a - B_Dinv_b;

  // Now we need to form S, the schur complement.
  ie_Mat Dinv_C = C;
  for (int level = lvls - 1; level >= 0; level--) {  // level>=0; level--){
    QuadTreeLevel* current_level = levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->schur_updated) {
        continue;
      }
      std::vector<unsigned int> small_redundants = big_to_small(
            current_node->src_dof_lists.redundant, red_big2small);

      apply_diag_inv_matrix(current_node->X_rr, &Dinv_C,
                            small_redundants);
    }
  }
  ie_Mat B_Dinv_C(B.height(), C.width());
  ie_Mat::gemm(NORMAL, NORMAL, 1., B, Dinv_C, 0., &B_Dinv_C);
  ie_Mat S = A - B_Dinv_C;

  ie_Mat x_vec(S.height(), 1);
  S.left_multiply_inverse(first_paren, &x_vec);

  *alpha = x_vec(allskel.size(), x_vec.height(), 0, 1);
  ie_Mat Cx(C.height(), 1);
  ie_Mat::gemm(NORMAL, NORMAL, 1., C, x_vec, 0., &Cx);
  ie_Mat second_paren = b_vec - Cx;
  ie_Mat y = second_paren;
  for (int level = lvls - 1; level >= 0; level--) {  // level>=0; level--){
    QuadTreeLevel* current_level = levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (current_node->src_dof_lists.redundant.size() == 0) continue;
      if (!current_node->schur_updated) {
        continue;
      }
      std::vector<unsigned int> small_redundants = big_to_small(
            current_node->src_dof_lists.redundant, red_big2small);

      apply_diag_inv_matrix(current_node->X_rr, &y,
                            small_redundants);
    }
  }
  x->set_submatrix(allredundant, 0, 1, y);
  x->set_submatrix(allskel, 0, 1, x_vec(0, allskel.size(), 0, 1));
  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = levels[level];
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
