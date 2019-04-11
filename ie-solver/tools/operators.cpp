// Copyright 2019 John Paul Ryan
#include <iostream>
#include <cassert>
#include "ie-solver/tools/ie_solver_tools.h"
namespace ie_solver {

// Sets vec(b) = vec(b) + mat*vec(a)
void IeSolverTools::apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
                                       const std::vector<unsigned int>& a,
                                       const std::vector<unsigned int>& b,
                                       bool transpose = false) {
  if (a.size()*b.size() == 0) return;
  if (transpose) {
    assert(mat.height() == a.size());
  } else {
    assert(mat.width() == a.size());
  }
  // This vector is just used for indexing an Vector
  std::vector<unsigned int> ZERO_VECTOR;
  ZERO_VECTOR.push_back(0);
  ie_Mat temp = (*vec)(a, ZERO_VECTOR);
  ie_Mat product(b.size(), 1);

  if (transpose) {
    ie_Mat::gemv(TRANSPOSE, 1., mat, temp, 0., &product);
  } else {
    ie_Mat::gemv(NORMAL, 1., mat, temp, 0., &product);
  }
  product += (*vec)(b, ZERO_VECTOR);
  vec->set_submatrix(b, ZERO_VECTOR, product);
}


// Sets vec(range) = mat * vec(range)
void IeSolverTools::apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
                                      const std::vector<unsigned int>& range) {
  if (range.size() == 0) return;
  std::vector<unsigned int> ZERO_VECTOR;
  ZERO_VECTOR.push_back(0);
  ie_Mat temp = (*vec)(range, ZERO_VECTOR);
  ie_Mat product(range.size(), 1);
  ie_Mat::gemv(NORMAL, 1., mat, temp, 0., &product);
  vec->set_submatrix(range, ZERO_VECTOR, product);
}


void IeSolverTools::b2d_apply_diag_matrix(const ie_Mat& mat,
    const std::vector<unsigned int>& tgt,
    const std::vector<unsigned int>& src,
    const ie_Mat& vec_in, ie_Mat* vec_out) {

  if (tgt.size()*src.size() == 0) return;
  std::vector<unsigned int> ZERO_VECTOR;
  ZERO_VECTOR.push_back(0);
  ie_Mat temp = vec_in(src, ZERO_VECTOR);

  ie_Mat product(tgt.size(), 1);
  ie_Mat::gemv(NORMAL, 1., mat, temp, 0., &product);
  product += (*vec_out)(tgt, ZERO_VECTOR);
  vec_out->set_submatrix(tgt, ZERO_VECTOR, product);
}


void IeSolverTools::apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat* vec,
    const std::vector<unsigned int>& range) {
  if (range.size() == 0) return;

  std::vector<unsigned int> ZERO_VECTOR;
  ZERO_VECTOR.push_back(0);
  ie_Mat temp = (*vec)(range, ZERO_VECTOR);
  ie_Mat product(range.size(), 1);

  mat.left_multiply_inverse(temp, &product);

  vec->set_submatrix(range, ZERO_VECTOR, product);
}


void IeSolverTools::b2dsparse_matvec(const Kernel& K, const QuadTree& tree,
                                     const ie_Mat& x,
                                     ie_Mat* b) {
  int lvls = tree.levels.size();

  for (int level_ = lvls - 2; level_ > 0; level_--) {
    ie_Mat add(b->height(), 1);
    ie_Mat start = x;

    // T _src inv
    for (int level = lvls - 1; level > level_; level--) {
      QuadTreeLevel* current_level = tree.levels[level];
      for (QuadTreeNode* current_node : current_level->nodes) {
        if (!current_node->schur_updated) continue;
        apply_sweep_matrix(current_node->src_T, &start,
                           current_node->src_dof_lists.redundant,
                           current_node->src_dof_lists.skel,
                           false);
      }
    }

    // K_r
    for (int level = lvls - 1; level > level_; level--) {
      QuadTreeLevel* current_level = tree.levels[level];
      for (QuadTreeNode* current_node : current_level->nodes) {
        if (!current_node->schur_updated) {
          continue;
        }
        if (current_node->X_rs.height() * current_node->X_rs.width() > 0) {
          b2d_apply_diag_matrix(current_node->X_rs,
                                current_node->tgt_dof_lists.redundant,
                                current_node->src_dof_lists.skel,
                                start, &add);
        }
        if (current_node->X_rr.height() * current_node->X_rr.width() > 0) {
          b2d_apply_diag_matrix(current_node->X_rr,
                                current_node->tgt_dof_lists.redundant,
                                current_node->src_dof_lists.redundant,
                                start, &add);
        }
        if (current_node->X_sr.height() * current_node->X_sr.width() > 0) {
          b2d_apply_diag_matrix(current_node->X_sr,
                                current_node->tgt_dof_lists.skel,
                                current_node->src_dof_lists.redundant,
                                start, &add);
        }
      }
    }
    std::cout << "Norm: " << add.frob_norm() << std::endl;

    // T_tgt inv
    for (int level = level_ + 1; level < lvls; level++) {
      QuadTreeLevel* current_level = tree.levels[level];
      for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
        QuadTreeNode* current_node = current_level->nodes[n];
        if (!current_node->schur_updated) continue;
        apply_sweep_matrix(current_node->tgt_T, &add,
                           current_node->tgt_dof_lists.skel,
                           current_node->tgt_dof_lists.redundant, true);
      }
    }
    (*b) += add;
  }


  ie_Mat final_start = x;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree.levels[level];
    for (QuadTreeNode* current_node : current_level->nodes) {
      if (!current_node->schur_updated) {
        continue;
      }
      apply_sweep_matrix(current_node->src_T, &final_start,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skel,
                         false);
    }
  }
// We need all of the skeleton indices. This is just the negation of
// [0,b.size()] and the redundant DoFs
// with that in mind...





  ie_Mat final_add(b->height(), 1);

  QuadTreeLevel* current_level = tree.levels[0];
  for (QuadTreeNode* current_node : current_level->nodes) {
    if (!current_node->schur_updated) {
      continue;
    }
    if (current_node->X_rs.height() * current_node->X_rs.width() > 0) {
      b2d_apply_diag_matrix(current_node->X_rs,
                            current_node->tgt_dof_lists.redundant,
                            current_node->src_dof_lists.skel,
                            final_start, &final_add);
    }
    if (current_node->X_rr.height() * current_node->X_rr.width() > 0) {
      b2d_apply_diag_matrix(current_node->X_rr,
                            current_node->tgt_dof_lists.redundant,
                            current_node->src_dof_lists.redundant,
                            final_start, &final_add);
    }
    if (current_node->X_sr.height() * current_node->X_sr.width() > 0) {
      b2d_apply_diag_matrix(current_node->X_sr,
                            current_node->tgt_dof_lists.skel,
                            current_node->src_dof_lists.redundant,
                            final_start, &final_add);
    }
  }


  std::vector<unsigned int> src_allskel =
    tree.root->src_dof_lists.active_box;

  std::vector<unsigned int> tgt_allskel =
    tree.root->tgt_dof_lists.active_box;

  if (src_allskel.size()*tgt_allskel.size() > 0) {
    ie_Mat allskel_mat = K.forward_get(tgt_allskel, src_allskel);
    b2d_apply_diag_matrix(allskel_mat,
                          tgt_allskel,
                          src_allskel,
                          final_start, &final_add);
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

      // Finally we need to apply U_T inverse
      // U_T inverse changes the redundant elements - it makes them equal
      // to T transpose times the skeleton elements + the redundant
      // elements
      apply_sweep_matrix(current_node->tgt_T, &final_add,
                         current_node->tgt_dof_lists.skel,
                         current_node->tgt_dof_lists.redundant, true);
    }
  }

  (*b) += final_add;

  for (int i = 0; i < tree.domain_points.size(); i += 2) {
    if (!tree.boundary->is_in_domain(
          Vec2(tree.domain_points[i], tree.domain_points[i + 1]))) {
      for (int j = 0; j < solution_dimension; j++) {
        b->set((i / 2) * solution_dimension + j, 0, 0.0);
      }
    }
  }





//***************************************************************************




  // ie_Mat s2s = x;
  // int lvls = tree.levels.size();

  // for (int level = lvls - 1; level >= 0; level--) {
  //   QuadTreeLevel* current_level = tree.levels[level];
  //   for (QuadTreeNode* current_node : current_level->nodes) {
  //     if (!current_node->schur_updated) {
  //       continue;
  //     }
  //     // First we need to apply L_T inverse
  //     // L_T inverse changes the skel elements - it makes them equal to
  //     // T times the redundant elements + the skeleton elements.
  //     apply_sweep_matrix(current_node->src_T, &s2s,
  //                        current_node->src_dof_lists.redundant,
  //                        current_node->src_dof_lists.skel,
  //                        false);
  //   }
  // }

  // ie_Mat s2t(b->height(), 1);
  // // ^ all zeros for now.

  // for (int level = lvls - 1; level >= 0; level--) {
  //   QuadTreeLevel* current_level = tree.levels[level];
  //   for (QuadTreeNode* current_node : current_level->nodes) {
  //     if (!current_node->schur_updated) {
  //       continue;
  //     }
  //     if (current_node->X_rs.height() * current_node->X_rs.width() > 0) {
  //       b2d_apply_diag_matrix(current_node->X_rs,
  //                             current_node->tgt_dof_lists.redundant,
  //                             current_node->src_dof_lists.skel,
  //                             s2s, b);
  //     }
  //     if (current_node->X_rr.height() * current_node->X_rr.width() > 0) {
  //       b2d_apply_diag_matrix(current_node->X_rr,
  //                             current_node->tgt_dof_lists.redundant,
  //                             current_node->src_dof_lists.redundant,
  //                             s2s, b);
  //     }
  //     if (current_node->X_sr.height() * current_node->X_sr.width() > 0) {
  //       b2d_apply_diag_matrix(current_node->X_sr,
  //                             current_node->tgt_dof_lists.skel,
  //                             current_node->src_dof_lists.redundant,
  //                             s2s, b);
  //     }
  //   }
  // }
  // std::cout<<"Norm: "<<b->frob_norm()<<std::endl;

  // // We need all of the skeleton indices. This is just the negation of
  // // [0,b.size()] and the redundant DoFs
  // // with that in mind...

  // std::vector<unsigned int> src_allskel =
  //   tree.root->src_dof_lists.active_box;

  // std::vector<unsigned int> tgt_allskel =
  //   tree.root->tgt_dof_lists.active_box;

  // if (src_allskel.size()*tgt_allskel.size() > 0) {
  //   ie_Mat allskel_mat = K.forward_get(tgt_allskel, src_allskel);
  //   b2d_apply_diag_matrix(allskel_mat,
  //                         tgt_allskel,
  //                         src_allskel,
  //                         s2s, b);
  // }


  // for (int level = 0; level < lvls; level++) {
  //   QuadTreeLevel* current_level = tree.levels[level];
  //   // TODO(John) record in notes and explore the following observation:
  //   // changing the order of this for loop affects the accuracy of the
  //   // sparse_mat_vec on a random vector, BUT NOT THE SOLUTION ERROR
  //   for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
  //     QuadTreeNode* current_node = current_level->nodes[n];
  //     if (!current_node->schur_updated) {
  //       continue;
  //     }

  //     // Finally we need to apply U_T inverse
  //     // U_T inverse changes the redundant elements - it makes them equal
  //     // to T transpose times the skeleton elements + the redundant
  //     // elements
  //     apply_sweep_matrix(current_node->tgt_T, b,
  //                        current_node->tgt_dof_lists.skel,
  //                        current_node->tgt_dof_lists.redundant, true);
  //   }
  // }


  // for (int i = 0; i < tree.domain_points.size(); i += 2) {
  //   if (!tree.boundary->is_in_domain(
  //         Vec2(tree.domain_points[i], tree.domain_points[i + 1]))) {
  //     for (int j = 0; j < solution_dimension; j++) {
  //       b->set((i / 2) * solution_dimension + j, 0, 0.0);

  //     }
  //   }
  // }




}


void IeSolverTools::sparse_matvec(const Kernel& K, const QuadTree& tree,
                                  const ie_Mat& x,
                                  ie_Mat* b) {
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

    allskel_mat *= -1;
    allskel_mat += K(allskel, allskel);
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
                          const ie_Mat& b) {
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
      current_node->L *= -1;
      current_node->T *= -1;
      // Finally we need to apply U_T inverse
      // U_T inverse changes the redundant elements - it makes them equal
      // to T transpose times the skeleton elements + the redundant
      // elements
      apply_sweep_matrix(current_node->T, x,
                         current_node->src_dof_lists.skel,
                         current_node->src_dof_lists.redundant, true);
      apply_sweep_matrix(current_node->L, x,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skelnear, false);

      current_node->L *= -1;
      current_node->T *= -1;
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
      // if (cond > 1000) {
      //   std::cout << "Node " << current_node->id << " solve X_rr inv -- ";
      //   std::cout << "Inverting w/ condition number " << cond << std::endl;
      // }
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
    allskel_mat *= -1;
    allskel_mat += K(allskel, allskel);
    double cond2 = allskel_mat.condition_number();
    // if (cond2 > 1000) {
    //   std::cout << "Allskel inv -- ";
    //   std::cout << "Inverting w/ condition number " << cond2 << std::endl;
    // }
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
      current_node->U *= -1;
      current_node->T *= -1;
      // First we need to apply L_T inverse
      // L_T inverse changes the skel elements - it makes them equal to T
      // times the redundant elements + the skeleton elements.
      // Next we need to apply U inverse
      // U inverse changes the redundant elements - it makes them equal to
      // L transpose times the skelnear elements + the redundant elements
      apply_sweep_matrix(current_node->U, x,
                         current_node->src_dof_lists.skelnear,
                         current_node->src_dof_lists.redundant, false);
      apply_sweep_matrix(current_node->T, x,
                         current_node->src_dof_lists.redundant,
                         current_node->src_dof_lists.skel,
                         false);

      current_node->U *= -1;
      current_node->T *= -1;
    }
  }
}


// void IeSolverTools::Conjugate_Gradient(ie_Mat& K, QuadTree& tree,
// ie_Mat& phi,
//   ie_Mat& f){
//   int iteration = 0;
//   Zero(phi);
//   ie_Mat r(f.height(), f.width());
//   Copy(f, r);
//   ie_Mat p(f.height(), f.width());
//   Copy(f,p);

//   ie_Mat Ap(p.height(), 1);


//   double r_dot = Dot(r, r);

//   for(; iteration<100; iteration++){
//     //r^T r

//     sparse_mat_vec(F, tree, p, Ap);
//     double p_dot = Dot(p, Ap);

//     double alpha = r_dot / p_dot;

//     //fix this stupidity
//     p *= alpha;
//     phi += p;
//     p*= 1.0/alpha;

//     Ap *= alpha;
//     r -= Ap;
//     Ap *= 1.0/alpha;



//     double new_rdot = Dot(r, r);

//     if(new_rdot < 0.0000001){
//       printf("Done in %d iterations\n", iteration);
//       break;
//     }
//     double beta = new_rdot / r_dot;

//     r_dot = new_rdot;

//     p *= beta;
//     p += r;


//   }

// }

}  // namespace ie_solver
