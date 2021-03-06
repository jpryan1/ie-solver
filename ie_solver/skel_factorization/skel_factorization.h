// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_SKEL_FACTORIZATION_SKEL_FACTORIZATION_H_
#define IE_SOLVER_SKEL_FACTORIZATION_SKEL_FACTORIZATION_H_

#include <atomic>
#include <vector>
#include "ie_solver/ie_mat.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/kernel/kernel.h"

#define MIN_DOFS_TO_COMPRESS 16
#define NODE_CAP INFINITY
#define LEVEL_CAP INFINITY

namespace ie_solver {

// Perhaps putting these in a class isn't best, maybe just a namespace.
// Or maybe make everything static.
class SkelFactorization {
 public:
  double id_tol, make_mat_time, id_time, schur_time;
  int num_threads;
  bool strong_admissibility;
  // // int solution_dimension, domain_dimension;
  // ie_Mat allskel_mat, allskel_mat_lu, U, Psi, S_LU;
  // std::vector<lapack_int> allskel_mat_piv, S_piv;
  // std::atomic<bool> kill_factorizer;
  // std::atomic<QuadTreeNode*> block_to_factorize;

  SkelFactorization() {}
  SkelFactorization(double id_tol, bool strong_admissibility, int nt);
  ~SkelFactorization() {}

  void decouple(const Kernel& K, QuadTreeNode* node);
  int id_compress(const Kernel& K, const QuadTree* tree,
                  QuadTreeNode* node);
  // int b2dinterpolative_decomposition(const Kernel& K, const QuadTree* tree,
  //                                   QuadTreeNode* node);

  void get_all_schur_updates(ie_Mat* updates,
                             const std::vector<int>& BN,
                             const QuadTreeNode* node,
                             bool get_neighbors) const;
  void get_descendents_updates(ie_Mat* updates,
                               const std::vector<int>& BN,
                               const QuadTreeNode* node) const;
  void get_update(ie_Mat* updates, const std::vector<int>& BN,
                  const QuadTreeNode* node) const;

  // void b2d_apply_diag_matrix(const ie_Mat& mat,
  //                           const std::vector<int>& tgt,
  //                           const std::vector<int>& src,
  //                           const ie_Mat& vec_in, ie_Mat* vec_out);

  void skeletonize(const Kernel& K, QuadTree* tree);

  // void diag_block_factorizer();
  // void b2dskeletonize(const Kernel& K, QuadTree* tree);

  // void b2dsparse_matvec(const Kernel& K, const QuadTree& tree,
  // const ie_Mat& x, ie_Mat* b)

  void make_id_mat(const Kernel& K, ie_Mat* pxy, const QuadTree* tree,
                   const QuadTreeNode* node);
  ie_Mat make_proxy_mat(const Kernel& kernel, double cntr_x,
                        double cntr_y, double r,
                        const QuadTree* tree,
                        const std::vector<int>& box_indices);

  // void make_src_id_mat(const Kernel& K, ie_Mat* pxy, const QuadTree* tree,
  //                      const QuadTreeNode* node);

  // void make_tgt_id_mat(const Kernel& K, ie_Mat* pxy, const QuadTree* tree,
  //                      const QuadTreeNode* node);

  ////////////////
  // SOLVING
  ////////////////

  void apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
                          const std::vector<int>& a,
                          const std::vector<int>& b,
                          bool transpose) const;
  void apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
                         const std::vector<int>& range) const;
  void apply_diag_inv_matrix(const ie_Mat& mat,
                             const std::vector<lapack_int>& piv, ie_Mat* vec,
                             const std::vector<int>& range) const;
  void sparse_matvec(const QuadTree& quadtree, const ie_Mat& x,
                     ie_Mat* b) const;

  void solve(const QuadTree& quadtree, ie_Mat* x, const ie_Mat& b) const;
  void multiply_connected_solve(const QuadTree& quadtree, ie_Mat* x,
                                ie_Mat* alpha, const ie_Mat& b) const;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_SKEL_FACTORIZATION_SKEL_FACTORIZATION_H_
