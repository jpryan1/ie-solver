// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_TOOLS_IE_SOLVER_TOOLS_H_
#define IE_SOLVER_TOOLS_IE_SOLVER_TOOLS_H_

#include <vector>
#include "ie-solver/ie_mat.h"
#include "ie-solver/quadtree.h"
#include "ie-solver/kernel.h"

#define MIN_DOFS_TO_COMPRESS 16

namespace ie_solver {

// Perhaps putting these in a class isn't best, maybe just a namespace.
// Or maybe make everything static.
class IeSolverTools {
 public:
  double id_tol;
  bool strong_admissibility;
  int solution_dimension, domain_dimension;
  IeSolverTools() {}
  IeSolverTools(double id_tol, bool strong_admissibility,
                int solution_dimension, int domain_dimension);
  ~IeSolverTools() {}

  void get_x_matrices(ie_Mat* K, const ie_Mat& Z, ie_Mat* Xrr,
                      const std::vector<unsigned int>& r,
                      const std::vector<unsigned int>& s,
                      const std::vector<unsigned int>& n);

  void schur_update(const Kernel& K, QuadTreeNode* node);
  int interpolative_decomposition(const Kernel& K, const QuadTree* tree,
                                  QuadTreeNode* node);

  void get_all_schur_updates(ie_Mat* updates,
                             const std::vector<unsigned int>& BN,
                             const QuadTreeNode* node, bool get_neighbors);
  void get_descendents_updates(ie_Mat* updates,
                               const std::vector<unsigned int>& BN,
                               const QuadTreeNode* node);
  void get_update(ie_Mat* updates, const std::vector<unsigned int>& BN,
                  const QuadTreeNode* node);

  void apply_sweep_matrix(const ie_Mat& mat, ie_Mat* vec,
                          const std::vector<unsigned int>& a,
                          const std::vector<unsigned int>& b, bool transpose);
  void apply_diag_matrix(const ie_Mat& mat, ie_Mat* vec,
                         const std::vector<unsigned int>& range);
  void apply_diag_inv_matrix(const ie_Mat& mat, ie_Mat* vec,
                             const std::vector<unsigned int>& range);

  void skeletonize(const Kernel& K, QuadTree* tree);

  void sparse_matvec(const Kernel& K, const QuadTree& tree, const ie_Mat& x,
                     ie_Mat* b);
  void solve(const Kernel& K, const QuadTree& tree, ie_Mat* x, const ie_Mat& b);

  void check_factorization_against_kernel(const Kernel& kernel,
                                          QuadTree* tree);

  void populate_all_active_boxes(QuadTree* tree);
  void populate_active_box(QuadTreeNode* node);

  void make_id_mat(const Kernel& K, ie_Mat* pxy, const QuadTree* tree,
                   const QuadTreeNode* node);
  void make_proxy_mat(const Kernel& kernel, ie_Mat* pxy, double cntr_x,
                      double cntr_y, double r,
                      const QuadTree* tree,
                      const std::vector<unsigned int>& box_indices);
 
  void set_rs_ranges(InteractionLists* src_dof_lists,
                     const std::vector<unsigned int>& prm,
                     unsigned int sk, unsigned int rd);
  void set_skelnear_range(InteractionLists* src_dof_lists);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_TOOLS_IE_SOLVER_TOOLS_H_
