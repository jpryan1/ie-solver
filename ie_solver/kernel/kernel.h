// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_KERNEL_KERNEL_H_
#define IE_SOLVER_KERNEL_KERNEL_H_

#include <vector>
#include "ie_solver/boundaries/boundary.h"
#include "ie_solver/ie_mat.h"
#include "ie_solver/io/io.h"
#include "ie_solver/quadtree/quadtree.h"

#define NUM_PROXY_POINTS 128
#define RADIUS_RATIO 1.5

namespace ie_solver {

struct Dof {
  Vec2 point;
  Vec2 normal;
  double curvature;
  double weight;
  bool is_boundary = true;
};  // struct


struct Kernel {
  static int IMPROVE_CONDITION;

  Kernel(int solution_dimension_, int domain_dimension_,
         ie_solver_config::Pde pde_, Boundary* boundary_,
         std::vector<double> domain_points_) :
    solution_dimension(solution_dimension_),
    domain_dimension(domain_dimension_),
    pde(pde_),
    boundary(boundary_),
    domain_points(domain_points_) {}

  int solution_dimension, domain_dimension;
  ie_solver_config::Pde pde;
  Boundary * boundary;
  std::vector<double> domain_points;
  // double get(unsigned int i, unsigned int j) const;
  // double forward_get(unsigned int i, unsigned int j) const;
  ie_Mat get(const Dof & a, const Dof & b) const;

  ie_Mat stokes_kernel(const Dof & a, const Dof & b) const;
  ie_Mat laplace_kernel(const Dof & a, const Dof & b) const;
  ie_Mat laplace_neumann_kernel(const Dof & a, const Dof & b) const;
  ie_Mat laplace_neumann_kernel_forward(const Dof & a, const Dof & b) const;

  ie_Mat operator()(const std::vector<unsigned int> & I_,
                    const std::vector<unsigned int> & J_,
                    double * timing = nullptr) const;

  // ie_Mat forward_get(const std::vector<unsigned int>& I_,
  //                    const std::vector<unsigned int>& J_) const;
  ie_Mat fast_laplace_get(const std::vector<unsigned int> & I_,
                          const std::vector<unsigned int> & J_,
                          double * timing) const;
  ie_Mat fast_laplace_neumann_get(const std::vector<unsigned int> & I_,
                                  const std::vector<unsigned int> & J_,
                                  double * timing) const;
  ie_Mat fast_stokes_get(const std::vector<unsigned int> & I_,
                         const std::vector<unsigned int> & J_,
                         double * timing) const;

  ie_Mat get_id_mat(const QuadTree * tree,
                    const QuadTreeNode * node,
                    bool strong_admissibility) const;
  ie_Mat get_proxy_mat(double cntr_x, double cntr_y,
                       double r, const QuadTree * tree,
                       const std::vector<unsigned int> & box_inds) const;

  ie_Mat fast_laplace_proxy_get(const std::vector<double> & pxy_p,
                                const std::vector<double> & pxy_n,
                                double pxy_w,
                                const std::vector<unsigned int> & box_inds)
  const;
  ie_Mat fast_laplace_neumann_proxy_get(const std::vector<double> & pxy_p,
                                        const std::vector<double> & pxy_n,
                                         double pxy_w,
                                        const std::vector<unsigned int> &
                                        box_inds) const;
  ie_Mat fast_stokes_proxy_get(const std::vector<double> & pxy_p,
                               const std::vector<double> & pxy_n,
                              double pxy_w,
                               const std::vector<unsigned int> & box_inds)
  const;
};  // struct

}  // namespace ie_solver

#endif  // IE_SOLVER_KERNEL_KERNEL_H_
