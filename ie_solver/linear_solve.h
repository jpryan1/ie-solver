// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_LINEAR_SOLVE_H_
#define IE_SOLVER_LINEAR_SOLVE_H_

#include <vector>
#include <string>
#include "ie_solver/io/io.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/boundaries/boundary.h"
#include "ie_solver/kernel/kernel.h"
#include "ie_solver/skel_factorization/skel_factorization.h"


#define TIMING_ITERATIONS 5

namespace ie_solver {

ie_Mat boundary_integral_solve(const ie_solver_config& config,
                               QuadTree* quadtree,
                               const std::vector<double>& domain_points);
void bie_time_trial(const ie_solver_config & config,
                    QuadTree * quadtree, double* avg_skel_time,
                    double* avg_solve_time);
// void get_circle_stokes_solution(double min, double max, ie_Mat& domain,
//     bool (*is_in_domain)(Vec2& a));
ie_Mat initialize_U_mat(const ie_solver_config::Pde pde,
                        const std::vector<Hole>& holes,
                        const std::vector<double>& tgt_points);
ie_Mat initialize_Psi_mat(const ie_solver_config::Pde pde,
                          const std::vector<Hole>& holes, Boundary * boundary);
void get_domain_points(int domain_size, std::vector<double>* points,
                       double x_min, double x_max, double y_min, double y_max);
void linear_solve(const SkelFactorization& skel_factorization,
                  const QuadTree& quadtree, const ie_Mat& f, ie_Mat* mu,
                  ie_Mat* alpha = nullptr);

// void linear_solve(const SkelFactorization& skel_factorization,
//                   const QuadTree& quadtree, const ie_Mat& f, ie_Mat* mu,
//                   double* c);


void schur_solve(const SkelFactorization & skel_factorization,
                 const QuadTree & quadtree, const ie_Mat & U,
                 const ie_Mat & Psi,
                 const ie_Mat & f, const ie_Mat & K_domain,
                 const ie_Mat & U_forward,  ie_Mat * solution);

void check_factorization_against_kernel(const Kernel& kernel, QuadTree* tree);


}  // namespace ie_solver

#endif  // IE_SOLVER_LINEAR_SOLVE_H_
