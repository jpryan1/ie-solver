// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_LINEAR_SOLVE_H_
#define IE_SOLVER_LINEAR_SOLVE_H_

#include <vector>
#include <string>
#include "ie_solver/io/io.h"
#include "ie_solver/quadtree/quadtree.h"
#include "ie_solver/boundaries/boundary.h"

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

void get_domain_points(unsigned int domain_size, std::vector<double>* points,
                       double min, double max);

void check_factorization_against_kernel(const Kernel& kernel, QuadTree* tree);


}  // namespace ie_solver

#endif  // IE_SOLVER_LINEAR_SOLVE_H_
