// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_HELPERS_H_
#define IE_SOLVER_HELPERS_H_

#include <vector>
#include <string>
#include "ie-solver/ie_solver_config.h"
#include "ie-solver/boundaries/boundary.h"

#define TEST_SIZE 20
#define TIMING_ITERATIONS 5

namespace ie_solver {

double boundary_integral_solve(ie_solver_config& config,
                               std::vector<double>* skel_times = nullptr);

// void get_circle_stokes_solution(double min, double max, ie_Mat& domain,
//     bool (*is_in_domain)(Vec2& a));
void write_boundary_to_file(const std::vector<double>& points);
void write_potential_to_file();
void write_times_to_files(int* scale_n, const std::vector<double>& n_times,
                          double* scale_eps,
                          const std::vector<double>& eps_times);
void write_solution_to_file(const std::string& filename, const ie_Mat& domain,
                            const std::vector<double>&
                            domain_points, int solution_dimension);
void get_domain_points(std::vector<double>* points, double min,
                       double max);
double laplace_error(const ie_Mat& domain, double id_tol,
                     const std::vector<double>& domain_points,
                     Boundary* boundary);
double stokes_error(const ie_Mat& domain, double id_tol,
                    const std::vector<double>& domain_points,
                    Boundary* boundary);
int parse_input_into_config(int argc, char** argv, ie_solver_config* config);

}  // namespace ie_solver

#endif  // IE_SOLVER_HELPERS_H_
