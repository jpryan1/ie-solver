#ifndef _HELPERS_H_
#define _HELPERS_H_

#include <vector>
#include <memory>
#include <fstream>
#include "ie_mat.h"
#include "ie-solver/boundaries/boundary.h"
#include "ie-solver/boundaries/circle.h"
#include "ie-solver/boundaries/rounded_square.h"
#include "ie-solver/boundaries/rounded_square_with_bump.h"

#define DEFAULT_NUM_DISCRETIZATION_POINTS 1000
#define DEFAULT_ID_TOL 1e-6
#define TEST_SIZE 100

namespace ie_solver{

struct ie_solver_config{
	enum Pde {
		LAPLACE,
		STOKES
	};
	enum Admissibility {
		WEAK,
		STRONG
	};
	int N = DEFAULT_NUM_DISCRETIZATION_POINTS;
	double id_tol = DEFAULT_ID_TOL;
	Pde pde = LAPLACE;
	Admissibility admissibility = WEAK;
	std::unique_ptr<Boundary> boundary;
	bool scaling = false;
	std::ofstream n_scaling_output, error_scaling_output;

};

// void get_circle_stokes_solution(double min, double max, ie_Mat& domain, 
// 	bool (*is_in_domain)(Vec2& a));
void write_boundary_to_file(std::vector<double>& points);
void write_potential_to_file();
void write_times_to_files(int* scale_n, std::vector<double>& n_times, 
	double* scale_eps, std::vector<double>& eps_times);
void write_solution_to_file(ie_Mat& domain, std::vector<double>& domain_points,
	bool is_stokes);
void get_domain_points(std::vector<double>& points, double min, double max);
void check_laplace_solution(ie_Mat& domain, double id_tol, 
	std::vector<double>& domain_points, Boundary* boundary);
int parse_input_into_config(int argc, char** argv, ie_solver_config& config);

} // namespace ie_solver

#endif