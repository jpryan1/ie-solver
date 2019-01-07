#include <fstream>
#include <memory>
#include <omp.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <cassert>
#include "ie_mat.h"
#include "initialization.h"
#include "ie-solver/tools/ie_solver_tools.h"
#include "quadtree.h"
#include "kernel.h"
#include "log.h"
#include "helpers.h"

#define TIMING_ITERATIONS 5

namespace ie_solver{

void boundary_integral_solve(const ie_solver_config& config, 
	std::vector<double>* skel_times = nullptr){

	bool is_stokes = (config.pde == ie_solver_config::STOKES);
	bool is_time_trial = (skel_times != nullptr);
	int N = config.N;
	double id_tol = config.id_tol;

	// // TODO insert comment here explaining why this is necessary
	// if(!timing && N > 10000){
	// 	printf("Turn down N or disable accuracy checking please\n");
	// 	return;
	// }
	config.boundary->initialize(N);

	if(!is_time_trial){
		write_boundary_to_file(config.boundary->points);
	}
	// TODO why not just N here?
	int dofs = config.boundary->points.size() / 2;
	QuadTree quadtree;
	quadtree.initialize_tree(config.boundary.get(), is_stokes); 
	if(!is_time_trial){
		quadtree.write_quadtree_to_file();
	}

	// Consider making init instead of constructor for readability
	bool strong_admissibility = 
		(config.admissibility == ie_solver_config::STRONG);
	IeSolverTools ie_solver_tools(id_tol, strong_admissibility, is_stokes);

	Kernel kernel;

	kernel.load(config.boundary.get(), is_stokes);

	// Now we calculate the solution to the PDE inside the domain by setting
	// up the relevant linear system. 

	int dim = is_stokes ? 2 : 1;
	
	std::vector<double> domain_points;
	get_domain_points(domain_points, quadtree.min, quadtree.max);
	
	Initialization init;
	
	ie_Mat f(dim*dofs, 1);
	// TODO get rid of these damn if(is_stokes) statements, push them to the 
	// functions
 	if(is_stokes){
 		init.Stokes_InitializeBoundary(f, config.boundary->normals); 
 		// notice here we are passing the normals since the flow will just be 
 		// unit tangent to the boundary. 
 	}else{
 		init.InitializeBoundary(f, config.boundary->points);
 	}
	
 	ie_Mat phi(dim*dofs, 1);

	if(is_time_trial){
		double elapsed = 0;
		for(int i=0; i < TIMING_ITERATIONS; i++){
			double start = omp_get_wtime();
			ie_solver_tools.skeletonize(kernel, quadtree);
			double end = omp_get_wtime();
			elapsed += (end-start);

			quadtree.reset();
		}
		skel_times->push_back(elapsed/TIMING_ITERATIONS);
		return;
	}
	
	ie_solver_tools.skeletonize(kernel, quadtree);

	ie_solver_tools.solve(kernel, quadtree, phi, f);

	// This will be done as a sparse mat vec in the future, for now we do 
	// dense matvec
	ie_Mat K_domain = ie_Mat(dim*TEST_SIZE*TEST_SIZE, dim*dofs);
	init.InitializeDomainKernel(K_domain, config.boundary->points, 
		config.boundary->normals, config.boundary->weights, 
 		domain_points, TEST_SIZE, config.boundary.get(), is_stokes);
	ie_Mat domain(dim*TEST_SIZE*TEST_SIZE, 1);
	ie_Mat::gemv(NORMAL, 1., K_domain, phi, 0., domain);
	write_solution_to_file(domain, domain_points, is_stokes);

	check_laplace_solution(domain, config.id_tol, domain_points, 
		config.boundary.get());
	
}

} // namespace ie_solver

int main(int argc, char** argv){
	
	// TODO allow for command line args for setting parameters

	ie_solver::ie_solver_config config;
	if(!ie_solver::parse_input_into_config(argc, argv, config)){
		return 1;
	}
	if(!config.scaling){
		ie_solver::boundary_integral_solve(config);
	}

	if(config.scaling){
		std::vector<double> n_times, eps_times;
		int scale_n[] = {5000, 6000, 7000, 8000, 9000};
		double scale_eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};

		for(int n : scale_n){
			config.N = n;
			config.id_tol = DEFAULT_ID_TOL;
			ie_solver::boundary_integral_solve(config, &n_times);
		}

		for(double eps : scale_eps){
			config.N = DEFAULT_NUM_DISCRETIZATION_POINTS*10;
			config.id_tol = eps;
			ie_solver::boundary_integral_solve(config, &eps_times);
		}

		ie_solver::write_times_to_files(scale_n, n_times, scale_eps, eps_times);
	}

	return 0;
}

