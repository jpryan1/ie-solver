#include <fstream>

#include "common.h"
#include "initialization.h"
#include "skelfac.h"
#include "quadtree.h"
#include "circle.h"

#define DEFAULT_NUM_DISCRETIZATION_POINTS 128
#define DEFAULT_ID_TOL 1e-10
#define TEST_SIZE 50

namespace ie_solver{

LOG::LOG_LEVEL LOG::log_level_ = LOG::LOG_LEVEL::WARNING_;

// TODO underscores on member variables

// TODO isn't this a common function? why not put it there?
double vec_norm(ie_Mat& vec){
	double sum=0;
	for(unsigned int i=0; i<vec.height(); i++){
		sum += pow(vec.get(i, 0),2);
	}
	return sqrt(sum);
}

void get_circle_stokes_solution(double min, double max, ie_Mat& domain, int (*out_of_shape)(Vec2& a)){
	//min and max describe the box inside which our sample is taken. 

	// double total_err = 0;
	// double max_err   = 0;
	// double true_norm = 0;
	
	for(int i=0; i<TEST_SIZE*TEST_SIZE; i++){

		//find out the point corresponding to this index
		double x0 = i/TEST_SIZE;
		x0 = min + (x0*(max-min))/TEST_SIZE;
		double y0 = i%TEST_SIZE;
		y0 = min + (y0*(max-min))/TEST_SIZE;
		Vec2 v(x0,y0);

		if(out_of_shape(v)){
			continue;
		}
		//now we get a radius associated with this index
		x0-=0.5;
		y0-=0.5;
		
		// double r = sqrt(x0*x0+y0*y0);

		domain.set(2*i  , 0, -4*y0/(1));
		domain.set(2*i+1, 0,  4*x0/(1));
		
	}
	
}


void boundary_integral_solve(int N, double id_tol, void (*make_shape) 
	(int, std::vector<double>&, std::vector<double>&, std::vector<double>&, 
		std::vector<double>&), int (*out_of_shape)(Vec2& a), bool is_stokes, 
	int verbosity){

	// TODO replace Clock with more specialized timing module, do printing of 
	// stats at very end
	Clock clock;
	bool timing = false;

	// TODO insert comment here explaining why this is necessary
	if(!timing && N>10000){
		printf("Turn down N or disable accuracy checking please\n");
		return;
	}

	std::vector<double> points, normals, curvatures, weights;
	make_shape(N, points, normals, curvatures, weights);
	int dofs = points.size() / 2;

	QuadTree quadtree;
	quadtree.initialize_tree(points, is_stokes); 

	Skelfac skelfac(id_tol, points, normals, weights, is_stokes);
	skelfac.verbosity = verbosity;
	// TODO why 1,1?
	ie_Mat K(1,1);
	K.is_stokes = is_stokes;
	K.is_dynamic = true;
	K.load(&points, &normals, &curvatures, &weights);
	
	ie_Mat K_copy, K_domain, true_domain;
	Initialization init;
		

	// Dimension of the solution:
	int dim = is_stokes ? 2 : 1;
	if(!timing){
		K_copy = ie_Mat(dim*dofs, dim*dofs);
		init.InitializeKernel(K_copy, points, normals, curvatures, weights, 
			is_stokes);
		
		K_domain = ie_Mat(dim*TEST_SIZE*TEST_SIZE, dim*dofs);
		init.InitializeDomainKernel(K_domain, points, normals, weights, 
	 		quadtree.min, quadtree.max, TEST_SIZE, out_of_shape, is_stokes);
		
		true_domain = ie_Mat(dim*TEST_SIZE*TEST_SIZE, 1);
		if(is_stokes){
			get_circle_stokes_solution(quadtree.min, quadtree.max, true_domain, 
				out_of_shape);
		}else{
	 		init.DomainSolution(true_domain, TEST_SIZE, quadtree.min, quadtree.max, out_of_shape);
		}
	}

	ie_Mat f(dim*dofs, 1);
	// TODO get rid of these damn if(is_stokes) statements, push them to the functions
 	if(is_stokes){
 		init.Stokes_InitializeBoundary(f, normals); //notice here we are passing the normals 
 		//since the flow will just be unit tangent to the boundary. 
 	}else{
 		init.InitializeBoundary(f, points);
 	}
 	ie_Mat rand_vec(dim*dofs,1);
	rand_vec.rand_vec(dim*dofs);

	ie_Mat result_r, result_f;
	//	STEP 1 - DENSE MATVEC (if not timing)
	if(!timing){
		result_r = ie_Mat(dim*dofs,1);
		Matmul::ie_gemv(NORMAL_, 1., K_copy, rand_vec, 0.0, result_r);
	}

	//	STEP 2 - SKELETONIZE TIMING
	clock.tic();
	skelfac.Skeletonize(K, quadtree);
	clock.toc("Factor");

	//	STEP 3 - SPARSE MATVEC TIMING AND ERROR CHECK
 	clock.tic();
 	ie_Mat result_skel_r(rand_vec.height(), 1);
	skelfac.SparseMatVec(K, quadtree, rand_vec, result_skel_r);
	clock.toc("Sparse Mat Vec");
	
	if(!timing){
		result_skel_r -= result_r;
		double smver = vec_norm(result_skel_r)/vec_norm(result_r);
		printf("Sparse Mat Vec Error: %.10f \n", smver);
	}

	//	STEP 4 - LINEAR SOLVE AND ERROR CHECK
 	ie_Mat phi(dim*dofs, 1);
 	clock.tic();
	skelfac.Solve(K, quadtree, phi, f);
	clock.toc("Solve");
	if(!timing){
		result_f = ie_Mat(dim*dofs,1);
		Matmul::ie_gemv(NORMAL_, 1., K_copy, phi, 0., result_f);
		result_f -= f;
		double lser = vec_norm(result_f)/vec_norm(f);
		printf("Solve Error: %.10f\n", lser);

	//	STEP 5 - BIE SOLVE AND OUTPUT
		ie_Mat domain(dim*TEST_SIZE*TEST_SIZE, 1);
		Matmul::ie_gemv(NORMAL_, 1., K_domain, phi, 0., domain);

		std::ofstream output;
		// TODO this filename should depend on what pde we're solving
		output.open("ie-solver-output.txt");
		// TODO move all this nonsense to a print function in the matrix class
		if(output.is_open()){
			for(unsigned int i=0; i<domain.height(); i += dim){
				output<<domain.get(i, 0);
				if(is_stokes) output << " " << domain.get(i+1,0);
				output << std::endl;
			}
			output.close();
		}else{
			printf("Failed to open output file!\n");
		}
		// STEP 6 - CHECK AGAINST TRUE ANSWER TO PHYSICAL PROBLEM
		 domain -= true_domain;
		 double der = vec_norm(domain)/vec_norm(true_domain);
		 printf("Error in Solution: %.10f\n", der);
	}
}


} // namespace ie_solver

int main(int argc, char** argv){
	
	// TODO incorporate logging struct instead of using verbosity variable.
	int verbosity = 0;
	// TODO this should obviously be a command line arg
	bool is_stokes = true;
	double id_tol = DEFAULT_ID_TOL;
	int num_discretization_points = DEFAULT_NUM_DISCRETIZATION_POINTS;
	// TODO allow for command line args for setting parameters

	ie_solver::LOG::INFO("Testing logging");
	ie_solver::boundary_integral_solve(num_discretization_points, id_tol, 
		ie_solver::circle, ie_solver::out_of_circle, is_stokes, verbosity);

	return 0;
}

