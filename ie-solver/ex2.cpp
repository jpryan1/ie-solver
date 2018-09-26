#include <fstream>

#include "common.h"
#include "initialization.h"
#include "skelfac.h"
#include "quadtree.h"
#include "circle.h"
#define DEFAULT_N 128
#define DEFAULT_ID_TOL 1e-6
#define TEST_SIZE 100

//TODO put everything in a class so the functions can have access to shared variables for things like blocksize
namespace ie_solver{

LOG::LOG_LEVEL LOG::log_level_ = LOG::LOG_LEVEL::WARNING_;


void get_circle_stokes_solution(double min, double max, ie_Mat& domain, int (*out_of_shape)(Vec2& a)){
	//min and max describe the box inside which our sample is taken. 

	double total_err = 0;
	double max_err   = 0;
	double true_norm = 0;
	
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
		
		double r = sqrt(x0*x0+y0*y0);

		domain.set(2*i  , 0, -4*y0/(1));
		domain.set(2*i+1, 0,  4*x0/(1));
		
	}
	
}


double vec_norm(ie_Mat& vec){
	double sum=0;
	for(int i=0; i<vec.height(); i++){
		sum += pow(vec.get(i, 0),2);
	}return sqrt(sum);
}


void stokes_integral_solve(int N, int verbosity, double id_tol,
	void (*make_shape) (int, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&),
	int (*out_of_shape)(Vec2& a)){

	int timing = 0;

	if(!timing && N>10000){
		printf("Turn down N or disable accuracy checking please\n");
		return;
	}
	omp_set_num_threads(1);

	Clock clock;
	bool is_stokes = true;


	std::vector<double> points, normals, curvatures, weights;
	make_shape(N, points, normals, curvatures, weights);
	int dofs = points.size()/2;
	
	
	printf("N= %d\n", dofs);
	printf("e= %.10f\n", id_tol);

	//clock.tic();
	QuadTree quadtree;
	quadtree.initialize_tree(points, is_stokes); 
	//clock.toc("Tree make");

	Skelfac skelfac(id_tol, points, normals, weights, is_stokes);
	skelfac.verbosity = verbosity;
	

	ie_Mat K(1,1);
	K.is_stokes = is_stokes;
	K.is_dynamic = true;
	K.load(&points, &normals, &curvatures, &weights);
	
	//ie_Mat K_copy(2*dofs,2*dofs);
	ie_Mat K_copy, K_domain, true_domain;
	Initialization stokes;
	//stokes.Stokes_InitializeKernel(K_copy, points, normals, curvatures, weights);

	if(!timing){
		K_copy = ie_Mat(2*dofs, 2*dofs);
		stokes.Stokes_InitializeKernel(K_copy, points, normals, curvatures, weights);

		K_domain = ie_Mat(2*TEST_SIZE*TEST_SIZE, 2*dofs);
		stokes.Stokes_InitializeDomainKernel(K_domain, points, normals, weights, 
			quadtree.min, quadtree.max, TEST_SIZE, out_of_shape);

		true_domain = ie_Mat(2*TEST_SIZE*TEST_SIZE, 1);
		get_circle_stokes_solution(quadtree.min, quadtree.max, true_domain, out_of_shape);
	}

	ie_Mat f(2*dofs, 1);
 	stokes.Stokes_InitializeBoundary(f, normals); //notice here we are passing the normals 
 		//since the flow will just be unit tangent to the boundary. 

	ie_Mat rand_vec(2*dofs,1);
	rand_vec.rand_vec(2*dofs);

	//	STEP 1 - DENSE MATVEC
	ie_Mat result_r, result_f;
	if(!timing){
		result_r = ie_Mat(2*dofs,1);
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
 	ie_Mat phi(2*dofs, 1);
 	clock.tic();
	skelfac.Solve(K, quadtree, phi, f);
	clock.toc("Solve");
	if(!timing){
		result_f = ie_Mat(2*dofs,1);
		Matmul::ie_gemv(NORMAL_, 1., K_copy, phi, 0., result_f);
		result_f -= f;
		double lser = vec_norm(result_f)/vec_norm(f);
		printf("Solve Error: %.10f\n", lser);

		//	STEP 5 - BIE SOLVE AND OUTPUT
		ie_Mat domain(2*TEST_SIZE*TEST_SIZE, 1);
		Matmul::ie_gemv(NORMAL_, 1., K_domain, phi, 0., domain);

		std::ofstream output;
		output.open("stokes.txt");

		if(output.is_open()){
			for(int i=0; i<domain.height(); i+=2){
				output<<domain.get(i, 0)<<" "<<domain.get(i+1,0)<<std::endl;
			}
			output.close();
		}else{
			printf("Failed to open output file!\n");
		}

		// STEP 6 - CHECK AGAINST TRUE ANSWER TO PHYSICAL PROBLEM
		 domain-=true_domain;
		 double der = vec_norm(domain)/vec_norm(true_domain);
		 printf("\n\nError in Solution: %.10f\n", der);
	}
}

} // namespace

int main(int argc, char** argv){

	int verbosity = 0;
	double id_tol = DEFAULT_ID_TOL;
	int N         = DEFAULT_N;
	//Set parameters based on arguments

	//boundary_integral_solve(verbosity, id_tol, squiggly, out_of_squiggly);
	ie_solver::stokes_integral_solve(N, verbosity, id_tol, ie_solver::circle, 
		ie_solver::out_of_circle);

	return 0;
}


