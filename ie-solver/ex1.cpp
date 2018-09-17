
#include <fstream>

#include "common.h"
#include "initialization.h"
#include "skelfac.h"
#include "quadtree.h"
#include "circle.h"
#define DEFAULT_N 128
#define DEFAULT_ID_TOL 1e-6
#define TEST_SIZE 50


//TODO put everything in a class so the functions can have access to shared variables for things like blocksize









double vec_norm(ie_Mat& vec){
	double sum=0;
	for(int i=0; i<vec.height(); i++){
		sum += pow(vec.get(i, 0),2);
	}return sqrt(sum);
}


void boundary_integral_solve(int N, int verbosity, double id_tol,
	void (*make_shape) (int, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&),
	int (*out_of_shape)(Vec2& a)){
	


	int timing = 1;

	if(!timing && N>10000){
		printf("Turn down N or disable accuracy checking please\n");
		return;
	}
	omp_set_num_threads(1);

	Clock clock;
	int is_stokes = 0;


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
	K.dynamic = 1;
	K.load(&points, &normals, &curvatures, &weights);
	

	ie_Mat K_copy, K_domain, true_domain;
	Initialization init;
		

	if(!timing){
		K_copy = ie_Mat(dofs, dofs);
		init.InitializeKernel(K_copy, points, normals, curvatures, weights);
		
		K_domain = ie_Mat(TEST_SIZE*TEST_SIZE, dofs);
		init.InitializeDomainKernel(K_domain, points, normals, weights, 
	 		quadtree.min, quadtree.max, TEST_SIZE, out_of_shape);
		
		true_domain = ie_Mat(TEST_SIZE*TEST_SIZE, 1);
	 	init.DomainSolution(true_domain, TEST_SIZE, quadtree.min, quadtree.max, out_of_shape);
	}


	ie_Mat f(dofs, 1);
 	init.InitializeBoundary(f, points);


	
	
 	ie_Mat rand_vec(dofs,1);
	rand_vec.rand_vec( dofs);
	


	ie_Mat result_r, result_f;
	//	STEP 1 - DENSE MATVEC (if not timing)
	if(!timing){
		result_r = ie_Mat(dofs,1);
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
 	ie_Mat phi(dofs, 1);
 	clock.tic();
	skelfac.Solve(K, quadtree, phi, f);
	clock.toc("Solve");
	if(!timing){
		result_f = ie_Mat(dofs,1);
		Matmul::ie_gemv(NORMAL_, 1., K_copy, phi, 0., result_f);
		result_f -= f;
		double lser = vec_norm(result_f)/vec_norm(f);
		printf("Solve Error: %.10f\n", lser);




	//	STEP 5 - BIE SOLVE AND OUTPUT
		 ie_Mat domain(TEST_SIZE*TEST_SIZE, 1);
		 Matmul::ie_gemv(NORMAL_, 1., K_domain, phi, 0., domain);

		std::ofstream output;
		output.open("laplace.txt");

		

		if(output.is_open()){
			for(int i=0; i<domain.height(); i++){
				output<<domain.get(i, 0)<<std::endl;
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





int main(int argc, char** argv){


	int verbosity = 0;
	double id_tol = DEFAULT_ID_TOL;
	int c;
	int N         = DEFAULT_N;
	//Set parameters based on arguments

	//boundary_integral_solve(verbosity, id_tol, squiggly, out_of_squiggly);
	boundary_integral_solve(N, verbosity, id_tol, circle, out_of_circle);

	return 0;
}


