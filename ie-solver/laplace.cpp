#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <omp.h>

#include "common.h"
#include "mat.h"
#include "ie.h"
#include <El.hpp>

#define DEFAULT_ID_TOL 0.001
#define DEFAULT_BLOCKSIZE 10
#define DEFAULT_DOFs 1000
using namespace El;

//TODO put everything in a class so the functions can have access to shared variables for things like blocksize


void Laplace_Kernel(Matrix<double>& K, int DOFs){

}

void Laplace_Tree(Tree& tree){

}

void Laplace_Boundary(Matrix<double> f){

}








int main(int argc, char** argv){
	

	int verbosity = 0;
	int blocksize = DEFAULT_BLOCKSIZE;
	int dofs = DEFAULT_DOFs;
	double id_tol = DEFAULT_ID_TOL;
	int c;
	while ((c = getopt (argc, argv, "b:d:i:vh")) != -1) {
		switch (c)
		{
		case 'b':
			blocksize = std::strtol(optarg, NULL, 0);
			break;
		case 'd':
			dofs = std::strtol(optarg, NULL, 0);
			break;
		case 'i':
			id_tol = std::strtod(optarg, NULL);
			break;
		case 'v':
			verbosity++;
			break;
		case 'h':
			printf("IE-Solver 0.1a\nOptions: bdi [vh]");
			return 0;
		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
						"Unknown option character `\\x%x'.\n",
						optopt);
			return 1;
		default:
			abort ();
		}
	}
	// Need to push a block size prevent segfault on inverse
	PushBlocksizeStack(128);

	IE ie(blocksize, id_tol, dofs);
	ie.verbosity = verbosity;
	
	//K_ij contains the potential between electron's i and j
	Factorization F;
	F.K = Matrix<double>(dofs,dofs);

	Laplace_Kernel(F.K, dofs);
	Tree tree(dofs/blocksize);
	Laplace_Tree(tree);
	
	ie.Skeletonize(F, tree);

	Matrix<double> f, phi;

	Laplace_Boundary(f);

	ie.Conjugate_Gradient(F,tree, phi, f);

	// Matrix<double> K_D(dofs*dofs, dofs);
	// Matrix<double> u(dofs*dofs, dofs*dofs);
	// Gemv(NORMAL, 1., K_D, phi, 0., u);








	return 0;
}