#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <string>
#include <cblas.h>
#include <lapacke.h>

#include "vec2.h"
#include "log.h"
#include "clock.h"
#include "interaction_lists.h"

#define NORMAL CblasNoTrans
#define TRANSPOSE CblasTrans

namespace ie_solver{


struct ie_Mat{

	// storage is column major, so by default lda is the height.
	double *mat;
	unsigned int lda_, height_, width_;
	bool is_stokes, is_dynamic;
	const std::vector<double> *points, *normals, *curvatures, *weights;

	double diag_00, diag_11, diag_01;
	double scale = 1.0/(2.0*M_PI);

	ie_Mat();
	~ie_Mat();
	ie_Mat(unsigned int h, unsigned int w);
	ie_Mat& operator=(const ie_Mat& copy);

	void resize(unsigned int h, unsigned int w);
	void copy(ie_Mat& copy) const;

	double get(unsigned int i, unsigned int j) const;
	void set(unsigned int i, unsigned int j, double a);
	void addset(unsigned int i, unsigned int j, double a);
	void set_submatrix(const std::vector<unsigned int>& I_, 
		const std::vector<unsigned int>& J_, ie_Mat& A);
	void transpose(ie_Mat& A) const;

	unsigned int height() const;
	unsigned int width() const;

	ie_Mat& operator-=(const ie_Mat& o);
	ie_Mat& operator+=(const ie_Mat& o);
	ie_Mat& operator*=(double o);
	ie_Mat operator()(const std::vector<unsigned int>& I_, 
		const std::vector<unsigned int>& J_) const;

	double norm2() const;

	// This function stores the DoF data, and calculates the diagonals of the 
	// mat.
	void load(const std::vector<double>* p, const std::vector<double>* n, 
		const std::vector<double>* c, const std::vector<double>* w);

	void rand_vec(unsigned  dofs);
	
	double stokes_kernel(unsigned int i, unsigned int j) const;
	double laplace_kernel(unsigned int i, unsigned int j) const;

	void left_multiply_inverse(const ie_Mat& K, ie_Mat& U) const;
	void right_multiply_inverse(const ie_Mat& K, ie_Mat& L) const;
	int id(std::vector<unsigned int>& p, ie_Mat& Z, double tol) const;

	void print() const;
};


struct Matmul{
	static void ie_gemv(CBLAS_TRANSPOSE trans0, double alpha, const ie_Mat& A, 
		const ie_Mat& x, double beta, ie_Mat& b);
	static void ie_gemm(CBLAS_TRANSPOSE trans0, CBLAS_TRANSPOSE trans1, 
		double alpha, const ie_Mat& A, const ie_Mat& B, double beta, ie_Mat& C);
};

} // namespace ie_solver

#endif