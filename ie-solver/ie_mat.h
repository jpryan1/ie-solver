#ifndef _IE_MAT_H_
#define _IE_MAT_H_

#include <vector>
#include <cblas.h>

#define NORMAL CblasNoTrans
#define TRANSPOSE CblasTrans

namespace ie_solver{


struct ie_Mat{

	// storage is column major, so by default lda is the height.
	double *mat;
	unsigned int lda_, height_, width_;

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
	void set_submatrix(int row_s, int row_e, int col_s, int col_e, ie_Mat& A);
	
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
	void rand_vec(unsigned  dofs);
	
	void left_multiply_inverse(const ie_Mat& K, ie_Mat& U) const;
	void right_multiply_inverse(const ie_Mat& K, ie_Mat& L) const;
	int id(std::vector<unsigned int>& p, ie_Mat& Z, double tol) const;

	void print() const;

	static void gemv(CBLAS_TRANSPOSE trans0, double alpha, const ie_Mat& A, 
		const ie_Mat& x, double beta, ie_Mat& b);
	static void gemm(CBLAS_TRANSPOSE trans0, CBLAS_TRANSPOSE trans1, 
		double alpha, const ie_Mat& A, const ie_Mat& B, double beta, ie_Mat& C);
};

} // namespace ie_solver

#endif