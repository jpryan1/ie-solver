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

#include <cblas.h>
#include <lapacke.h>


#define NORMAL_ 0
#define TRANSPOSE_ 1






struct Clock{
	double time;
	double elapsed_;
	void tic(){
		time = omp_get_wtime();
	}
	void toc(){
		elapsed_ += (omp_get_wtime()-time);
	}
	void elapsed(const char* s){
		printf("%s: %f seconds elapsed\n", s, elapsed_);
	}
	// For main.cpp
	void toc(const char* s){
		printf("%s: %f seconds\n", s, (omp_get_wtime() - time));
	}
};





struct Box {
// Note, maybe redundant and skel could be Range<int>s, that would make things faster maybe
	std::vector<int> box_range;
	std::vector<int> redundant_range;
	std::vector<int> skel_range;
	std::vector<int> near_range;
	std::vector<int> far_range;
	std::vector<int> skelnear_range;

	std::vector<int> p;

};




struct Vec2{

	double a[2];
	Vec2(double* arr){
		a[0]=arr[0];
		a[1]=arr[1];
	}
	Vec2(){
		a[0]=0;
		a[1]=0;
	}
	Vec2(double m, double n){
		a[0]=m;
		a[1]=n;
	}
	double norm(){
		return sqrt(a[0]*a[0] + a[1]*a[1]);
	}

	double dot(const Vec2& o){
		return a[0]*o.a[0] + a[1]*o.a[1];
	}

	Vec2 operator-(const Vec2 &o){
		return Vec2(a[0]-o.a[0], a[1]-o.a[1]);
	}
	Vec2 operator*(const double d){
		return Vec2(a[0]*d, a[1]*d);
	}
	void print(){
		std::cout<<"Vec2: "<<a[0]<<" "<<a[1]<<std::endl;
	}


};

struct ie_Mat{

	// storage is column major, so by default lda is the height
	// TODO preprocessor to get rid of ugly std::vectors?
	double *mat;
	int lda_, height_, width_;
	int dynamic;
	int is_stokes;
	std::vector<double>* points;
	std::vector<double>* normals;
	std::vector<double>* curvatures;
	std::vector<double>* weights;

	double diag_00, diag_11, diag_01;
	double scale = 1.0/(2.0*M_PI);


	ie_Mat(){
		dynamic   = 0;
		is_stokes = 0;
		mat       = NULL;
		lda_      = 0;
		height_   = 0;
		width_    = 0;
	}
	~ie_Mat(){
		if(mat) delete[] mat;
	}
	ie_Mat(int h, int w){
		assert(h>0 && w>0);
		dynamic   = 0;
		is_stokes = 0;
		lda_      = h;
		height_   = h;
		width_    = w;
		mat       = new double[height_*width_];
		memset(mat, 0, height_*width_*sizeof(double));
	}



	// This function stores the DoF data,  and calculates the diagonals of the mat
	void load(std::vector<double>* p, std::vector<double>* n, 
		std::vector<double>* c, std::vector<double>* w);





	void resize(int h, int w);


	void copy(ie_Mat& copy);

	void rand_vec(int dofs);
	

	double stokes_kernel(int i, int j);

	double laplace_kernel(int i, int j);

	double get(int i, int j);


	void set(int i, int j, double a);
	
	void addset(int i, int j, double a);




	void inverse();

	void left_multiply_inverse(ie_Mat& K, ie_Mat& U);
	void right_multiply_inverse(ie_Mat& K, ie_Mat& L);

	void set_submatrix(std::vector<int> I_, std::vector<int> J_, ie_Mat& A);
	

	void transpose(ie_Mat& A){
	

		if(height_ != A.width_ || width_ != A.height_){
			if(A.mat) delete[] A.mat;
			A.lda_    = width_;
			A.height_ = width_;
			A.width_  = height_;
			A.mat     = new double[height_*width_];
		}

		for(int i=0; i<height_; i++){
			for(int j=0; j<width_; j++){
				A.set(j,i, get(i,j));
			}
		}

		//Transpose(mat, A.mat);
	}
	double norm2(){
		double sum=0;
		for(int i=0; i<height_; i++){
			for(int j=0; j<width_; j++){
				sum += pow( get(i,j), 2);
			}
		}
		return sqrt(sum);
	}

	




	int id( std::vector<int>& p, ie_Mat& Z, double tol);


	int height(){
		return height_;
	}
	int width(){
		return width_;
	}

	ie_Mat& operator-=(const ie_Mat& o);

	ie_Mat& operator+=(const ie_Mat& o);


	ie_Mat& operator*=(const double o);


	ie_Mat& operator()(std::vector<int> I_, std::vector<int> J_);


	void operator=(const ie_Mat& copy);

	void print();

};



class Matmul{
public:

	static void ie_gemv(int trans0, double alpha, ie_Mat& A, ie_Mat& x, double beta, ie_Mat& b){


		
	//assumption is that b is the right size
	//please fix this stupid if/else nonsense
	if(trans0){
		cblas_dgemv(CblasColMajor, CblasTrans, A.height(), A.width(), alpha, A.mat, A.lda_, x.mat, 1, beta, b.mat, 1);
	}else{
		cblas_dgemv(CblasColMajor, CblasNoTrans, A.height(), A.width(), alpha, A.mat, A.lda_, x.mat, 1, beta, b.mat, 1);
	}


	}
	static void ie_gemm(int trans0, int trans1, double alpha, ie_Mat& A, ie_Mat& B, double beta, ie_Mat& C){
	//assumption is that b is the right size
//	please fix this stupid if/else nonsense
		if(trans0){
			if(trans1){
				cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, C.height(), C.width(), A.height(),
				alpha, A.mat, A.height(), B.mat, B.height(), beta, 
				C.mat, C.height());
			}
			else{
				cblas_dgemm(CblasColMajor, CblasTrans,CblasNoTrans, C.height(), C.width(), A.height(),
				alpha, A.mat, A.height(), B.mat, B.height(), beta, 
				C.mat, C.height());
			}
		}else{
			if(trans1){
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, C.height(), C.width(), A.width(),
				alpha, A.mat, A.height(), B.mat, B.height(), beta, 
				C.mat, C.height());
			}
			else{
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, C.height(), C.width(), A.width(),
				alpha, A.mat, A.height(), B.mat, B.height(), beta, 
				C.mat, C.height());
			
			}
		}
	}

};































#endif
