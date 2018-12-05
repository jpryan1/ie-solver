#include "common.h"

namespace ie_solver{

void Matmul::ie_gemv(CBLAS_TRANSPOSE trans0, double alpha, const ie_Mat& A, 
	const ie_Mat& x, double beta, ie_Mat& b){
	assert(A.height()*A.width()*x.height()*x.width() !=0 && 
		"gemv needs positive dimensions only.");
	assert(x.width() == 1 && "gemv only works on a column vector.");
	assert(A.mat != nullptr && "gemv fails on null A mat.");
	assert(x.mat != nullptr && "gemv fails on null x mat.");
	assert(b.mat != nullptr && "gemv fails on null b mat.");


	cblas_dgemv(CblasColMajor, trans0, A.height(), A.width(), alpha, A.mat, 
		A.lda_, x.mat, 1, beta, b.mat, 1);
}

void Matmul::ie_gemm(CBLAS_TRANSPOSE trans0, CBLAS_TRANSPOSE trans1, 
	double alpha, const ie_Mat& A, const ie_Mat& B, double beta, ie_Mat& C){
	assert(A.height()*B.height()*A.width()*B.width()!=0 &&
		"gemm needs positive dimensions only.");
	assert(A.mat != nullptr && "gemm fails on null A mat.");
	assert(B.mat != nullptr && "gemm fails on null B mat.");
	assert(C.mat != nullptr && "gemm fails on null C mat.");
	assert(A.mat != nullptr && "gemm fails on null A mat.");
	assert(B.mat != nullptr && "gemm fails on null B mat.");
	assert(C.mat != nullptr && "gemm fails on null C mat.");

	unsigned int k = (trans0 == CblasTrans) ? A.height() : A.width();
	cblas_dgemm(CblasColMajor, trans0, trans1, C.height(), C.width(), 
		k, alpha, A.mat, A.height(), B.mat, B.height(), beta, C.mat, 
		C.height());
}


} // namespace