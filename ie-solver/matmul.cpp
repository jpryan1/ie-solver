#include "common.h"

namespace ie_solver{

void Matmul::ie_gemv(CBLAS_TRANSPOSE trans0, double alpha, const ie_Mat& A, 
	const ie_Mat& x, double beta, ie_Mat& b){

	assert(A.height()*A.width()*x.height()*x.width() !=0);
	assert(x.width() == 1);
	cblas_dgemv(CblasColMajor, trans0, A.height(), A.width(), alpha, A.mat, 
		A.lda_, x.mat, 1, beta, b.mat, 1);
}

void Matmul::ie_gemm(CBLAS_TRANSPOSE trans0, CBLAS_TRANSPOSE trans1, 
	double alpha, const ie_Mat& A, const ie_Mat& B, double beta, ie_Mat& C){

	assert(A.height()*B.height()*A.width()*B.width()!=0);
	assert(A.mat != nullptr && B.mat != nullptr);

	unsigned int k = (trans0 == CblasTrans) ? A.height() : A.width();
	cblas_dgemm(CblasColMajor, trans0, trans1, C.height(), C.width(), 
		k, alpha, A.mat, A.height(), B.mat, B.height(), beta, C.mat, 
		C.height());
}


} // namespace