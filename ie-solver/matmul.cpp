#include "common.h"

namespace ie_solver{

void Matmul::ie_gemv(int trans0, double alpha, const ie_Mat& A, 
	const ie_Mat& x, double beta, ie_Mat& b){
	//assumption is that b is the right size
	//please fix this stupid if/else nonsense
	if(trans0){
		cblas_dgemv(CblasColMajor, CblasTrans, A.height(), A.width(), alpha,
		 A.mat, A.lda_, x.mat, 1, beta, b.mat, 1);
	}else{
		cblas_dgemv(CblasColMajor, CblasNoTrans, A.height(), A.width(), 
			alpha, A.mat, A.lda_, x.mat, 1, beta, b.mat, 1);
	}


}
void Matmul::ie_gemm(int trans0, int trans1, double alpha,
	const ie_Mat& A, const ie_Mat& B, double beta, ie_Mat& C){
	//assumption is that b is the right size
	//	please fix this stupid if/else nonsense
	if(trans0){
		if(trans1){
			cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, C.height(), 
				C.width(), A.height(),
			alpha, A.mat, A.height(), B.mat, B.height(), beta, 
			C.mat, C.height());
		}
		else{
			cblas_dgemm(CblasColMajor, CblasTrans,CblasNoTrans, C.height(), 
				C.width(), A.height(),
			alpha, A.mat, A.height(), B.mat, B.height(), beta, 
			C.mat, C.height());
		}
	}else{
		if(trans1){
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, C.height(),
			 C.width(), A.width(),
			alpha, A.mat, A.height(), B.mat, B.height(), beta, 
			C.mat, C.height());
		}
		else{
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
				C.height(), C.width(), A.width(), alpha, A.mat, A.height(), 
				B.mat, B.height(), beta, C.mat, C.height());
		
		}
	}
}


} // namespace