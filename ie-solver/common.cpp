#include "common.h"

namespace ie_solver{

// Performs interpolative decomposition, and eturns number of skeleton columns. 
// Takes double /tol/, tolerance factorfor error in CPQR factorization. 
// Populates /p/ with permutation, Z with linear transformation.
int ie_Mat::id( std::vector<unsigned int>& p, ie_Mat& Z, double tol){

	ie_Mat cpy;
	this->copy(cpy);

	lapack_int pvt[width_];
	memset(pvt, 0, width_*sizeof(lapack_int));
	
	// /tau/ will contain an output from dgeqp3 that we don't need.
	double tau[width_];
	LAPACKE_dgeqp3(CblasColMajor, height_, width_, cpy.mat, lda_, pvt, tau);

	unsigned int skel = 0;
	double thresh = fabs(tol * cpy.get(0,0));
	
	for (unsigned int i = 1; i < width_; i++){
		// check if R_{i,i} / R_{0,0} < tol
		if(fabs(cpy.get(i,i)) < thresh){
			skel = i;
			break;
		}
	}
	
	if(skel == 0){
		//no compression to be done :/
		return 0;
	}

	for (unsigned int i = 0; i < width_; i++){
		p.push_back(pvt[i]-1);
	}

	int redund = width_ - skel;
	
	// set Z to be R_11^-1 R_12. Note 'U' (above diagonal) part of cp.mat
	// is the R matrix from dgeqp3.
	LAPACKE_dtrtrs(CblasColMajor, 'U', 'N', 'N', skel, redund, cpy.mat, 
		cpy.lda_, cpy.mat + cpy.lda_ * skel, cpy.lda_);
	
	std::vector<unsigned int> I_;
	std::vector<unsigned int> J_;
	for(unsigned int i = 0; i < skel; i++) I_.push_back(i);
	for(unsigned int i = skel; i < skel + redund; i++) J_.push_back(i);
	cpy(I_, J_).copy(Z);

	return skel;
}


void ie_Mat::left_multiply_inverse(ie_Mat& K, ie_Mat& U){
	// X^-1K = U
	//aka, XU = K

	//TODO insert asserts for these functions

	ie_Mat X_copy(height_, width_);
	copy(X_copy);

	ie_Mat K_copy(K.height_, K.width_);

	K.copy(K_copy);

	lapack_int ipiv[height_];
	memset(ipiv,0, height_*sizeof(lapack_int));
	
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, X_copy.height_, X_copy.width_, X_copy.mat, X_copy.lda_, ipiv );


	int status = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , X_copy.height_ ,  K_copy.width_ , X_copy.mat ,
	 X_copy.lda_ , ipiv , K_copy.mat, K_copy.lda_ );
	
	assert(status==0);
	K_copy.copy(U);
}


void ie_Mat::right_multiply_inverse(ie_Mat& K, ie_Mat& L){
	// KX^-1 = L
	// aka X_T L^T = K^T
	
	ie_Mat X_copy(height_, width_);
	transpose(X_copy);

	ie_Mat K_copy(K.width_, K.height_);

	K.transpose(K_copy);

	lapack_int ipiv[height_];
	memset(ipiv,0, height_*sizeof(lapack_int));
	
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, X_copy.height_, X_copy.width_, X_copy.mat, X_copy.lda_, ipiv );


	int err = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , X_copy.height_ ,  K_copy.width_ , X_copy.mat ,
	 X_copy.lda_ , ipiv , K_copy.mat, K_copy.lda_ );
	
	assert(err==0);


	K_copy.transpose(L);
}

// This function stores the DoF data,  and calculates the diagonals of the mat
void ie_Mat::load(std::vector<double>* p, std::vector<double>* n, 
	std::vector<double>* c, std::vector<double>* w){
	assert(is_dynamic);
	points     = p;
	normals    = n;
	curvatures = c;
	weights    = w;

	if(is_stokes){
		scale = 1/(M_PI);
	}else{
		scale = 1/(2*M_PI);
	}

	if(is_stokes){
		double avg = 0;
		for(double weight : *weights){
			avg += weight;
		}
		avg /= weights->size();

		double alpha = avg/2.0;
		double beta  = alpha*alpha;
		diag_00 = 8 * beta + 2 * beta * log(1/(2*beta)) - beta*M_PI;
		diag_01 = 0;
		diag_11 = diag_00;
	
		//printf("Diagonals are %f\n", diag_00);
	}	
}


void ie_Mat::resize(unsigned int h, unsigned int w){
	if(mat) delete[] mat;
	lda_ = h;
	height_ = h;
	width_ = w;
	mat = new double[height_*width_];
	memset(mat, 0, height_*width_*sizeof(double));
}


void ie_Mat::copy(ie_Mat& copy){

	assert(height_>0 && width_>0 && mat != NULL);
	//check if we need a resize
	if(height_ != copy.height_ || width_ != copy.width_ || lda_ != copy.lda_){
		if(copy.mat) delete[] copy.mat;
		copy.lda_    = lda_;
		copy.height_ = height_;
		copy.width_  = width_;
		copy.mat     = new double[height_*width_];
	}

	for(unsigned int i=0; i<copy.height(); i++){
		for(unsigned int j=0; j<copy.width(); j++){
			copy.set(i,j, get(i,j));
		}
	}
}

void ie_Mat::rand_vec(unsigned int dofs){
	//check if we need a resize
	if (width_ != 1 || height_ != dofs || lda_ != dofs){
		if(mat) delete[] mat;
		lda_    = dofs;
		height_ = dofs;
		width_  = 1;
		mat     = new double[height_*width_];
	}
	for(unsigned int i=0; i<dofs; i++){
		mat[i] = rand()/(0.0+RAND_MAX);
	}
}


double ie_Mat::stokes_kernel(unsigned int i, unsigned int j){
	// So this is much more awkwardly written than the function in stokes_init.cpp
	// that just writes the entire matrix at once, but the cost of writing the whole
	// matrix is just stupid. So we power through this function
	//below commented is single layer
	// // We first need to ascertain which DoFs are relevant here. 
	// int dof_i = i/2;
	// int dof_j = j/2;
	// // If the DoFs are the same, follow the singular procedure

	// if(dof_i==dof_j){
	// 	if(i!=j) return scale*diag_01;
	// 	return scale*diag_00;
	// }
	// // If they are different, figure out which of the three tensor components
	// // is desired, then return it. 

	// Vec2 x((*points)[2*dof_i], (*points)[2*dof_i+1]);
	// Vec2 y((*points)[2*dof_j], (*points)[2*dof_j+1]);
	// Vec2 r = x-y;  

	// double r0 = r.a[0];
	// double r1 = r.a[1];
	
	// if(i%2==0 && j%2==0){
	// 	return (*weights)[dof_j]*scale*(log(1.0/r.norm()) + (r0*r0/r.dot(r)));
	// }
	// else if(i%2==1 && j%2==1){
	// 	return (*weights)[dof_j]*scale*(log(1.0/r.norm()) + (r1*r1/r.dot(r)));
	// }
	// else{
	// 	return (*weights)[dof_j]*scale*(                    (r1*r0/r.dot(r)));
	// }


	int dof_i = i/2;
	int dof_j = j/2;
	
	if(dof_i==dof_j){

		double t0 = -(*normals)[2*dof_i+1];
		double t1 =  (*normals)[2*dof_i];

		if(i%2==0 && j%2==0){
			return -0.5 - 0.5*(*curvatures)[dof_i]*(*weights)[dof_i]*scale*t0*t0;
		}
		else if(i%2==1 && j%2==1){
			return -0.5 - 0.5*(*curvatures)[dof_i]*(*weights)[dof_i]*scale*t1*t1;
		}
		else{
			return -0.5 * (*curvatures)[dof_i]*(*weights)[dof_i]*scale*t0*t1;
		}

	}

	Vec2 x((*points)[2*dof_i], (*points)[2*dof_i+1]);
	Vec2 y((*points)[2*dof_j], (*points)[2*dof_j+1]);
	Vec2 r = x-y;  

	double r0 = r.a[0];
	double r1 = r.a[1];
	Vec2 n((*normals)[2*dof_j], (*normals)[2*dof_j+1]);
	double potential = (*weights)[dof_j]*scale*(r.dot(n))/(pow(r.dot(r),2));
		
		
	if(i%2==0 && j%2==0){
		return potential*r0*r0;
	}
	else if(i%2==1 && j%2==1){
		return  potential*r1*r1;
	}
	else{
		return potential*r0*r1;
	}
}

double ie_Mat::laplace_kernel(unsigned int i, unsigned int j){
	
	if(i==j){
		return 0.5 + 0.5*(*curvatures)[i]*(*weights)[i]*scale;
		
	 }
	
	Vec2 x((*points)[2*i], (*points)[2*i+1]);
	Vec2 y((*points)[2*j], (*points)[2*j+1]);
	Vec2 r = x-y;  

	Vec2 n((*normals)[2*j], (*normals)[2*j+1]);
	double potential = -(*weights)[j]*scale*(r.dot(n))/(r.dot(r));
	
	return potential;
}


double ie_Mat::get(unsigned int i, unsigned int j){
	if(is_dynamic){
		if(is_stokes){
			return stokes_kernel(i,j);
		}else{
			return laplace_kernel(i,j);
		}
	}
	assert(i < height_ && j < width_ && mat !=NULL);
	
	return mat[i + lda_*j];
}


void ie_Mat::set(unsigned int i, unsigned int j, double a){
	
	assert(!is_dynamic && i < height_ && j < width_ && mat !=NULL);
	mat[i + lda_*j] = a;
}


void ie_Mat::addset(unsigned int i, unsigned int j, double a){
	assert(!is_dynamic);
	assert(i < height_ && j < width_ && mat !=NULL);
	mat[i + lda_*j] += a;
}


void ie_Mat::set_submatrix(std::vector<unsigned int> I_, 
	std::vector<unsigned int> J_, ie_Mat& A){
	assert(I_.size() == A.height_ && J_.size() == A.width_ && !is_dynamic);
	for(unsigned int i = 0; i < I_.size(); i++){
		for(unsigned int j = 0; j < J_.size(); j++){
			set(I_[i], J_[j], A.get(i,j));
		}
	}
}


void ie_Mat::inverse(){
	assert(mat!= NULL && height_>0 && height_==width_);
	if(!mat) return;


	lapack_int ipiv[height_*2];
	//for some reason, the pivots are written into the array every OTHER element, even with the
	// correct data type (lapack_int) used. To prevent buffer overflow, we double the size of the
	// pivot array
	// TODO investigate this further. 

	LAPACKE_dgetrf(LAPACK_COL_MAJOR, height_, width_, mat, lda_, ipiv );
	LAPACKE_dgetri(LAPACK_COL_MAJOR, height_, mat, lda_, ipiv);
}


ie_Mat& ie_Mat::operator-=(const ie_Mat& o){
	assert(o.height_ == height_ && o.width_ == width_ && !is_dynamic);

	for(unsigned int i=0; i<height_; i++){
		for(unsigned int j=0; j<width_; j++){
			mat[i + lda_*j] =  mat[i + lda_*j] - o. mat[i + lda_*j];
		}
	}
	return *this;
}


ie_Mat& ie_Mat::operator+=(const ie_Mat& o){

	assert(o.height_== height_ && o.width_== width_ && !is_dynamic);
	for(unsigned int i = 0; i < height_; i++){
		for(unsigned int j = 0; j < width_; j++){
			 mat[i + lda_*j] =  mat[i + lda_*j] + o. mat[i + lda_*j];
		}
	}
	return *this;
}


ie_Mat& ie_Mat::operator*=(const double o){
	assert(!is_dynamic);
	for(unsigned int i = 0; i < height_; i++){
		for(unsigned int j = 0; j < width_; j++){
			 mat[i + lda_*j] =  mat[i + lda_*j] *o;
		}
	}
	return *this;
}

//TODO shouldn't this->I have the underscore after it, not this arg?
ie_Mat& ie_Mat::operator()(std::vector<unsigned int> I_, 
	std::vector<unsigned int> J_){

	int olda_ = I_.size();
	double* omat = new double[I_.size()*J_.size()];
	for(unsigned int i = 0; i < I_.size(); i++){
		for(unsigned int j = 0; j < J_.size(); j++){
			omat[i + olda_*j] = get(I_[i],J_[j]);
		}
	}
	ie_Mat* ret  = new ie_Mat();
	ret->mat     = omat;
	ret->height_ = I_.size();
	ret->width_  = J_.size();
	ret->lda_    = ret->height_;
	return *ret;
}


void ie_Mat::operator=(const ie_Mat& copy){
	assert(!is_dynamic);
	
	if(height_ != copy.height_ || width_ != copy.width_ || lda_ != copy.lda_){
		if(mat) delete[] mat;
		lda_ = copy.lda_;
		height_ = copy.height_;
		width_ = copy.width_;
		mat = new double[height_*width_];
	}
	memcpy(mat, copy.mat, width_*height_*sizeof(double));
}

// TODO fix this mess
void ie_Mat::print(){
	// for(int i=0; i<height_; i++){
	// 	for(int j=0; j<width_; j++){
	// 		sprintf("%.7f ",get(i,j));
	// 	}sprintf("\n");
	// }
}

} // namespace ie_solver