#include "common.h"

namespace ie_solver{


ie_Mat::ie_Mat(){
	is_dynamic   = false;
	is_stokes = false;
	mat       = NULL;
	lda_      = 0;
	height_   = 0;
	width_    = 0;
}


ie_Mat::~ie_Mat(){
	if(mat) delete[] mat;
}

ie_Mat::ie_Mat(unsigned int h, unsigned int w){
	is_dynamic   = false;
	is_stokes = false;
	lda_      = h;
	height_   = h;
	width_    = w;
	mat       = new double[height_*width_];
	memset(mat, 0, height_*width_*sizeof(double));
}


ie_Mat& ie_Mat::operator=(const ie_Mat& copy){
	assert(!is_dynamic);
	
	if(height_ != copy.height_ || width_ != copy.width_ || lda_ != copy.lda_){
		if(mat) delete[] mat;
		lda_ = copy.lda_;
		height_ = copy.height_;
		width_ = copy.width_;
		mat = new double[height_*width_];
	}
	memcpy(mat, copy.mat, width_*height_*sizeof(double));
	return *this;
}


void ie_Mat::resize(unsigned int h, unsigned int w){
	if(mat) delete[] mat;
	lda_ = h;
	height_ = h;
	width_ = w;
	mat = new double[height_*width_];
	memset(mat, 0, height_*width_*sizeof(double));
}


void ie_Mat::copy(ie_Mat& copy) const{

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


double ie_Mat::get(unsigned int i, unsigned int j) const{
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


void ie_Mat::set_submatrix(const std::vector<unsigned int>& I_, 
	const std::vector<unsigned int>& J_, ie_Mat& A){
	assert(I_.size() == A.height_ && J_.size() == A.width_ && !is_dynamic);
	for(unsigned int i = 0; i < I_.size(); i++){
		for(unsigned int j = 0; j < J_.size(); j++){
			set(I_[i], J_[j], A.get(i,j));
		}
	}
}


void ie_Mat::transpose(ie_Mat& A) const {
	if(height_ != A.width_ || width_ != A.height_){
		if(A.mat) delete[] A.mat;
		A.lda_    = width_;
		A.height_ = width_;
		A.width_  = height_;
		A.mat     = new double[height_*width_];
	}
	for(unsigned int i = 0; i < height_; i++){
		for(unsigned int j = 0; j < width_; j++){
			A.set(j,i, get(i,j));
		}
	}
}


unsigned int ie_Mat::height() const {
	return height_;
}


unsigned int ie_Mat::width() const {
	return width_;
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


ie_Mat& ie_Mat::operator*=(double o){
	assert(!is_dynamic);
	for(unsigned int i = 0; i < height_; i++){
		for(unsigned int j = 0; j < width_; j++){
			 mat[i + lda_*j] =  mat[i + lda_*j] *o;
		}
	}
	return *this;
}


//TODO shouldn't this->I have the underscore after it, not this arg?
ie_Mat ie_Mat::operator()(const std::vector<unsigned int>& I_, 
	const std::vector<unsigned int>& J_) const {

	ie_Mat ret(I_.size(), J_.size());

	int olda_ = I_.size();
	for(unsigned int i = 0; i < I_.size(); i++){
		for(unsigned int j = 0; j < J_.size(); j++){
			ret.mat[i + olda_*j] = get(I_[i],J_[j]);
		}
	}
	return ret;
}


double ie_Mat::norm2() const{
	double sum=0;
	for(unsigned int i = 0; i < height_; i++){
		for(unsigned int j = 0; j < width_; j++){
			sum += pow( get(i,j), 2);
		}
	}
	return sqrt(sum);
}


// This function stores the DoF data,  and calculates the diagonals of the mat
void ie_Mat::load(const std::vector<double>* p, const std::vector<double>* n, 
	const std::vector<double>* c, const std::vector<double>* w){
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


double ie_Mat::stokes_kernel(unsigned int i, unsigned int j) const {
	// So this is much more awkwardly written than the function in 
	// stokes_init.cpp that just writes the entire matrix at once, but the cost 
	// of writing the whole matrix is just stupid. So we power through this 
	// function
	// below commented is single layer
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
			return -0.5 - 
			0.5*(*curvatures)[dof_i]*(*weights)[dof_i]*scale*t0*t0;
		}
		else if(i%2==1 && j%2==1){
			return -0.5 - 
			0.5*(*curvatures)[dof_i]*(*weights)[dof_i]*scale*t1*t1;
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

double ie_Mat::laplace_kernel(unsigned int i, unsigned int j) const {
	
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


void ie_Mat::left_multiply_inverse(const ie_Mat& K, ie_Mat& U) const {
	// X^-1K = U
	//aka, XU = K

	//TODO insert asserts for these functions

	ie_Mat X_copy(height_, width_);
	copy(X_copy);

	ie_Mat K_copy(K.height_, K.width_);

	K.copy(K_copy);

	lapack_int ipiv[height_];
	memset(ipiv,0, height_*sizeof(lapack_int));
	
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, X_copy.height_, X_copy.width_, X_copy.mat,
	 X_copy.lda_, ipiv );

	int status = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , X_copy.height_ , 
	 K_copy.width_ , X_copy.mat , X_copy.lda_ , ipiv , K_copy.mat, K_copy.lda_);
	
	assert(status==0);
	K_copy.copy(U);
}


void ie_Mat::right_multiply_inverse(const ie_Mat& K, ie_Mat& L) const {
	// KX^-1 = L
	// aka X_T L^T = K^T
	
	ie_Mat X_copy(height_, width_);
	transpose(X_copy);

	ie_Mat K_copy(K.width_, K.height_);

	K.transpose(K_copy);

	lapack_int ipiv[height_];
	memset(ipiv,0, height_*sizeof(lapack_int));
	
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, X_copy.height_, X_copy.width_, X_copy.mat, 
		X_copy.lda_, ipiv );

	int err = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', X_copy.height_, 
		K_copy.width_, X_copy.mat, X_copy.lda_, ipiv, K_copy.mat, K_copy.lda_);
	
	assert(err==0);

	K_copy.transpose(L);
}


// Performs interpolative decomposition, and eturns number of skeleton columns. 
// Takes double /tol/, tolerance factorfor error in CPQR factorization. 
// Populates /p/ with permutation, Z with linear transformation.
int ie_Mat::id( std::vector<unsigned int>& p, ie_Mat& Z, double tol) const {
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
	Z = cpy(I_, J_);	
	return skel;
}


void ie_Mat::print() const{
	std::string message = "\n(Printing matrix*1000)\n";
	for(unsigned int i=0; i<height_; i++){
		for(unsigned int j=0; j<width_; j++){
			message += std::to_string(1000*get(i,j)) + " ";
		}
		message += "\n";
	}
	LOG::INFO(message);
}

} // namespace ie_solver