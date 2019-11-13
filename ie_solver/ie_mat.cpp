// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include "ie_solver/ie_mat.h"
#include "ie_solver/log.h"

namespace ie_solver {

double ie_Mat::kernel_time = 0.;
double ie_Mat::proxy_time = 0.;
// TODO(John) honestly the lda parameter is completely unused, we should get
// rid of it
ie_Mat::ie_Mat() {
  mat       = NULL;
  lda_      = 0;
  height_   = 0;
  width_    = 0;
}


ie_Mat::~ie_Mat() {
  if (mat) delete[] mat;
}

// Copy constructor
ie_Mat::ie_Mat(const ie_Mat &copy) {
  lda_ = copy.lda_;
  height_ = copy.height_;
  width_ = copy.width_;
  mat = new double[height_ * width_];
  memcpy(mat, copy.mat, width_ * height_ * sizeof(double));
}

// Copy assignment
ie_Mat& ie_Mat::operator=(const ie_Mat& copy) {
  if (height_ != copy.height_ || width_ != copy.width_ || lda_ != copy.lda_) {
    if (mat) delete[] mat;
    lda_ = copy.lda_;
    height_ = copy.height_;
    width_ = copy.width_;
    mat = new double[height_ * width_];
  }
  memcpy(mat, copy.mat, width_ * height_ * sizeof(double));
  return *this;
}


// Move constructor
ie_Mat::ie_Mat(ie_Mat&& move) {
  lda_ = move.lda_;
  height_ = move.height_;
  width_ = move.width_;
  mat = move.mat;
  move.mat = nullptr;
}


// Move assignment
ie_Mat& ie_Mat::operator= (ie_Mat&& move) {
  if (mat) delete[] mat;
  if (height_ != move.height_ || width_ != move.width_ || lda_ != move.lda_) {
    lda_ = move.lda_;
    height_ = move.height_;
    width_ = move.width_;
  }
  mat = move.mat;
  move.mat = nullptr;
  return *this;
}


ie_Mat::ie_Mat(unsigned int h, unsigned int w) {
  lda_      = h;
  height_   = h;
  width_    = w;
  mat       = new double[height_ * width_];
  memset(mat, 0, height_ * width_ * sizeof(double));
}


double ie_Mat::get(unsigned int i, unsigned int j) const {
  assert(i < height_ && j < width_ && mat != NULL);

  return mat[i + lda_ * j];
}


void ie_Mat::set(unsigned int i, unsigned int j, double a) {
  assert(i < height_ && j < width_ && mat != NULL);
  mat[i + lda_ * j] = a;
}


void ie_Mat::addset(unsigned int i, unsigned int j, double a) {
  assert(i < height_ && j < width_ && mat != NULL);
  mat[i + lda_ * j] += a;
}


void ie_Mat::set_submatrix(const std::vector<unsigned int>& I_,
                           const std::vector<unsigned int>& J_,
                           const ie_Mat& A, bool transpose_A) {
  if (transpose_A) {
    assert(I_.size() == A.width_ && J_.size() == A.height_);
    for (unsigned int i = 0; i < I_.size(); i++) {
      for (unsigned int j = 0; j < J_.size(); j++) {
        set(I_[i], J_[j], A.get(j, i));
      }
    }
  } else {
    assert(I_.size() == A.height_ && J_.size() == A.width_);
    for (unsigned int i = 0; i < I_.size(); i++) {
      for (unsigned int j = 0; j < J_.size(); j++) {
        set(I_[i], J_[j], A.get(i, j));
      }
    }
  }
}


void ie_Mat::set_submatrix(unsigned int row_s, unsigned int row_e,
                           unsigned int col_s, unsigned int col_e,
                           const ie_Mat& A, bool transpose_A, bool timing) {
  if (transpose_A) {
    for (unsigned int i = 0; i < row_e - row_s; i++) {
      for (unsigned int j = 0; j < col_e - col_s; j++) {
        set(i + row_s, j + col_s, A.get(j, i));
      }
    }
    assert(row_e - row_s == A.width_ && col_e - col_s == A.height_);
  } else {
    assert(row_e - row_s == A.height_ && col_e - col_s == A.width_);
    for (unsigned int j = 0; j < col_e - col_s; j++) {
      memcpy(&(mat[row_s + lda_ * (j + col_s)]), &(A.mat[A.lda_ * j]),
             (row_e - row_s)*sizeof(double));
    }
  }
}


void ie_Mat::set_submatrix(const std::vector<unsigned int>& I_,
                           unsigned int col_s, unsigned int col_e,
                           const ie_Mat& A, bool transpose_A) {
  if (transpose_A) {
    assert(I_.size() == A.width_ &&  col_e - col_s  == A.height_);
    for (unsigned int i = 0; i < I_.size(); i++) {
      for (unsigned int j = 0; j < col_e - col_s; j++) {
        set(I_[i], j + col_s, A.get(j, i));
      }
    }
  } else {
    assert(I_.size() == A.height_ &&  col_e - col_s  == A.width_);
    for (unsigned int i = 0; i < I_.size(); i++) {
      for (unsigned int j = 0; j < col_e - col_s; j++) {
        set(I_[i], j + col_s, A.get(i, j));
      }
    }
  }
}


void ie_Mat::set_submatrix(unsigned int row_s, unsigned int row_e,
                           const std::vector<unsigned int>& J_,
                           const ie_Mat& A, bool transpose_A) {
  if (transpose_A) {
    assert(row_e - row_s == A.width_ && J_.size() == A.height_);
    for (unsigned int i = 0; i < row_e - row_s; i++) {
      for (unsigned int j = 0; j < J_.size(); j++) {
        set(i + row_s, J_[j], A.get(j, i));
      }
    }
  } else {
    assert(row_e - row_s == A.height_ && J_.size() == A.width_);
    for (unsigned int j = 0; j < J_.size(); j++) {
      memcpy(&(mat[row_s + lda_ * J_[j] ]), &(A.mat[A.lda_ * j]),
             (row_e - row_s)*sizeof(double));
    }
  }
}


ie_Mat ie_Mat::operator()(unsigned int row_s, unsigned int row_e,
                          unsigned int col_s, unsigned int col_e) const {
  ie_Mat submatrix(row_e - row_s, col_e - col_s);
  for (unsigned int i = 0; i < row_e - row_s; i++) {
    for (unsigned int j = 0; j < col_e - col_s; j++) {
      submatrix.set(i, j, this->get(i + row_s, j + col_s));
    }
  }
  return submatrix;
}


void ie_Mat::transpose_into(ie_Mat* transpose) const {
  if (height_ != transpose->width_ || width_ != transpose->height_) {
    if (transpose->mat) delete[] transpose->mat;
    transpose->lda_    = width_;
    transpose->height_ = width_;
    transpose->width_  = height_;
    transpose->mat     = new double[height_ * width_];
  }
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      transpose->mat[j + i * width_] = mat[i + lda_ * j];
    }
  }
}


void ie_Mat::eye(unsigned int n) {
  if (width_ != n || height_ != n || lda_ != n) {
    if (mat) delete[] mat;
    lda_    = n;
    height_ = n;
    width_  = n;
    mat     = new double[height_ * width_];
  }
  for (unsigned int i = 0; i < n; i++) {
    set(i, i, 1.0);
  }
}


unsigned int ie_Mat::height() const {
  return height_;
}


unsigned int ie_Mat::width() const {
  return width_;
}


ie_Mat& ie_Mat::operator-=(const ie_Mat& o) {
  assert(o.height_ == height_ && o.width_ == width_);

  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] - o. mat[i + lda_ * j];
    }
  }
  return *this;
}


ie_Mat& ie_Mat::operator+=(const ie_Mat& o) {
  assert(o.height_ == height_ && o.width_ == width_);
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] + o.mat[i + lda_ * j];
    }
  }
  return *this;
}


ie_Mat& ie_Mat::operator*=(double o) {
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] * o;
    }
  }
  return *this;
}

ie_Mat& ie_Mat::operator/=(double o) {
  assert(std::abs(o) > 1e-8 && "Error: divide matrix by 0.");
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] / o;
    }
  }
  return *this;
}


ie_Mat ie_Mat::operator-() const {
  ie_Mat result(height_, width_);
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      result.set(i, j, -this->get(i, j));
    }
  }
  return result;
}


ie_Mat ie_Mat::operator-(const ie_Mat& o) const {
  assert(o.height_ == height_ && o.width_ == width_);
  ie_Mat result(height_, width_);
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      result.set(i, j, this->get(i, j) - o.get(i, j));
    }
  }
  return result;
}


ie_Mat ie_Mat::operator+(const ie_Mat& o) const {
  assert(o.height_ == height_ && o.width_ == width_);
  ie_Mat sum(height_, width_);
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      sum.set(i, j, this->get(i, j) + o.get(i, j));
    }
  }
  return sum;
}


ie_Mat ie_Mat::operator*(double o) const {
  ie_Mat result(height_, width_);
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      result.set(i, j, this->get(i, j) *o);
    }
  }
  return result;
}


ie_Mat ie_Mat::operator/(double o) const {
  assert(std::abs(o) > 1e-8 && "Error: divide matrix by 0.");
  ie_Mat result(height_, width_);
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      result.set(i, j, this->get(i, j) / o);
    }
  }
  return result;
}


// TODO(John) shouldn't this->I have the underscore after it, not this arg?
ie_Mat ie_Mat::operator()(const std::vector<unsigned int>& I_,
                          const std::vector<unsigned int>& J_) const {
  ie_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (unsigned int i = 0; i < I_.size(); i++) {
    for (unsigned int j = 0; j < J_.size(); j++) {
      assert(I_[i] < height() && J_[j] < width());
      ret.mat[i + olda_ * j] = get(I_[i], J_[j]);
    }
  }
  return ret;
}


ie_Mat ie_Mat::operator()(const std::vector<unsigned int>& I_,
                          unsigned int col_s, unsigned int col_e) const {
  ie_Mat ret(I_.size(), col_e - col_s);
  int olda_ = I_.size();
  for (unsigned int i = 0; i < I_.size(); i++) {
    for (unsigned int j = 0; j < col_e - col_s; j++) {
      assert(I_[i] < height() && col_s + j < width());
      ret.mat[i + olda_ * j] = get(I_[i], col_s + j);
    }
  }
  return ret;
}


ie_Mat ie_Mat::operator()(unsigned int row_s, unsigned int row_e,
                          const std::vector<unsigned int>& J_) const {
  ie_Mat ret(row_e - row_s, J_.size());
  int olda_ = row_e - row_s;
  for (unsigned int i = 0; i < row_e - row_s; i++) {
    for (unsigned int j = 0; j < J_.size(); j++) {
      assert(row_s + i < height() && J_[j] < width());
      ret.mat[i + olda_ * j] = get(row_s + i, J_[j]);
    }
  }
  return ret;
}


double ie_Mat::frob_norm() const {
  double sum = 0;
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      sum += pow(get(i, j), 2);
    }
  }
  return sqrt(sum);
}

double ie_Mat::max_entry_magnitude() const {
  double max = 0;
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      max = std::max(max, std::abs(get(i, j)));
    }
  }
  return max;
}


double ie_Mat::one_norm() const {
  double top = 0;
  for (unsigned int j = 0; j < width_; j++) {
    double sum = 0;
    for (unsigned int i = 0; i < height_; i++) {
      sum += fabs(get(i, j));
    }
    if (sum > top) {
      top = sum;
    }
  }
  return top;
}


double ie_Mat::vec_two_norm() const {
  assert(width() == 1);
  return sqrt(frob_norm());
}


void ie_Mat::rand_vec(unsigned int dofs) {
  // check if we need a resize
  if (width_ != 1 || height_ != dofs || lda_ != dofs) {
    if (mat) delete[] mat;
    lda_    = dofs;
    height_ = dofs;
    width_  = 1;
    mat     = new double[height_ * width_];
  }
  for (unsigned int i = 0; i < dofs; i++) {
    mat[i] = rand() / (0.0 + RAND_MAX);
  }
}

void ie_Mat::left_multiply_pseudoinverse(const ie_Mat& K, ie_Mat* U) const {
  ie_Mat U_(height(), width()), V(height(), width());

  ie_Mat cpy = *this;
  std::vector<double> superb(height());
  std::vector<double> sing(height());
  lapack_int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A',
                                   height(), width(), cpy.mat,
                                   lda_, &(sing[0]), U_.mat,
                                   height(), V.mat, width(),
                                   &(superb[0]));
  assert(info == 0);

  ie_Mat UT_K(height(), K.width());
  ie_Mat::gemm(TRANSPOSE, NORMAL, 1., U_, K, 0., &UT_K);

  for (int row = 0; row < UT_K.height(); row++) {
    double sing_val;
    if (sing[row] > 1e-8) {
      sing_val = 1.0 / sing[row];
    } else {
      sing_val = 0.0;
    }
    for (int col = 0; col < UT_K.width(); col++) {
      UT_K.set(row, col, sing_val * UT_K.get(row, col));
    }
  }
  ie_Mat::gemm(TRANSPOSE, NORMAL, 1., V, UT_K, 0., U);
}


void ie_Mat::LU_factorize(ie_Mat* K_LU, std::vector<lapack_int>* piv) const {
  *K_LU = *this;
  *piv = std::vector<lapack_int>(height_);

  LAPACKE_dgetrf(LAPACK_COL_MAJOR, K_LU->height_, K_LU->width_, K_LU->mat,
                 K_LU->lda_, &(*piv)[0]);
}

void ie_Mat::left_multiply_inverse(const ie_Mat& K, ie_Mat* U) const {
  // X^-1K = U
  // aka, XU = K

  ie_Mat X_copy;
  *U = K;
  std::vector<lapack_int> ipiv;
  LU_factorize(&X_copy, &ipiv);

  int status = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , X_copy.height_ ,
                              U->width_ , X_copy.mat , X_copy.lda_ ,
                              &ipiv[0] , U->mat, U->lda_);
  assert(status == 0);
}


void ie_Mat::right_multiply_inverse(const ie_Mat& K, ie_Mat* L) const {
  ie_Mat K_copy(K.width_, K.height_);
  K.transpose_into(&K_copy);
  std::vector<lapack_int> ipiv;

  // KX^-1 = L
  // aka X_T L^T = K^T
  ie_Mat X_copy;
  LU_factorize(&X_copy, &ipiv);

  int err2 = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'T', X_copy.height_,
                            K_copy.width_, X_copy.mat, X_copy.lda_, &ipiv[0],
                            K_copy.mat, K_copy.lda_);
  assert(err2 == 0);
  K_copy.transpose_into(L);
}


void ie_Mat::left_multiply_inverse(const ie_Mat& K,
                                   const std::vector<lapack_int>& piv,
                                   ie_Mat* U) const {
  // X^-1K = U
  // aka, XU = K

  ie_Mat X_copy = *this;
  *U = K;

  int status = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , this->height_ ,
                              U->width_ , this->mat , this->lda_ ,
                              &piv[0] , U->mat, U->lda_);
  assert(status == 0);
}


void ie_Mat::right_multiply_inverse(const ie_Mat& K,
                                    const std::vector<lapack_int>& piv,
                                    ie_Mat* L) const {
  ie_Mat K_copy(K.width_, K.height_);
  K.transpose_into(&K_copy);

  // KX^-1 = L
  // aka X_T L^T = K^T

  int err2 = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'T', this->height_,
                            K_copy.width_, this->mat, this->lda_, &piv[0],
                            K_copy.mat, K_copy.lda_);
  assert(err2 == 0);
  K_copy.transpose_into(L);
}

// TODO(John) considering LAPACK has transpose options, should this ever be
// called?
ie_Mat ie_Mat::transpose() const {
  ie_Mat transpose(width(), height());
  transpose_into(&transpose);
  return transpose;
}


double ie_Mat::condition_number() const {
  // TODO(John) this should be done without using LU factorization if possible?
  ie_Mat X_copy = *this;
  std::vector<lapack_int> ipiv(height_);
  memset(&ipiv[0], 0, height_ * sizeof(lapack_int));
  double one_norm =  X_copy.one_norm();
  LAPACKE_dgetrf(LAPACK_COL_MAJOR, X_copy.height_, X_copy.width_, X_copy.mat,
                 X_copy.lda_, &ipiv[0]);
  double rcond;
  LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', X_copy.height(), X_copy.mat,
                 X_copy.lda_, one_norm, &rcond);
  return 1.0 / rcond;
}


// Performs interpolative decomposition, and returns number of skeleton columns.
// Takes double /tol/, tolerance factorfor error in CPQR factorization.
// Populates /p/ with permutation, Z with linear transformation.
int ie_Mat::id(std::vector<unsigned int>* p, ie_Mat* Z, double tol) const {
  ie_Mat cpy = *this;
  assert(height() != 0 &&  width() != 0);
  std::vector<lapack_int> pvt(width_);
  memset(&pvt[0], 0, width_ * sizeof(lapack_int));

  // /tau/ will contain an output from dgeqp3 that we don't need.
  std::vector<double> tau(width_);
  int info1 = LAPACKE_dgeqp3(CblasColMajor, height_, width_, cpy.mat, lda_,
                             &pvt[0], &tau[0]);
  assert(info1 == 0);
  unsigned int skel = 0;

  double thresh = fabs(tol * cpy.get(0, 0));
  for (unsigned int i = 1; i < width_; i++) {
    // check if R_{i,i} / R_{0,0} < tol
    if (fabs(cpy.get(i, i)) < thresh) {
      skel = i;
      break;
    }
  }
  if (skel == 0) {
    // no compression to be done :/
    return 0;
  }
  for (unsigned int i = 0; i < width_; i++) {
    p->push_back(pvt[i] - 1);
  }

  int redund = width_ - skel;

  // set Z to be R_11^-1 R_12. Note 'U' (above diagonal) part of cp.mat
  // is the R matrix from dgeqp3.
  int info2 = LAPACKE_dtrtrs(CblasColMajor, 'U', 'N', 'N', skel, redund,
                             cpy.mat, cpy.lda_, cpy.mat + cpy.lda_ * skel,
                             cpy.lda_);

  assert(info2 == 0);
  *Z = cpy(0, skel, skel, skel + redund);
  return skel;
}


std::vector<double> ie_Mat::real_eigenvalues() {
  std::vector<double> eigvs;
  double *eigs = new double[width_];
  double *imags = new double[width_];
  int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', width_, mat,
                           width_, eigs, imags,
                           nullptr, 1, nullptr, 1);
  assert(info == 0);
  for (unsigned int i = 0; i < width_; i++) {
    if (fabs(imags[i]) < 1e-14) {
      eigvs.push_back(eigs[i]);
    }
  }

  delete[] imags;
  delete[] eigs;
  return eigvs;
}


void ie_Mat::print() const {
  std::string message = "\n(Printing matrix*1000)\n";
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      message += std::to_string(1000 * get(i, j)) + " ";
    }
    message += "\n";
  }
  LOG::INFO(message);
}


void ie_Mat::write_singular_values_to_file(const std::string& filename) const {
  std::ofstream output;
  output.open(filename);

  ie_Mat cpy = *this;
  std::vector<double> superb(height());
  std::vector<double> sing(height());
  ie_Mat U(height(), height());
  lapack_int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'N', 'N',
                                   height(), width(), cpy.mat,
                                   lda_, &(sing[0]), U.mat,
                                   height(), nullptr, width(),
                                   &(superb[0]));
  assert(info == 0);
  if (output.is_open()) {
    for (unsigned int i = 0; i < sing.size(); i++) {
      output << sing[i] << std::endl;
    }
    output.close();
  } else {
    printf("Failed to open singular values output file!\n");
  }
}


void ie_Mat::gemv(CBLAS_TRANSPOSE trans0, double alpha, const ie_Mat& A,
                  const ie_Mat& x, double beta, ie_Mat* b) {
  assert(A.height() != 0 && x.height() != 0 && A.width() != 0  &&
         "gemv needs positive dimensions only.");
  assert(x.width() == 1 && "gemv only works on a column vector.");
  assert(A.mat != nullptr && "gemv fails on null A mat.");
  assert(x.mat != nullptr && "gemv fails on null x mat.");
  assert(b->mat != nullptr && "gemv fails on null b mat.");


  cblas_dgemv(CblasColMajor, trans0, A.height(), A.width(), alpha, A.mat,
              A.lda_, x.mat, 1, beta, b->mat, 1);
}


void ie_Mat::gemm(CBLAS_TRANSPOSE trans0, CBLAS_TRANSPOSE trans1,
                  double alpha, const ie_Mat& A, const ie_Mat& B, double beta,
                  ie_Mat* C) {
  assert(A.height() != 0 && B.height() != 0 && A.width() != 0
         && B.width() != 0 && "gemm needs positive dimensions only.");
  assert(A.mat != nullptr && "gemm fails on null A mat.");
  assert(B.mat != nullptr && "gemm fails on null B mat.");
  assert(C->mat != nullptr && "gemm fails on null C mat.");

  unsigned int k = (trans0 == CblasTrans) ? A.height() : A.width();
  cblas_dgemm(CblasColMajor, trans0, trans1, C->height(), C->width(),
              k, alpha, A.mat, A.height(), B.mat, B.height(), beta, C->mat,
              C->height());
}

}  // namespace ie_solver
