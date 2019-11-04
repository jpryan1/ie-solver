// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_IE_MAT_H_
#define IE_SOLVER_IE_MAT_H_

#include <cblas.h>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#include <vector>
#include <string>

#define NORMAL CblasNoTrans
#define TRANSPOSE CblasTrans

namespace ie_solver {


struct ie_Mat {
  // storage is column major, so by default lda is the height.
  double *mat;
  static double proxy_time, kernel_time;

  unsigned int lda_, height_, width_;
  ie_Mat();
  ~ie_Mat();

  ie_Mat(unsigned int h, unsigned int w);

  // Copy constructor
  ie_Mat(const ie_Mat &o);
  // Copy assignment
  ie_Mat& operator=(const ie_Mat& copy);
  // Move constructor
  ie_Mat(ie_Mat&& move);
  // Move assignment
  ie_Mat& operator=(ie_Mat&& move);

  double get(unsigned int i, unsigned int j) const;
  void set(unsigned int i, unsigned int j, double a);
  void addset(unsigned int i, unsigned int j, double a);
  void set_submatrix(const std::vector<unsigned int>& I_,
                     const std::vector<unsigned int>& J_, const ie_Mat& A,
                     bool transpose_A = false);
  void set_submatrix(unsigned int row_s, unsigned int row_e, unsigned int col_s,
                     unsigned int col_e,
                     const ie_Mat& A, bool transpose_A = false,
                     bool timing = false);
  void set_submatrix(const std::vector<unsigned int>& I_, unsigned int col_s,
                     unsigned int col_e,
                     const ie_Mat& A, bool transpose_A = false);
  void set_submatrix(unsigned int row_s, unsigned int row_e,
                     const std::vector<unsigned int>& J_,
                     const ie_Mat& A, bool transpose_A = false);

  void transpose_into(ie_Mat* transpose) const;
  void eye(unsigned int n);
  ie_Mat transpose() const;

  unsigned int height() const;
  unsigned int width() const;

  ie_Mat& operator-=(const ie_Mat& o);
  ie_Mat& operator+=(const ie_Mat& o);
  ie_Mat& operator*=(double o);
  ie_Mat& operator/=(double o);

  ie_Mat operator-() const;
  ie_Mat operator-(const ie_Mat& o) const;
  ie_Mat operator+(const ie_Mat& o) const;
  ie_Mat operator*(double o) const;
  ie_Mat operator/(double o) const;

  ie_Mat operator()(const std::vector<unsigned int>& I_,
                    const std::vector<unsigned int>& J_) const;
  ie_Mat operator()(const std::vector<unsigned int>& I_,
                    unsigned int col_s, unsigned int col_e) const;
  ie_Mat operator()(unsigned int row_s, unsigned int row_e,
                    const std::vector<unsigned int>& J_) const;
  ie_Mat operator()(unsigned int row_s, unsigned int row_e, unsigned int col_s,
                    unsigned int col_e) const;


  double one_norm() const;
  double frob_norm() const;
  void write_singular_values_to_file(const std::string& filename) const;


  // This function stores the DoF data, and calculates the diagonals of the
  // mat.
  void rand_vec(unsigned  dofs);
  double condition_number() const;

  void LU_factorize(ie_Mat* K_LU, std::vector<lapack_int>* piv) const;
  void left_multiply_inverse(const ie_Mat& K, ie_Mat* U) const;
  void right_multiply_inverse(const ie_Mat& K, ie_Mat* L) const;
  void left_multiply_inverse(const ie_Mat& K,
                             const std::vector<lapack_int>& piv,
                             ie_Mat* U) const;
  void right_multiply_inverse(const ie_Mat& K,
                              const std::vector<lapack_int>& piv,
                              ie_Mat* L) const;
  void left_multiply_pseudoinverse(const ie_Mat& K, ie_Mat* U) const;

  int id(std::vector<unsigned int>* p, ie_Mat* Z, double tol) const;
  std::vector<double> real_eigenvalues();

  void print() const;

  static void gemv(CBLAS_TRANSPOSE trans0, double alpha, const ie_Mat& A,
                   const ie_Mat& x, double beta, ie_Mat* b);
  static void gemm(CBLAS_TRANSPOSE trans0, CBLAS_TRANSPOSE trans1,
                   double alpha, const ie_Mat& A, const ie_Mat& B, double beta,
                   ie_Mat* C);
};

}  // namespace ie_solver

#endif  // IE_SOLVER_IE_MAT_H_
