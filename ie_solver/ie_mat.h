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

  int lda_, height_, width_;
  ie_Mat();
  ~ie_Mat();

  ie_Mat(int h, int w);

  // Copy constructor
  ie_Mat(const ie_Mat &o);
  // Copy assignment
  ie_Mat& operator=(const ie_Mat& copy);
  // Move constructor
  ie_Mat(ie_Mat&& move);
  // Move assignment
  ie_Mat& operator=(ie_Mat&& move);

  double get(int i, int j) const;
  void set(int i, int j, double a);
  void addset(int i, int j, double a);
  void set_submatrix(const std::vector<int>& I_,
                     const std::vector<int>& J_, const ie_Mat& A,
                     bool transpose_A = false);
  void set_submatrix(int row_s, int row_e, int col_s,
                     int col_e,
                     const ie_Mat& A, bool transpose_A = false,
                     bool timing = false);
  void set_submatrix(const std::vector<int>& I_, int col_s,
                     int col_e,
                     const ie_Mat& A, bool transpose_A = false);
  void set_submatrix(int row_s, int row_e,
                     const std::vector<int>& J_,
                     const ie_Mat& A, bool transpose_A = false);

  void transpose_into(ie_Mat* transpose) const;
  void eye(int n);
  ie_Mat transpose() const;

  int height() const;
  int width() const;

  ie_Mat& operator-=(const ie_Mat& o);
  ie_Mat& operator+=(const ie_Mat& o);
  ie_Mat& operator*=(double o);
  ie_Mat& operator/=(double o);

  ie_Mat operator-() const;
  ie_Mat operator-(const ie_Mat& o) const;
  ie_Mat operator+(const ie_Mat& o) const;
  ie_Mat operator*(const ie_Mat& o) const;
  ie_Mat operator*(double o) const;
  ie_Mat operator/(double o) const;

  ie_Mat operator()(const std::vector<int>& I_,
                    const std::vector<int>& J_) const;
  ie_Mat operator()(const std::vector<int>& I_,
                    int col_s, int col_e) const;
  ie_Mat operator()(int row_s, int row_e,
                    const std::vector<int>& J_) const;
  ie_Mat operator()(int row_s, int row_e, int col_s,
                    int col_e) const;


  double one_norm() const;
  double vec_two_norm() const;
  double frob_norm() const;
  double max_entry_magnitude() const;
  void write_singular_values_to_file(const std::string& filename) const;


  // This function stores the DoF data, and calculates the diagonals of the
  // mat.
  void rand_vec(int dofs);
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

  int id(std::vector<int>* p, ie_Mat* Z, double tol) const;
  std::vector<double> real_eigenvalues();

  void print() const;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_IE_MAT_H_
