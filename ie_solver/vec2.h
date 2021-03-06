// Copyright 2019 John Paul Ryan
#ifndef IE_SOLVER_VEC2_H_
#define IE_SOLVER_VEC2_H_

namespace ie_solver {

struct Vec2 {
  double a[2];

  Vec2();
  Vec2(double m, double n);
  double norm() const;
  double dot(const Vec2& o) const;
  Vec2 operator-(const Vec2 &o) const;
  Vec2 operator*(const double d) const;
  void print() const;
};

}  // namespace ie_solver

#endif  // IE_SOLVER_VEC2_H_
