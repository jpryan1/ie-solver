#ifndef _VEC2_H_
#define _VEC2_H_

#include "log.h"

namespace ie_solver{

struct Vec2{
	double a[2];
	Vec2(double* arr){
		a[0] = arr[0];
		a[1] = arr[1];
	}
	Vec2(){
		a[0] = 0;
		a[1] = 0;
	}
	Vec2(double m, double n){
		a[0] = m;
		a[1] = n;
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
		const std::string message = "Vec2: " + std::to_string(a[0]) + " " + 
			std::to_string(a[1]);
		LOG::INFO(message);
	}
};

} // namespace

#endif