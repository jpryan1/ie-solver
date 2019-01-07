#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "ie-solver/boundaries/boundary.h"
#include "ie_mat.h"

namespace ie_solver{

struct Kernel{


	double scale, diag_00, diag_01, diag_11;
	bool is_stokes;
	// TODO don't have the kernel store the boundary
	Boundary* boundary;
	double electric_kernel(unsigned int i, unsigned int j) const;

	double get(unsigned int i, unsigned int j) const;

	double stokes_kernel(unsigned int i, unsigned int j) const ;
	double laplace_kernel(unsigned int i, unsigned int j) const ;

	// This function stores the DoF data,  and calculates the diagonals of the mat
	void load(Boundary* boundary, bool is_stokes_);
//TODO shouldn't this->I have the underscore after it, not this arg?
	ie_Mat operator()(const std::vector<unsigned int>& I_, 
		const std::vector<unsigned int>& J_) const;

}; // struct

} // namespace ie_solver

#endif