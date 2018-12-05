#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_

#include "common.h"
#include "quadtree.h"

namespace ie_solver{
	
struct Initialization{

	Initialization(){}
	~Initialization(){}
	
	void InitializeKernel(ie_Mat& K, std::vector<double>&, std::vector<double>&, std::vector<double>&,
		std::vector<double>&, bool);
	void InitializeDomainKernel(ie_Mat& K, std::vector<double>&, std::vector<double>&,
		std::vector<double>&, std::vector<double>&, int, int (*out_of_domain)(Vec2& a), bool);

	void InitializeBoundary(ie_Mat& f, std::vector<double>& points);
	void DomainSolution(ie_Mat& K, int, double, double, int (*out_of_domain)(Vec2& a));

	void Stokes_InitializeKernel(ie_Mat& K, std::vector<double>&, std::vector<double>&, std::vector<double>&,
		std::vector<double>&);
	void Stokes_InitializeDomainKernel(ie_Mat& K, std::vector<double>&, std::vector<double>&,
		std::vector<double>&, std::vector<double>&, int, int (*out_of_domain)(Vec2& a));

	void Stokes_InitializeBoundary(ie_Mat& f, std::vector<double>& points);
	
};

} // namespace ie_solver

#endif
