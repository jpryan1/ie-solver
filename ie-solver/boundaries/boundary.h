#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "ie-solver/vec2.h"
#include "ie-solver/ie_mat.h"
#include <vector> 

namespace ie_solver{

class Boundary {
	public:
		std::vector<double> points, normals, curvatures, weights;
		ie_Mat boundary_condition;
		virtual void initialize(int n, int bc_enum) = 0;
		virtual bool is_in_domain(Vec2& a) = 0;

};

} // namespace ie_solver

#endif