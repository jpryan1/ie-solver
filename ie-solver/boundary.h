#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "ie-solver/vec2.h"
#include <vector> 

namespace ie_solver{

class Boundary {
	public:
		std::vector<double> points, normals, curvatures, weights;
		virtual void initialize(int n) = 0;
		virtual bool is_in_domain(Vec2& a) = 0;

};

} // namespace ie_solver

#endif