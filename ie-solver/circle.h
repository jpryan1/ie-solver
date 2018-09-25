#ifndef _CIRCLE_H_
#define _CIRCLE_H_


#include "common.h"

namespace ie_solver{

void circle(int, std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights);

int out_of_circle(Vec2& a);

} // namespace ie_solver

#endif