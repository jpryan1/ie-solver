#ifndef _ROUNDEDSQUARE_H_
#define _ROUNDEDSQUARE_H_

#include "common.h"

namespace ie_solver{

void rounded_square(int N, std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights);

int out_of_rounded_square(Vec2& a);

} // namespace

#endif