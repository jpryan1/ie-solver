#ifndef _SQUIGGLY_H_
#define _SQUIGGLY_H_
#include <boost/math/special_functions/ellint_2.hpp>
#include "common.h"


void squiggly(int, std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights);

int out_of_squiggly(Vec2& a);

#endif
