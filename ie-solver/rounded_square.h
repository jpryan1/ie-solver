#ifndef _ROUNDEDSQUARE_H_
#define _ROUNDEDSQUARE_H_

#include "common.h"


void rounded_square(std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights);

int out_of_rounded_square(Vec2& a);

#endif
