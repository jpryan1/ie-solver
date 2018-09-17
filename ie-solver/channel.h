#ifndef _CHANNEL_H_
#define _CHANNEL_H_



#include <boost/math/special_functions/ellint_2.hpp>
#include "common.h"


void channel(std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights);

int out_of_channel(Vec2& a);

#endif
