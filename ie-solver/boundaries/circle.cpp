#include "circle.h"
#include <cmath>

namespace ie_solver{

void Circle::initialize(int N){

	for(int i = 0; i < N; i++){
		double ang = i*2.0*M_PI/N;
		points.push_back(0.5 + 0.25*cos(ang));
		points.push_back(0.5 + 0.25*sin(ang));
		normals.push_back(cos(ang));
		normals.push_back(sin(ang));
		curvatures.push_back(4); // 1/r, r=0.25
		weights.push_back(M_PI/(N*2));
	}

}

bool Circle::is_in_domain(Vec2& a){

	double x = a.a[0]-0.5;
	double y = a.a[1]-0.5;
	double eps = 1e-2;

	double dist = sqrt(pow(x,2)+pow(y,2));
	if(dist + eps> 0.25) return false;
	return true;
}

} // namespace ie_solver