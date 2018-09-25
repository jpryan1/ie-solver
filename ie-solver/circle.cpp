#include "circle.h"

namespace ie_solver{

void circle(int N, std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights){



	for(int i = 0; i < N; i++){
		double ang = i*2.0*M_PI/N;
		points.push_back(0.5 + 0.25*cos(ang));
		points.push_back(0.5 + 0.25*sin(ang));
		normals.push_back(cos(ang));
		normals.push_back(sin(ang));
		curvatures.push_back(4);
		weights.push_back(M_PI/(N*2));
	}





	for(int i = 0; i < N; i++){
		double ang = i*2.0*M_PI/N;
		points.push_back(0.5 + 0.10*cos(ang));
		points.push_back(0.5 + 0.10*sin(ang));
		normals.push_back(-cos(ang));
		normals.push_back(-sin(ang));
		curvatures.push_back(10);
		weights.push_back(M_PI/(N*5));
	}

//for exterior
	// for(int i=0; i<normals.size(); i++){
	// 	normals[i] = -1*normals[i];
	// }
}





int out_of_circle(Vec2& a){
	

	double x = a.a[0]-0.5;
	double y = a.a[1]-0.5;
	double eps = 1e-2;

	double dist = sqrt(pow(x,2)+pow(y,2));
	 if(dist + eps> 0.25) return 1;
	if(dist - eps < 0.1) return 1;
	 return 0;


//for exterior
	// if(dist - eps> 0.25) return 0;
	// return 1;

}

} // namespace ie_solver