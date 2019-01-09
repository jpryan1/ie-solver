#include "initialization.h"
#include <cmath>

namespace ie_solver{


//TODO now points vec might need to be boundary_points instead
void Initialization::InitializeDomainKernel(ie_Mat& K, std::vector<double>& points,
										std::vector<double>& normals, 
										std::vector<double>& weights, 
										std::vector<double>& domain_points, int test_size,
										Boundary* boundary,
										bool is_stokes){
	// is stokes TODO

	if(is_stokes){
		Stokes_InitializeDomainKernel(K, points, normals, weights, 
			domain_points, test_size, boundary);
		return;
	}

	//columns for phi (aka dofs), rows for spatial domain
	int dofs = points.size()/2;
	double scale = 1.0 / (2*M_PI);
	// omp_set_num_threads(4);
	// #pragma omp parallel for 	
	for(int i = 0; i < test_size*test_size; i++){

		Vec2 x(domain_points[2*i], domain_points[2*i+1]);
		for(int j = 0; j < dofs; j++){
			
			Vec2 y(points[2*j], points[2*j+1]);

			if(!boundary->is_in_domain(x)){
				K.set(i, j, 0);
				continue;
			}
			Vec2 r = x-y;
			Vec2 n(normals[2*j], normals[2*j+1]);
			double potential = -weights[j]*scale*(r.dot(n))/(r.dot(r));
			K.set(i, j, potential);
		}
	}
}


void Initialization::DomainSolution(ie_Mat& K, int test_size, 
					double min, double max,int (*out_of_domain)(Vec2& a)){
	//columns for phi (aka dofs), rows for spatial domain
	// #pragma omp parallel for 	
	for(int i = 0; i < test_size*test_size; i++){
		double x0 = i/test_size;
		x0 = min + (x0*(max-min))/test_size;
		double x1 = i%test_size;
		x1 = min + (x1*(max-min))/test_size;

		Vec2 v(x0, x1);
		if(out_of_domain(v)){
			K.set(i, 0, 0);
			continue;
		}
		double potential = log(sqrt( pow(x0+2,2)+pow(x1+2,2)))/(2*M_PI);	
		//double potential = -weights[ind_j]*scale*(r.dot(n))/(r.dot(r));
		K.set(i, 0, potential);
	}
}


void Initialization::Electric_InitializeBoundary(ie_Mat& f, std::vector<double>& points){
	for(unsigned int i = 0; i < f.height(); i++){
		f.set(i, 0, 1);	
	}
}

} // namespace ie_solver