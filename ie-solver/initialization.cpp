#include "initialization.h"

namespace ie_solver{

//TODO - subclasses for different PDES
void Initialization::InitializeKernel(ie_Mat& K, std::vector<double>& points,
										std::vector<double>& normals, 
										std::vector<double>& curvatures,
										std::vector<double>& weights){
	int dofs = points.size() / 2;
	double scale = 1.0 / (2*M_PI);

	for(int i = 0; i < dofs; i++){
	
		for(int j = 0; j < dofs; j++){
			
			if(i==j){
				K.set(i, j, 0.5 + 0.5*curvatures[i]*weights[i]*scale);
				continue;
			 }

		 	Vec2 x(points[2*i], points[2*i+1]);
			Vec2 y(points[2*j], points[2*j+1]);
			Vec2 r = x-y;  

			Vec2 n(normals[2*j], normals[2*j+1]);
			double potential = -weights[j]*scale*(r.dot(n))/(r.dot(r));
			K.set(i, j, potential);
		}
	}
}


void Initialization::InitializeDomainKernel(ie_Mat& K, std::vector<double>& points,
										std::vector<double>& normals, 
										std::vector<double>& weights, 
										double min, double max, int test_size,
										int (*out_of_domain)(Vec2& a)){
	//columns for phi (aka dofs), rows for spatial domain
	int dofs = points.size()/2;
	double scale = 1.0 / (2*M_PI);
	omp_set_num_threads(4);
	// #pragma omp parallel for 	
	for(int i = 0; i < test_size*test_size; i++){
		double x0 = i/test_size;
		x0 = min + (x0*(max-min))/test_size;
		double x1 = i%test_size;
		x1 = min + (x1*(max-min))/test_size;
		Vec2 x(x0, x1);
		for(int j = 0; j < dofs; j++){
			
			Vec2 y(points[2*j], points[2*j+1]);

			if(out_of_domain(x)){
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
	omp_set_num_threads(4);
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


void Initialization::InitializeBoundary(ie_Mat& f, std::vector<double>& points){
	for(int i = 0; i < f.height(); i++){
		double x0 = points[2*i]+2;
		double y0 = points[2*i+1]+2;

		double potential = log(sqrt( pow(x0,2)+pow(y0,2)))/(2*M_PI);
		f.set(i, 0, potential);
	}
}

} // namespace ie_solver