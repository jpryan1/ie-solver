#include "initialization.h"
#include <cassert>
#include <cmath>

namespace ie_solver{
//TODO - subclasses for different PDES
	

void Initialization::Stokes_InitializeDomainKernel(ie_Mat& K, std::vector<double>& points,
										std::vector<double>& normals, 
										std::vector<double>& weights, 
										std::vector<double>& domain_points, int test_size,
										Boundary* boundary){
	//columns for phi (aka dofs), rows for spatial domain
	
	assert(points.size() == normals.size() && 
		"Points and normals must have same size in Stokes domain init.");
	assert(points.size() == 2*weights.size() &&
		"In dimension 2, pts must be 2*size of weights in stokes domain init.");
	
	// int dofs = points.size()/2;
	double scale = 1.0 / (M_PI);

	// min-=d/2;
	// max+=d/2;
	//#pragma omp parallel for 

	for(int i=0; i<test_size*test_size; i++){
		Vec2 x(domain_points[2*i], domain_points[2*i+1]);
		
		for(unsigned int j=0; j<points.size(); j+=2){
			int ind_j = j/2;
			
			Vec2 y(points[j], points[j+1]);

			if(!boundary->is_in_domain(x)){
				K.set(2*i  , 2*ind_j  , 0);
				K.set(2*i+1, 2*ind_j  , 0);
				K.set(2*i  , 2*ind_j+1, 0);
				K.set(2*i+1, 2*ind_j+1, 0);
				continue;
			}

			Vec2 r = x-y;
			Vec2 n(normals[j], normals[j+1]);
			double r0 = r.a[0];
			double r1 = r.a[1];

			// K.set(2*i  , 2*ind_j  , weights[ind_j]*scale*(log(1.0/r.norm()) + (r0*r0/r.dot(r)) ));
			// K.set(2*i+1, 2*ind_j  , weights[ind_j]*scale*(                    (r0*r1/r.dot(r)) ));
			// K.set(2*i  , 2*ind_j+1, weights[ind_j]*scale*(                    (r1*r0/r.dot(r)) ));
			// K.set(2*i+1, 2*ind_j+1, weights[ind_j]*scale*(log(1.0/r.norm()) + (r1*r1/r.dot(r)) ));
			

		//	DOUBLE LAYER POTENTIAL
			double potential = weights[ind_j]*scale*(r.dot(n))/(pow(r.dot(r),2));
			
			
			K.set(2*i  , 2*ind_j  , potential*r0*r0);
			K.set(2*i+1, 2*ind_j  , potential*r0*r1);
			K.set(2*i  , 2*ind_j+1, potential*r0*r1);
			K.set(2*i+1, 2*ind_j+1, potential*r1*r1);
			
		}
	}


}


void Initialization::Stokes_InitializeBoundary(ie_Mat& f, 
	std::vector<double>& normals){
	assert(f.height()==normals.size());
	for(unsigned int i=0; i<f.height()/2; i+=2){
		f.set(i,   0, 1);
		f.set(i+1, 0, 0);
	}
	for(unsigned int i=f.height()/2; i<f.height(); i+=2){
		f.set(i,   0, 0);
		f.set(i+1, 0, 0);
	}
}

} // namespace