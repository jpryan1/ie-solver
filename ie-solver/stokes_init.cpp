#include "initialization.h"

namespace ie_solver{
//TODO - subclasses for different PDES

void Initialization::Stokes_InitializeKernel(ie_Mat& K, std::vector<double>& points,
										std::vector<double>& normals, std::vector<double>& curvatures,
										std::vector<double>& weights){

	double avg = 0;
	for(unsigned int i=0; i<weights.size(); i++) avg += weights[i];
	avg /= weights.size();
	
	// int Q = floor(sqrt(K.height()));
	double alpha = avg/2.0;
	double beta = alpha*alpha;
	double singular_00= 8*beta + 2*beta*log(1/(2*beta)) - beta*M_PI;
	double singular_01=0;

	//singular_00 = 10000;

	double singular_11 = singular_00;

	printf("Singulars calculated as %f %f %f\n", singular_00, singular_01, singular_11);
	// exit(0);


// //SINGLE LAYER BELOW
	//double scale = 1.0 / (4*M_PI);

	// for(int i=0; i<points.size(); i+=2){
	// //	omp_set_num_threads(4);
	// //	#pragma omp parallel for 	
	
	// 	for(int j=0; j<points.size(); j+=2){
			
	// 		int ind_i = i/2;
			
	// 		if(i==j){
	// 			K.set(2*ind_i  , 2*ind_i  , scale*singular_00);
	// 			K.set(2*ind_i+1, 2*ind_i  , scale*singular_01);
	// 			K.set(2*ind_i  , 2*ind_i+1, scale*singular_01);
	// 			K.set(2*ind_i+1, 2*ind_i+1, scale*singular_11);
			

	// 			continue;
	// 		 }
		 	

	// 	 	int ind_j = j/2;
			
	// 		Vec2 x(points[i], points[i+1]);
	// 		Vec2 y(points[j], points[j+1]);
	// 		Vec2 r = x-y;  

	// 		double r0 = r.a[0];
	// 		double r1 = r.a[1];
			
			
	// 		K.set(2*ind_i  , 2*ind_j  , weights[ind_j]*scale*(log(1.0/r.norm()) + (r0*r0/r.dot(r)) ));
	// 		K.set(2*ind_i+1, 2*ind_j  , weights[ind_j]*scale*(                    (r0*r1/r.dot(r)) ));
	// 		K.set(2*ind_i  , 2*ind_j+1, weights[ind_j]*scale*(                    (r1*r0/r.dot(r)) ));
	// 		K.set(2*ind_i+1, 2*ind_j+1, weights[ind_j]*scale*(log(1.0/r.norm()) + (r1*r1/r.dot(r)) ));
		


	// 	}
	// }

//DOUBLE LAYER POTENTIAL BELOW
	double scale = 1.0 / (M_PI);

	for(unsigned int i=0; i<points.size(); i+=2){
	//	omp_set_num_threads(4);
	//	#pragma omp parallel for 	
	
		for(unsigned int j=0; j<points.size(); j+=2){
			
			int ind_i = i/2;
			
			if(i==j){

				double t0 = -normals[i+1];
				double t1 =  normals[i];
				K.set(2*ind_i  , 2*ind_i  , -0.5 - 0.5*curvatures[ind_i]*weights[ind_i]*scale*t0*t0);
				K.set(2*ind_i+1, 2*ind_i  ,      - 0.5*curvatures[ind_i]*weights[ind_i]*scale*t0*t1);
				K.set(2*ind_i  , 2*ind_i+1,      - 0.5*curvatures[ind_i]*weights[ind_i]*scale*t0*t1);
				K.set(2*ind_i+1, 2*ind_i+1, -0.5 - 0.5*curvatures[ind_i]*weights[ind_i]*scale*t1*t1);
			

				continue;
			 }

		 	int ind_j = j/2;
			
			Vec2 x(points[i], points[i+1]);
			Vec2 y(points[j], points[j+1]);
			Vec2 r = x-y;  

			double r0 = r.a[0];
			double r1 = r.a[1];
			Vec2 n(normals[j], normals[j+1]);
			double potential = weights[ind_j]*scale*(r.dot(n))/(pow(r.dot(r),2));
			
			
			K.set(2*ind_i  , 2*ind_j  , potential*r0*r0);
			K.set(2*ind_i+1, 2*ind_j  , potential*r0*r1);
			K.set(2*ind_i  , 2*ind_j+1, potential*r0*r1);
			K.set(2*ind_i+1, 2*ind_j+1, potential*r1*r1);
		


		}
	}
}


void Initialization::Stokes_InitializeDomainKernel(ie_Mat& K, std::vector<double>& points,
										std::vector<double>& normals, 
										std::vector<double>& weights, double min, double max, int test_size,
										int (*out_of_domain)(Vec2& a)){
	//columns for phi (aka dofs), rows for spatial domain
	

	// int dofs = points.size()/2;
	double scale = 1.0 / (M_PI);
	omp_set_num_threads(4);

	// min-=d/2;
	// max+=d/2;
	//#pragma omp parallel for 

	for(int i=0; i<test_size*test_size; i++){
		double x0 = i/test_size;
		x0 = min + (x0*(max-min))/test_size;
		double x1 = i%test_size;
		x1 = min + (x1*(max-min))/test_size;
		Vec2 x(x0, x1);
		for(unsigned int j=0; j<points.size(); j+=2){
			int ind_j = j/2;
			
			Vec2 y(points[j], points[j+1]);

			if(out_of_domain(x)){
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


void Initialization::Stokes_InitializeBoundary(ie_Mat& f, std::vector<double>& normals){
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