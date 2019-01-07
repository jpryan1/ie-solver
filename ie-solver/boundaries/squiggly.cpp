#include "squiggly.h"


void squiggly(int N, std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights){



	int side_scale = N;
	int side_points = 6*side_scale;

	std::vector<double> integrals(side_scale);


	// bottom
	for(int i=0; i<side_points; i++){
		//t ranges from 0 to 3pi
		double t = (3.0*M_PI*i)/side_points;
		
		points.push_back(t);
		points.push_back(-sin(t));


		double x0 = 1;
		double x1 = -cos(t);
		double nrm = sqrt(pow(x0,2)+pow(x1,2));
		normals.push_back(x1/nrm);
		normals.push_back(-x0/nrm);


		curvatures.push_back(sin(t)/pow(1+pow(cos(t),2) , 1.5));
		
	}


	// right
	for(int i=0; i<side_points; i++){
		//t ranges from 0 to 3pi
		double t = (3.0*M_PI*i)/side_points;
		
		points.push_back(3*M_PI + sin(t));
		points.push_back(t);

		double x0 = cos(t);
		double x1 = 1;
		double nrm = sqrt(pow(x0,2)+pow(x1,2));
		normals.push_back(x1/nrm);
		normals.push_back(-x0/nrm);

		curvatures.push_back(sin(t)/pow(1+pow(cos(t),2) , 1.5));
		
	}

	// top
	for(int i=0; i<side_points; i++){
		//t ranges from 0 to 3pi
		double t = (3.0*M_PI*i)/side_points;
		
		points.push_back(3*M_PI - t);
		points.push_back(3*M_PI + sin(3*M_PI-t));
	


		double x0 = -1;
		double x1 = -cos(3*M_PI-t);
		double nrm = sqrt(pow(x0,2)+pow(x1,2));
		normals.push_back(x1/nrm);
		normals.push_back(-x0/nrm);



		curvatures.push_back(sin(3*M_PI-t)/pow(1+pow(cos(3*M_PI-t),2) , 1.5));
		
	
	}

	// left
	for(int i=0; i<side_points; i++){
		//t ranges from 0 to 3pi
		double t = (3.0*M_PI*i)/side_points;
		
		points.push_back(-sin(3*M_PI-t));
		points.push_back(3*M_PI - t);
		


		double x0 = cos(3*M_PI-t);
		double x1 = -1;
		double nrm = sqrt(pow(x0,2)+pow(x1,2));
		normals.push_back(x1/nrm);
		normals.push_back(-x0/nrm);

		curvatures.push_back(sin(3*M_PI-t)/pow(1+pow(cos(3*M_PI-t),2) , 1.5));
	
	}






	// now we need to figure out some arclengths
	double k = 1.0/sqrt(2.0);

	std::vector<double> E(side_scale+1);
	E[0] = 0; //this isn't used, just to make the indexing a bit easier
	for(int i=1; i<side_scale+1; i++){
		double ang = 0.5*M_PI*((i+0.0) / side_scale);
		E[i] = sqrt(2)*boost::math::ellint_2(k, ang);

	}

	for(int k=0; k<12; k++){
		weights.push_back(E[1]);

		for(int i=1; i<side_scale; i++){
			weights.push_back((E[i+1] - E[i-1])/2.0);		
		}


		weights.push_back(E[side_scale] - E[side_scale-1]);
		for(int i=side_scale-1; i>=1; i--){
			weights.push_back((E[i+1] - E[i-1])/2.0);
		}	
	}









	// //for exterior problem
	// for(int i=0; i<normals.size(); i++){
	// 	normals[i] = -1*normals[i];
	// }
}





int out_of_squiggly(Vec2& a){
	double x = a.a[0];
	double y = a.a[1];
	double eps = 1;
	if( y - eps > -sin(x) 		   && 
		y + eps <  sin(x) + 3*M_PI && 
		x - eps > -sin(y) 		   && 
		x + eps <  sin(y) + 3*M_PI    ) return 0;

	return 1;



	// if(y+eps>-sin(x) && y -eps< sin(x)+3*M_PI && x+eps>-sin(y) && x-eps<sin(y)+3*M_PI) return 1;

	// return 0;
}
