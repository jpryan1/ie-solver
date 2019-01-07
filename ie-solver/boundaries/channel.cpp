#include "channel.h"


void channel(std::vector<double>& points, std::vector<double>& normals, std::vector<double>& curvatures,
	std::vector<double>& weights){

	int flat_points = 4;
	int sine_points = 16*flat_points;

	//top left flat

	for(int i=0; i<flat_points; i++){
		double t = (i+0.0)/flat_points;
		points.push_back( -0.25 + 0.25*t);
		points.push_back( 4.5);

		normals.push_back(0);
		normals.push_back(1);

		curvatures.push_back(0);

	}


	//top sine


	for(int i=0; i<sine_points; i++){
		//t ranges from 0 to 3pi
		double t = (2*M_PI*i)/sine_points;
		
		points.push_back(t);
		points.push_back(cos(t) + 3.5);


		double x0 = 1;
		double x1 = -sin(t);
		double nrm = sqrt(pow(x0,2)+pow(x1,2));
		normals.push_back(-x1/nrm);
		normals.push_back(x0/nrm);


		curvatures.push_back(-cos(t)/pow(1+pow(sin(t),2) , 1.5));
		
	}

	//top right flat
	// NOTE THE <=flat_points, this means we are now adding the top right corner point!
	for(int i=0; i<=flat_points; i++){

		double t = (i+0.0)/flat_points;

		points.push_back( 2*M_PI + 0.25*t);
		points.push_back( 4.5);

		normals.push_back(0);
		normals.push_back(1);

		curvatures.push_back(0);

	}

	// bottom left flat
	for(int i=0; i<flat_points; i++){
		double t = (i+0.0)/flat_points;
		points.push_back( -0.25 + 0.25*t);
		points.push_back(0);

		normals.push_back(0);
		normals.push_back(-1);

		curvatures.push_back(0);

	}


	//bottom sine


	for(int i=0; i<sine_points; i++){
		//t ranges from 0 to 3pi
		double t = (2*M_PI*i)/sine_points;
		
		points.push_back(t);
		points.push_back(-cos(t)+1);

		double x0 = 1;
		double x1 = sin(t);
		double nrm = sqrt(pow(x0,2)+pow(x1,2));
		normals.push_back(x1/nrm);
		normals.push_back(-x0/nrm);



		curvatures.push_back(sin(t)/pow(1+pow(cos(t),2) , 1.5));
		
	}

	//bottom right flat
	// NOTE THE <=flat_points, this means we are now adding the top right corner point!
	for(int i=0; i<=flat_points; i++){

		double t = (i+0.0)/flat_points;

		points.push_back( 2*M_PI + 0.25*t);
		points.push_back( 0);

		normals.push_back(0);
		normals.push_back(-1);

		curvatures.push_back(0);

	}


	//some init work for sine arclength calculation later on
	double k = 1.0/sqrt(2.0);
	int sine_slice = sine_points/4;
	std::vector<double> E(sine_slice);
	E[0] = 0; //this isn't used, just to make the indexing a bit easier
	for(int i=1; i<sine_slice+1; i++){
		double ang = 0.5*M_PI*((i+0.0) / sine_slice);
		E[i] = sqrt(2)*boost::math::ellint_2(k, ang);
	}





	//top left flat
	weights.push_back(0.125/flat_points);
	for(int i=1; i<flat_points; i++) weights.push_back(0.25/flat_points);


	//top sine

	weights.push_back( (E[sine_slice] - E[sine_slice-1] + (0.25/flat_points)) /2.0);
	for(int i=sine_slice-1; i>0; i--){
		weights.push_back((E[i+1] - E[i-1])/2.0);
	}
	weights.push_back(E[1]);
	for(int i=1; i<sine_slice; i++){
		weights.push_back((E[i+1] - E[i-1])/2.0);
	}	


	weights.push_back(E[sine_slice] - E[sine_slice-1] );
	
	// Now we have all the weights we need, for simplicity, we now enter the rest of the 
	// top boundary by referring to the first half
	int half_top = weights.size();
	for(int i=half_top-2; i>=0; i--) weights.push_back(weights[i]);
	// and now the bottom has the same weights as the top
	int top_size = weights.size();
	for(int i=0; i<top_size; i++) weights.push_back(weights[i]);



	assert(weights.size() == points.size()/2);



	// for(int i=2; i<points.size()-2; i+=2){
	// 	int w_ind = i/2;
	// 	double dist1 = sqrt(pow(points[i] - points[i-2],2) + pow(points[i+1]-points[i-1],2));
	// 	double dist2 = sqrt(pow(points[i] - points[i+2],2) + pow(points[i+1]-points[i+3],2));
	// 	double avg = (dist1+dist2)/2.0;
	// 	printf("%f %f\n", weights[w_ind] , avg);
	// }

}





int out_of_channel(Vec2& a){
	double x = a.a[0];
	double y = a.a[1];
	double eps = 1e-1;
	
	if(x < -0.25 || x > M_PI*2 + 0.25 || y < 0 || y > 4.5) return 1;

	if( x > 0 && x < 2*M_PI && (y < -cos(x) || y> cos(x)+3.5)) return 1;
	return 0;

}
