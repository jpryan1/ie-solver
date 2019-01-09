#include "rounded_square.h"
#include "ie-solver/log.h"

#include <cmath>

namespace ie_solver{

void RoundedSquare::initialize(int N, int bc_enum){
	if(bc_enum != BoundaryCondition::SINGLE_ELECTRON){
		LOG::ERROR("Circle boundary can only do single electron bc currently.");
	}
	boundary_condition = ie_Mat(N, 1);

	// This square will have side length 0.5 and will have BL corner 0.25, 0.25

	// Sides will go from 0.3 to 0.7, rounded corners will have the rest

	// each side will have 128 discretization points, each corner will have 16 points

	// radius of rounded corner is 0.05
	// TODO obviously N should be the total num of points, not that/40
	int SCALE = 20*(N/20); // Now SCALE is a multiple of 20 (?)

	int NUM_SIDE_POINTS = SCALE/5;
	int NUM_CORN_POINTS = SCALE/20;


	double side_AL = 0.8/(NUM_SIDE_POINTS+1);
	double corner_AL = 0.05*M_PI/(2*NUM_CORN_POINTS-2);
	
	int bc_idx = 0;

	// bottom side

	for(int i=1; i<NUM_SIDE_POINTS+1; i++){
		double x = 0.1+(i/(NUM_SIDE_POINTS+1.0))*0.8;
		double y = 0.05;
		points.push_back(x);
		points.push_back(y);
		normals.push_back(0);
		normals.push_back(-1);
		curvatures.push_back(0);
		weights.push_back(side_AL);

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

	//bottom right corner
	for(int i=0; i<NUM_CORN_POINTS; i++){
		
		double t = i*M_PI / (2*NUM_CORN_POINTS-2.0); //This ranges from 0 to pi/2, inclusive
		double ang = (3.0*M_PI/2.0)+t;
		double x = 0.9+0.05*cos(ang);
		double y = 0.1+0.05*sin(ang);
		points.push_back(x);
		points.push_back(y);

		normals.push_back(cos(ang));
		normals.push_back(sin(ang));

		curvatures.push_back(1/(0.05));
	
		if(i==0 || i==(NUM_CORN_POINTS)-1){
			weights.push_back(0.5*(side_AL+corner_AL));
		}else{
			weights.push_back(corner_AL);
		}

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

	// right side
	for(int i=1; i<NUM_SIDE_POINTS+1; i++){
		double x = 0.95;
		double y = 0.1+(i/(NUM_SIDE_POINTS+1.0))*0.8;
		points.push_back(x);
		points.push_back(y);

		normals.push_back(1);
		normals.push_back(0);

		curvatures.push_back(0);

		weights.push_back(side_AL);

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

	//top right corner
	for(int i=0; i<NUM_CORN_POINTS; i++){
		double t = i*M_PI / (2*NUM_CORN_POINTS-2.0); //This ranges from 0 to pi/2, inclusive
		
		double ang = t;
		double x = 0.9+0.05*cos(ang);
		double y = 0.9+0.05*sin(ang);
		points.push_back(x);
		points.push_back(y);

		normals.push_back(cos(ang));
		normals.push_back(sin(ang));

		curvatures.push_back(1/(0.05));

		if(i==0 || i==(NUM_CORN_POINTS)-1){

			weights.push_back(0.5*(side_AL+corner_AL));
		}else{
			weights.push_back(corner_AL);
		}
		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

	// top side
	for(int i=NUM_SIDE_POINTS; i>0; i--){
		double x = 0.1+(i/(NUM_SIDE_POINTS+1.0))*0.8;
		double y = 0.95;
		points.push_back(x);
		points.push_back(y);

		normals.push_back(0);
		normals.push_back(1);

		curvatures.push_back(0);

		weights.push_back(side_AL);

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

		//top left corner
	for(int i=0; i<NUM_CORN_POINTS; i++){
		double t = i*M_PI / (2*NUM_CORN_POINTS-2.0); //This ranges from 0 to pi/2, inclusive
		
		double ang = (M_PI/2.0)+t;
		double x = 0.1+0.05*cos(ang);
		double y = 0.9+0.05*sin(ang);
		points.push_back(x);
		points.push_back(y);

		normals.push_back(cos(ang));
		normals.push_back(sin(ang));

		curvatures.push_back(1/(0.05));

		if(i==0 || i==(NUM_CORN_POINTS)-1){
			weights.push_back(0.5*(side_AL+corner_AL));
		}else{
			weights.push_back(corner_AL);
		}

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}


	// left side
	for(int i=NUM_SIDE_POINTS; i>0; i--){
		double x = 0.05;
		double y = 0.1+(i/(NUM_SIDE_POINTS+1.0))*0.8;
		points.push_back(x);
		points.push_back(y);
		
		normals.push_back(-1);
		normals.push_back(0);

		curvatures.push_back(0);

		weights.push_back(side_AL);

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

	//bottom left corner
	for(int i=0; i<NUM_CORN_POINTS; i++){

		double t = i*M_PI / (2*NUM_CORN_POINTS-2.0); //This ranges from 0 to pi/2, inclusive
		double ang = M_PI+t;
		double x = 0.1+0.05*cos(ang);
		double y = 0.1+0.05*sin(ang);
		points.push_back(x);
		points.push_back(y);

		normals.push_back(cos(ang));
		normals.push_back(sin(ang));

		curvatures.push_back(1/(0.05));

		if(i==0 || i==(NUM_CORN_POINTS)-1){
			weights.push_back(0.5*(side_AL+corner_AL));
		}else{
			weights.push_back(corner_AL);
		}

		double potential = log(sqrt( pow(x+2,2)+pow(y+2,2)))/(2*M_PI);
		boundary_condition.set(bc_idx++, 0, potential);
	}

}


bool RoundedSquare::is_in_domain(Vec2& a){
	double* v = a.a;

	double eps = 1e-2;

	if( fabs(v[0] - 0.5) > 0.45 - eps 
		|| fabs(v[1] - 0.5) > 0.45 - eps){
		return false;
	}

	if(v[0]<=0.95 + eps && v[0] + eps>=0.05 && v[1] + eps>=0.1 && v[1]<=0.9 + eps) return true;
	if(v[0]<=0.9  + eps && v[0] + eps>=0.1 && v[1] + eps>=0.05 && v[1]<=0.95 + eps) return true;

	double min = 1;
	Vec2 bl(0.1, 0.1);
	Vec2 br(0.9, 0.1);
	Vec2 tl(0.1, 0.9);
	Vec2 tr(0.9, 0.9);

	min = fmin(min, (bl-a).norm());
	min = fmin(min, (tr-a).norm());
	min = fmin(min, (tl-a).norm());
	min = fmin(min, (br-a).norm());
	if(min + eps>0.05) return false;
	return true;
}

} // namespace
