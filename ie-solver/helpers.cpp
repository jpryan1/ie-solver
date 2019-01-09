#include "helpers.h"
#include <cassert>
#include "log.h"
#include <cmath>
#include <iostream>
#include <string.h>

namespace ie_solver{

// void get_circle_stokes_solution(double min, double max, ie_Mat& domain, 
// 	bool (*is_in_domain)(Vec2& a)){
// 	//min and max describe the box inside which our sample is taken. 

// 	// double total_err = 0;
// 	// double max_err   = 0;
// 	// double true_norm = 0;
	
// 	for(int i=0; i<TEST_SIZE*TEST_SIZE; i++){

// 		//find out the point corresponding to this index
// 		double x0 = i / TEST_SIZE;
// 		x0 = min + (x0 * (max - min)) / TEST_SIZE;
// 		double y0 = i%TEST_SIZE;
// 		y0 = min + (y0 * (max - min)) / TEST_SIZE;
// 		Vec2 v(x0,y0);

// 		if(!is_in_domain(v)){
// 			continue;
// 		}
// 		//now we get a radius associated with this index
// 		x0-=0.5;
// 		y0-=0.5;
		
// 		// double r = sqrt(x0*x0+y0*y0);

// 		domain.set(2*i  , 0, -4 * y0);
// 		domain.set(2*i+1, 0,  4 * x0);
// 	}	
// }


// TODO put write functions in another file or something
void write_boundary_to_file(std::vector<double>& points){
	assert(points.size() % 2 == 0);
	std::ofstream output;
	output.open("output/data/ie_solver_boundary.txt");
	if(output.is_open()){
		for(unsigned int i=0; i<points.size(); i += 2){
			output << points[i] << "," << points[i+1]<<std::endl;
		}
		output.close();
	}else{
		LOG::ERROR("Failed to open boundary output file!");
	}
}


void write_potential_to_file(){

}


void write_times_to_files(int* scale_n, std::vector<double>& n_times, 
	double* scale_eps, std::vector<double>& eps_times){

	std::ofstream n_output, e_output;
	n_output.open("output/data/ie_solver_n_scaling.txt");
	e_output.open("output/data/ie_solver_e_scaling.txt");


	if(n_output.is_open()){
		for(unsigned int i = 0; i < n_times.size(); i++){
			n_output << scale_n[i]<<","<<n_times[i]<<std::endl;
		}
		n_output.close();
	}else{
		printf("Failed to open n output file!\n");
	}

	if(e_output.is_open()){
		for(unsigned int i = 0; i < eps_times.size(); i++){
			e_output << scale_eps[i]<<","<<eps_times[i]<<std::endl;
		}
		e_output.close();
	}else{
		printf("Failed to open e output file!\n");
	}

}


void write_solution_to_file(ie_Mat& domain, std::vector<double>& domain_points,
	bool is_stokes){
	
	assert(domain.height() > 0 && domain.width() == 1);
	std::ofstream output;
	output.open("output/data/ie_solver_solution.txt");

	int dim = (is_stokes ? 2 : 1);
	// if it is 2D solution, each row is 2 numbers, ie a vector. else a scalar.

	int points_index = 0;

	if(output.is_open()){
		for(unsigned int i=0; i<domain.height(); i += dim){
			output << domain_points[points_index] << "," <<
				domain_points[points_index+1] << ",";
			points_index += 2;
			output<<domain.get(i, 0);
			if(is_stokes) output << "," << domain.get(i+1,0);
			output << std::endl;
		}
		output.close();
	}else{
		printf("Failed to open solution output file!\n");
	}
}


void get_domain_points(std::vector<double>& points, double min, double max){
	for(int i=0; i<TEST_SIZE; i++){
		double x = min + ((i+0.0)/TEST_SIZE)*(max-min);
		for(int j = 0; j < TEST_SIZE; j++){
			double y = min + ((j+0.0)/TEST_SIZE)*(max-min);
			points.push_back(x);
			points.push_back(y);
		}
	}
}


void check_laplace_solution(ie_Mat& domain, double id_tol, 
	std::vector<double>& domain_points, Boundary* boundary){
	double max = 0;
	for(unsigned int i=0; i<domain_points.size(); i+=2){
		double x0 = domain_points[i];
		double x1 = domain_points[i+1];
		Vec2 x(x0,x1);
		if (!boundary->is_in_domain(x)){
			continue;
		}
		double potential = log(sqrt( pow(x0+2,2)+pow(x1+2,2)))/(2*M_PI);
		if(std::isnan(domain.get(i/2, 0))){
			continue;
		}
		max = std::max(max, std::abs(std::abs(potential - domain.get(i/2, 0))));
	}
	if(max <= id_tol*100){
		std::cout << "Solution looks good."<< std::endl;
	}else{
		std::cout << "Error max err too big: " << max << std::endl;
	}
}


int parse_input_into_config(int argc, char** argv, ie_solver_config& config){
	// The default boundary is Circle, set it here.
	config.boundary.reset(new Circle());

	for(int i=1; i<argc; i++){
		if(!strcmp(argv[i], "-pde")){
			if(i < argc - 1){
				if(!strcmp(argv[i+1], "STOKES")){
					config.pde = ie_solver_config::STOKES;
				}
				else if(!strcmp(argv[i+1], "LAPLACE")){
					config.pde = ie_solver_config::LAPLACE;
				}
				else{
					LOG::ERROR("Unrecognized pde: " + std::string(argv[i+1])
						+ "\n Acceptable pdes: STOKES, LAPLACE");
					return 0;
				}
			}
			i++;
		}else if(!strcmp(argv[i], "-N")){
			if(i < argc - 1){
				config.N = std::stoi(argv[i+1]);
			}
			i++;
		}
		else if(!strcmp(argv[i], "-strong")){
			config.admissibility = ie_solver_config::STRONG;
		}
		else if(!strcmp(argv[i], "-scaling")){
			config.scaling = true;
		}
		else if(!strcmp(argv[i], "-h")){
			LOG::log_level_ = LOG::LOG_LEVEL::INFO_;
			LOG::INFO("\n\tusage: ./ie-solver -pde {LAPLACE|STOKES}"
				" -boundary {CIRCLE|ROUNDED_SQUARE} -N"
				" {number of nodes} -e {ID error tolerance} {-scaling}"
				" \nOmitting an arg triggers a default value.");
			return 0;
		}
		else if(!strcmp(argv[i], "-e")){
			if(i < argc - 1){
				config.id_tol = std::stof(argv[i+1]);
			}
			i++;
		}else if(!strcmp(argv[i], "-boundary")){
			if(i < argc - 1){
				if(!strcmp(argv[i+1], "CIRCLE")){
					config.boundary.reset(new Circle());
				}
				else if(!strcmp(argv[i+1], "ROUNDED_SQUARE")){
					config.boundary.reset(new RoundedSquare());
				}
				else if(!strcmp(argv[i+1], "ROUNDED_SQUARE_WITH_BUMP")){
					config.boundary.reset(new RoundedSquareWithBump());
				}
				else{
					LOG::ERROR("Unrecognized boundary: "+ std::string(argv[i+1])
						+ "\n Acceptable boundaries: CIRCLE, ROUNDED_SQUARE");
					return 0;
				}
			}
			i++;
		}else{
			LOG::ERROR("Unrecognized argument: " + std::string(argv[i]) + 
				"usage: ./ie-solver -pde {LAPLACE|STOKES}"
				" -boundary {CIRCLE|ROUNDED_SQUARE} -N"
				" {number of nodes} -e {ID error tolerance} {-scaling}"
				" \nOmitting an arg triggers a default value.");
			return 0;
		}
	}
	return 1;
	// TODO after parsing, print out input configuration
}

} // namespace ie_solver