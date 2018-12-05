#include <fstream>

#include "common.h"
#include "initialization.h"
#include "ie_solver_tools.h"
#include "quadtree.h"
#include "circle.h"
#include "rounded_square.h"

#define DEFAULT_NUM_DISCRETIZATION_POINTS 1000
#define DEFAULT_ID_TOL 1e-6
#define TEST_SIZE 100
#define TIMING_ITERATIONS 5
namespace ie_solver{

LOG::LOG_LEVEL LOG::log_level_ = LOG::LOG_LEVEL::INFO_;


struct ie_solver_config{
	enum Pde {
		LAPLACE,
		STOKES
	};
	enum Shape {
		CIRCLE,
		ROUNDED_SQUARE
	};

	int N = DEFAULT_NUM_DISCRETIZATION_POINTS;
	double id_tol = DEFAULT_ID_TOL;
	Pde pde = LAPLACE;
	Shape shape = CIRCLE;
	bool scaling = false;
	std::ofstream n_scaling_output, error_scaling_output;

};


void get_circle_stokes_solution(double min, double max, ie_Mat& domain, 
	int (*out_of_shape)(Vec2& a)){
	//min and max describe the box inside which our sample is taken. 

	// double total_err = 0;
	// double max_err   = 0;
	// double true_norm = 0;
	
	for(int i=0; i<TEST_SIZE*TEST_SIZE; i++){

		//find out the point corresponding to this index
		double x0 = i/TEST_SIZE;
		x0 = min + (x0*(max-min))/TEST_SIZE;
		double y0 = i%TEST_SIZE;
		y0 = min + (y0*(max-min))/TEST_SIZE;
		Vec2 v(x0,y0);

		if(out_of_shape(v)){
			continue;
		}
		//now we get a radius associated with this index
		x0-=0.5;
		y0-=0.5;
		
		// double r = sqrt(x0*x0+y0*y0);

		domain.set(2*i  , 0, -4*y0/(1));
		domain.set(2*i+1, 0,  4*x0/(1));	
	}	
}

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


void boundary_integral_solve(const ie_solver_config& config, 
	std::vector<double>* skel_times = nullptr){

	void (*make_shape) (int, std::vector<double>&, std::vector<double>&, 
		std::vector<double>&, std::vector<double>&);
	int (*out_of_shape)(Vec2& a);

	bool is_stokes = (config.pde == ie_solver_config::STOKES);
	bool is_time_trial = (skel_times != nullptr);
	int N = config.N;
	double id_tol = config.id_tol;
	switch(config.shape){
		case ie_solver_config::CIRCLE:
			make_shape = circle;
			out_of_shape = out_of_circle;
			break;
		case ie_solver_config::ROUNDED_SQUARE:
			make_shape = rounded_square;
			out_of_shape = out_of_rounded_square;
			break;
	}

	// // TODO insert comment here explaining why this is necessary
	// if(!timing && N > 10000){
	// 	printf("Turn down N or disable accuracy checking please\n");
	// 	return;
	// }

	std::vector<double> points, normals, curvatures, weights;
	make_shape(N, points, normals, curvatures, weights);
	if(!is_time_trial){
		write_boundary_to_file(points);
	}
	int dofs = points.size() / 2;

	QuadTree quadtree;
	quadtree.initialize_tree(points, is_stokes); 
	if(!is_time_trial){
		quadtree.write_quadtree_to_file();
	}
	
	// Consider making init instead of constructor for readability
	IeSolverTools ie_solver_tools(id_tol, points, normals, weights, is_stokes);

	ie_Mat K;
	K.is_stokes = is_stokes;
	K.is_dynamic = true;
	K.load(&points, &normals, &curvatures, &weights);

	if(is_time_trial){
		double elapsed = 0;
		for(int i=0; i < TIMING_ITERATIONS; i++){
			double start = omp_get_wtime();
			ie_solver_tools.Skeletonize(K, quadtree);
			double end = omp_get_wtime();
			elapsed += (end-start);

			quadtree.reset();
		}
		skel_times->push_back(elapsed/TIMING_ITERATIONS);
		return;
	}
	else{
		ie_solver_tools.Skeletonize(K, quadtree);
	}
	// Now we calculate the solution to the PDE inside the domain by setting
	// up the relevant linear system. 

	int dim = is_stokes ? 2 : 1;
	
	std::vector<double> domain_points;
	get_domain_points(domain_points, quadtree.min, quadtree.max);
	
	ie_Mat K_domain = ie_Mat(dim*TEST_SIZE*TEST_SIZE, dim*dofs);
	Initialization init;
	init.InitializeDomainKernel(K_domain, points, normals, weights, 
 		domain_points, TEST_SIZE, out_of_shape, is_stokes);
	
	
	ie_Mat f(dim*dofs, 1);
	// TODO get rid of these damn if(is_stokes) statements, push them to the 
	// functions
 	if(is_stokes){
 		init.Stokes_InitializeBoundary(f, normals); 
 		// notice here we are passing the normals since the flow will just be 
 		// unit tangent to the boundary. 
 	}else{
 		init.InitializeBoundary(f, points);
 	}
	
 	ie_Mat phi(dim*dofs, 1);
	ie_solver_tools.Solve(K, quadtree, phi, f);
	
	// This will be done as a sparse mat vec in the future, for now we do 
	// dense matvec

	ie_Mat domain(dim*TEST_SIZE*TEST_SIZE, 1);
	Matmul::ie_gemv(NORMAL, 1., K_domain, phi, 0., domain);
	write_solution_to_file(domain, domain_points, is_stokes);
	
}


int parse_input_into_config(int argc, char** argv, ie_solver_config& config){
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
		else if(!strcmp(argv[i], "-scaling")){
			config.scaling = true;
		}
		else if(!strcmp(argv[i], "-e")){
			if(i < argc - 1){
				config.id_tol = std::stof(argv[i+1]);
			}
			i++;
		}else if(!strcmp(argv[i], "-boundary")){
			if(i < argc - 1){
				if(!strcmp(argv[i+1], "CIRCLE")){
					config.shape = ie_solver_config::CIRCLE;
				}
				else if(!strcmp(argv[i+1], "ROUNDED_SQUARE")){
					config.shape = ie_solver_config::ROUNDED_SQUARE;
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
				" {number of nodes} -e {ID error tolerance}\nOmitting an arg"
				" triggers a default value.");
			return 0;
		}
	}
	return 1;
	// TODO after parsing, print out input configuration
}

} // namespace ie_solver


int main(int argc, char** argv){
	
	// TODO allow for command line args for setting parameters

	ie_solver::ie_solver_config config;
	if(!ie_solver::parse_input_into_config(argc, argv, config)){
		return 1;
	}

	if(!config.scaling){
		ie_solver::boundary_integral_solve(config);
	}

	if(config.scaling){
		std::vector<double> n_times, eps_times;
		int scale_n[] = {5000, 6000, 7000, 8000, 9000};
		double scale_eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
		for(int n : scale_n){
			config.N = n;
			config.id_tol = DEFAULT_ID_TOL;
			ie_solver::boundary_integral_solve(config, &n_times);
		}
		for(double eps : scale_eps){
			config.N = DEFAULT_NUM_DISCRETIZATION_POINTS*10;
			config.id_tol = eps;
			ie_solver::boundary_integral_solve(config, &eps_times);
		}
		ie_solver::write_times_to_files(scale_n, n_times, scale_eps, eps_times);
	}

	return 0;
}

