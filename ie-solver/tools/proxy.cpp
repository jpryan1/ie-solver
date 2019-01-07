#include "ie_solver_tools.h"
#include "ie-solver/kernel.h"
#include <cmath>

namespace ie_solver{

void IeSolverTools::make_id_mat(const Kernel& K, ie_Mat& mat, QuadTree& tree,
	QuadTreeNode* node){

	double cntr_x = node->corners[0] + node->side_length/2.0;
	double cntr_y = node->corners[1] + node->side_length/2.0;
	if(is_stokes){
		// // TODO fix this later, match laplace please, make cleaner please
		// mat = ie_Mat(200, node->interaction_lists.active_box.size());
		// make_stokes_proxy_mat(mat, cntr_x, cntr_y, node->side_length*2, 
		// 	node->interaction_lists.active_box );
	}
	else{
		ie_Mat box_near = K(node->interaction_lists.near, 
			node->interaction_lists.active_box);
		ie_Mat proxy = ie_Mat(100, node->interaction_lists.active_box.size());
		make_proxy_mat(proxy, cntr_x, cntr_y, node->side_length*2, tree,
			node->interaction_lists.active_box);
		mat = ie_Mat(node->interaction_lists.near.size() + 100, 
			node->interaction_lists.active_box.size());
		mat.set_submatrix(0, node->interaction_lists.near.size(), 
						  0, node->interaction_lists.active_box.size(), box_near);
		mat.set_submatrix(node->interaction_lists.near.size(), 
						100 + node->interaction_lists.near.size(),
						  0, node->interaction_lists.active_box.size(), proxy);

	}
}


void IeSolverTools::make_proxy_mat(ie_Mat& pxy, double cntr_x, double cntr_y, 
	double r, QuadTree& tree, const std::vector<unsigned int>& box_indices){
	
	//each row is a pxy point, cols are box dofs

	double scale = 1.0 / (2*M_PI);

	for(int i=0; i<100; i++){
		
		double ang = 2*M_PI*i*0.01;
		Vec2 p(cntr_x+r*cos(ang), cntr_y+r*sin(ang));
		for(unsigned int j=0; j < box_indices.size(); j++){
			unsigned int box_index = box_indices[j];
			Vec2 q(tree.boundary->points[2*box_index], tree.boundary->points[2*box_index+1]);

			Vec2 r = p - q;
			Vec2 n(tree.boundary->normals[2*box_index], tree.boundary->normals[2*box_index+1]);

			double potential = -tree.boundary->weights[box_index]*scale*(r.dot(n))/(r.dot(r));
			pxy.set(i, j, potential);
		}
	}
}


void IeSolverTools::make_stokes_proxy_mat(ie_Mat& pxy, double cntr_x, 
	double cntr_y, double r, const std::vector<unsigned int>& box_indices){

	// double scale = 1.0 / (M_PI);
	// for(int i = 0; i < 200; i += 2){
	// 	double ang = 2 * M_PI * i * 0.01;
	// 	Vec2 p(cntr_x + r*cos(ang), cntr_y + r*sin(ang));
	// 	for(unsigned int j = 0; j < box_indices.size(); j++){
			
	// 		unsigned int box_index = box_indices[j];
	// 		Vec2 q(points[2*box_index], points[2*box_index+1]);

	// 		Vec2 r = p-q;
	// 		Vec2 n(normals[2*box_index], normals[2*box_index+1]);
	// 		double r0 = r.a[0];
	// 		double r1 = r.a[1];
			
	// 		double potential = weights[box_index] * scale * 
	// 			(r.dot(n)) / (pow(r.dot(r),2));
		
	// 		if( box_index % 2 == 0 ){
	// 			pxy.set(i  , j, potential*r0*r0);
	// 			pxy.set(i+1, j, potential*r0*r1);
	// 		}else{
	// 			pxy.set(i  , j, potential*r0*r1);
	// 			pxy.set(i+1, j, potential*r1*r1);
	// 		}
	// 	}
	// }
}

} // namespace ie_solver
