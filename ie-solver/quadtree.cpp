#include "quadtree.h"
#include <fstream>
#include "log.h"
#include <cassert>
#include <cmath>

namespace ie_solver{

unsigned int QuadTreeNode::id_count = 0;

void QuadTree::initialize_tree(Boundary* boundary_, bool is_stokes_) {


	assert(boundary_->points.size() > 0 && "Number of boundary->points to init tree cannot be 0.");
	// TODO later we can assert that the boundary->points are a multiple of dimension
	this->boundary = boundary_;
	this->is_stokes = is_stokes_;

	min = boundary->points[0];
	max = boundary->points[0];

	for( double point : boundary->points ){ 
		if(point < min) min = point;
		if(point > max) max = point;
	}
	//this is a tad silly
	// double tree_min = min - (0.1 + rand()*(1.0/RAND_MAX)*1e-3);
	// double tree_max = max + (0.1 + rand()*(1.0/RAND_MAX)*1e-3);

	// JUST FOR DEBUGGING
	double tree_min = 0;//min - (0.1 + rand()*(1.0/RAND_MAX)*1e-3);
	double tree_max = 1;//max + (0.1 + rand()*(1.0/RAND_MAX)*1e-3);

	root = new QuadTreeNode();
	root->level = 0;
	root->parent = nullptr;
	//Root is the box [0,1]x[0,1]
	//BL
	root->corners[0] = tree_min;
	root->corners[1] = tree_min;
	//TL
	root->corners[2] = tree_min;
	root->corners[3] = tree_max;
	//TR
	root->corners[4] = tree_max;
	root->corners[5] = tree_max;
	//BR
	root->corners[6] = tree_max;
	root->corners[7] = tree_min;

	root->side_length = tree_max-tree_min;

	QuadTreeLevel* level_one = new QuadTreeLevel();

	level_one->nodes.push_back(root);
	levels.push_back(level_one);



	// all near and far ranges should be taken care of by this guy
	for(unsigned int i=0; i<boundary->points.size(); i+=2){
		recursive_add(this->root, boundary->points[i], boundary->points[i+1], i/2);
	}

	// This is a naive method and should probably be made smarter

	// After the quadtree is made, we go through the boundary->points and check if they 
	// should be in the near field of EVERY node


	// for(int i=0; i<boundary->points.size(); i+=2){
	// 	double x = boundary->points[i];
	// 	double y = boundary->points[i+1];
	// 	for(int j=0; j<levels.size(); j++){
	// 		QuadTreeLevel* current_level = levels[j];
	// 		for(int k=0; k< current_level->nodes.size(); k++){

	// 			QuadTreeNode* current_node = current_level->nodes[k];
				
	// 			int field = which_field(x,y,current_node);
	// 			if( field==1){
	// 				add_index(current_node->interaction_lists.near, i/2);
	// 			}
	// 			else if(field==2){
	// 				add_index(current_node->box.far_range, i/2);
	// 			}

	// 		}
	// 	}
	// }
	// //just a quick check on the validity of the tree
	// for(int j=0; j<levels.size(); j++){
	// 	QuadTreeLevel* current_level = levels[j];
	// 	for(int k=0; k< current_level->nodes.size(); k++){

	// 		QuadTreeNode* current_node = current_level->nodes[k];
	// 		if(is_stokes){
	// 			assert(current_node->interaction_lists.original_box.size()+
	// 			current_node->interaction_lists.near.size()+
	// 			current_node->box.far_range.size() == boundary->points.size());
	// 		}else{
	// 			assert(current_node->interaction_lists.original_box.size()+
	// 			current_node->interaction_lists.near.size()+
	// 			current_node->box.far_range.size() == boundary->points.size()/2);
	// 		}
	// 	}
	// }


	// make neighbor lists in a stupid way
	for(unsigned int j=0; j<levels.size(); j++){
		QuadTreeLevel* current_level = levels[j];
		for(unsigned int k=0; k< current_level->nodes.size(); k++){

			QuadTreeNode* node_a = current_level->nodes[k];
			for(unsigned int l=k+1; l<current_level->nodes.size(); l++){
				QuadTreeNode* node_b = current_level->nodes[l];

				double dist = sqrt( pow(node_a->corners[0]-node_b->corners[0],2) 
					+ pow(node_a->corners[1]-node_b->corners[1],2));	
				// just need to check if the distance of the BL corners 
				// is <=s*sqrt(2)
				if(  dist < node_a->side_length*sqrt(2)+1e-5){
					node_a->neighbors.push_back(node_b);
					node_b->neighbors.push_back(node_a);
				}
				// else if(  dist < node_a->side_length*2*sqrt(2)+1e-5){
				// 	node_a->far_neighbors.push_back(node_b);
				// 	node_b->far_neighbors.push_back(node_a);
				// }
			}		
		}
	}
	for(unsigned int j=0; j<levels.size(); j++){
		QuadTreeLevel* current_level = levels[j];
		for( QuadTreeNode* node_a : current_level->nodes ){
			for(QuadTreeNode* neighbor: node_a->neighbors){
				neighbor->interaction_lists.near.insert(
					neighbor->interaction_lists.near.end(), 
					node_a->interaction_lists.original_box.begin(), 
					node_a->interaction_lists.original_box.end());
			}		
		}
	}
}


void QuadTree::recursive_add(QuadTreeNode* node, double x, double y,
	unsigned int mat_ind){
	
	assert(node != nullptr && "recursive_add fails on null node.");

	add_index(node->interaction_lists.original_box, mat_ind);
	
	//figure out which child
	double midx = ((node->corners[6]-node->corners[0])/2.0 ) + node->corners[0];
	double midy = ((node->corners[3]-node->corners[1])/2.0 ) + node->corners[1];
	
	QuadTreeNode* child;
	if(x < midx && y< midy){
		child = node->bl;
	}
	else if(x < midx && y>= midy){
		child = node->tl;
	}
	else if(x >= midx && y< midy){
		child = node->br;
	}
	else{
		child = node->tr;
	}
	//does that child exist?
	if(child){
		recursive_add(child, x, y, mat_ind);
	}
	else{

		//do we need one? 
		//If this node is exploding and needs children
		if(node->is_leaf && node->interaction_lists.original_box.size() > MAX_LEAF_DOFS){
			node_subdivide(node);
		}
	}
}

// 1) gives node its four children
// 2) puts these children in their proper level
// 3) gives these children their corners
void QuadTree::node_subdivide(QuadTreeNode* node){

	assert(node != nullptr && "node_subdivide fails on null node.");
	node->is_leaf = false;
			
	QuadTreeNode *bl, *br, *tl, *tr;
	double midx = ((node->corners[6]-node->corners[0])/2.0 ) + node->corners[0];
	double midy = ((node->corners[3]-node->corners[1])/2.0 ) + node->corners[1];



	bl = new QuadTreeNode();
	bl->corners[0] = node->corners[0];
	bl->corners[1] = node->corners[1];
	bl->corners[2] = node->corners[2];
	bl->corners[3] = midy;
	bl->corners[4] = midx;
	bl->corners[5] = midy;
	bl->corners[6] = midx;
	bl->corners[7] = node->corners[7];



	tl = new QuadTreeNode();
	tl->corners[0] = node->corners[0];
	tl->corners[1] = midy;
	tl->corners[2] = node->corners[2];
	tl->corners[3] = node->corners[3];
	tl->corners[4] = midx;
	tl->corners[5] = node->corners[5];
	tl->corners[6] = midx;
	tl->corners[7] = midy;

	tr = new QuadTreeNode();
	tr->corners[0] = midx;
	tr->corners[1] = midy;
	tr->corners[2] = midx;
	tr->corners[3] = node->corners[3];
	tr->corners[4] = node->corners[4];
	tr->corners[5] = node->corners[5];
	tr->corners[6] = node->corners[6];
	tr->corners[7] = midy;


	br = new QuadTreeNode();
	br->corners[0] = midx;
	br->corners[1] = node->corners[1];
	br->corners[2] = midx;
	br->corners[3] = midy;
	br->corners[4] = node->corners[4];
	br->corners[5] = midy;
	br->corners[6] = node->corners[6];
	br->corners[7] = node->corners[7];


	bl->level = node->level+1;
	br->level = node->level+1;
	tl->level = node->level+1;
	tr->level = node->level+1;

	bl->side_length = node->side_length/2.0;
	br->side_length = node->side_length/2.0;
	tl->side_length = node->side_length/2.0;
	tr->side_length = node->side_length/2.0;

	bl->parent = node;
	br->parent = node;
	tl->parent = node;
	tr->parent = node;


	 if(levels.size()<node->level+2){
		QuadTreeLevel* new_level = new QuadTreeLevel();
		levels.push_back(new_level);
		new_level->nodes.push_back(bl);
		new_level->nodes.push_back(br);
		new_level->nodes.push_back(tr);
		new_level->nodes.push_back(tl);
	}else{
		
		levels[node->level+1]->nodes.push_back(bl);
		levels[node->level+1]->nodes.push_back(br);
		levels[node->level+1]->nodes.push_back(tr);
		levels[node->level+1]->nodes.push_back(tl);
	}
	node->bl = bl;
	node->tl = tl;
	node->tr = tr;
	node->br = br;
	node->children[0] = bl;
	node->children[1] = tl;
	node->children[2] = tr;
	node->children[3] = br;

	for(unsigned int dof = 0; dof < node->interaction_lists.original_box.size(); dof++){

		unsigned int ind = node->interaction_lists.original_box[dof];
		//theres a non-trivial difference between the stokes
		// and laplace case here, so we split this up
		if(is_stokes){
			if(ind%2==1) continue;
			 //because in stokes case, two "dof"s correspond to one point
			double x = boundary->points[ind];
			double y = boundary->points[ind+1];
		
			if(x < midx && y< midy){
				bl->interaction_lists.original_box.push_back(ind);
				bl->interaction_lists.original_box.push_back(ind+1);
			}
			else if(x < midx && y>= midy){
				tl->interaction_lists.original_box.push_back(ind);
				tl->interaction_lists.original_box.push_back(ind+1);
			}else if(x >= midx && y< midy){
				br->interaction_lists.original_box.push_back(ind);
				br->interaction_lists.original_box.push_back(ind+1);
			}else{
				tr->interaction_lists.original_box.push_back(ind);
				tr->interaction_lists.original_box.push_back(ind+1);
			}
		}
		else{
			double x = boundary->points[2*ind];
			double y = boundary->points[2*ind+1];
	

			if(x < midx && y< midy){
				bl->interaction_lists.original_box.push_back(ind);
			}
			else if(x < midx && y>= midy){
				tl->interaction_lists.original_box.push_back(ind);
			}else if(x >= midx && y< midy){
				br->interaction_lists.original_box.push_back(ind);
			}else{
				tr->interaction_lists.original_box.push_back(ind);
			}
		}
	}

	if(bl->interaction_lists.original_box.size() > MAX_LEAF_DOFS){
		node_subdivide(bl);
	}
	if(tl->interaction_lists.original_box.size() > MAX_LEAF_DOFS){
		node_subdivide(tl);
	}
	if(tr->interaction_lists.original_box.size() > MAX_LEAF_DOFS){
		node_subdivide(tr);
	}
	if(br->interaction_lists.original_box.size() > MAX_LEAF_DOFS){
		node_subdivide(br);
	}
}


void QuadTree::add_index(std::vector<unsigned int>& r, unsigned int ind){

	if(is_stokes){
		r.push_back(2*ind);
		r.push_back(2*ind+1);
	}else{
		r.push_back(ind);
	}
}


int QuadTree::which_field(double x, double y, QuadTreeNode* node){
	assert( node != nullptr && "which_field fails on null node.");

	double bl_x = node->corners[0];
	double bl_y = node->corners[1];
	double side_length = node->corners[3]-node->corners[1];
	if(x>= bl_x && x < bl_x+side_length && y>= bl_y && y< bl_y+side_length){
		return 0;//inside box
	} 
	if(x>= bl_x-side_length && x < bl_x+2*side_length && y>= bl_y-side_length 
		&& y< bl_y+2*side_length){
		return 1;//near field
	}
	 return 2;//far field
}


void QuadTree::write_quadtree_to_file(){

	std::ofstream output;
	output.open("output/data/ie_solver_tree.txt");
	if(output.is_open()){
		for(QuadTreeLevel* level : levels){
			for(QuadTreeNode* node : level->nodes){
				if(node->is_leaf){
					output << node->corners[0] << "," << node->corners[1] 
						<< "," << node->side_length << std::endl;
				}
			}
		}
		output.close();
	}else{
		LOG::ERROR("Failed to open output file!");
	}
}

void QuadTree::perturb(){
	QuadTreeLevel* last_level = levels[levels.size()-1];
	QuadTreeNode* perturbed;
	// Find a leaf with more than 10 DoFs to perturb.
	for(QuadTreeNode* node: last_level->nodes){
		if(node->interaction_lists.original_box.size()>10){
			perturbed = node;
			break;
		}
	}
	double min_x = perturbed->corners[0];
	double min_y = perturbed->corners[1];
	// Pick 10 new random boundary->points for the first 10 in the leaf.
	for(int i=0; i<10; i++){
		int pt_index = perturbed->interaction_lists.original_box[i];
		double randx = (rand()+0.0) / RAND_MAX;
		double randy = (rand()+0.0) / RAND_MAX;
		randx = randx * perturbed->side_length + min_x;
		randy = randy * perturbed->side_length + min_y;
		boundary->points[2*pt_index] = randx;
		boundary->points[2*pt_index+1] = randy;
	}

	// Go up tree, marking parents and their neighbors.
	QuadTreeNode* current = perturbed;
	while(current != nullptr){
		current->schur_updated = false;
		current->interaction_lists.active_box.clear();
		current->interaction_lists.skel.clear();
		current->interaction_lists.skelnear.clear();
		current->interaction_lists.redundant.clear();

		for(QuadTreeNode* neighbor : current->neighbors){
			neighbor->schur_updated = false;
			neighbor->interaction_lists.active_box.clear();
			neighbor->interaction_lists.skel.clear();
			neighbor->interaction_lists.skelnear.clear();
			neighbor->interaction_lists.redundant.clear();
		}
		current = current->parent;
	}
}

void QuadTree::reset(){
	if(root){
		delete root;
	}
	for(QuadTreeLevel* level : levels){
		if(level){
			delete level;
		}
	}
	levels.clear();
	QuadTreeNode::id_count = 0;
	initialize_tree(boundary, is_stokes);	
}


} // namespace ie_solver