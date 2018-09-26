#include "quadtree.h"

namespace ie_solver{

int QuadTreeNode::id_count = 0;

void QuadTree::initialize_tree(std::vector<double> points, bool is_stokes) {
		
	is_stokes_ = is_stokes;

	min = points[0];
	max = points[1];
	for( double point : points ){
		if(point < min) min = point;
		if(point > max) max = point;
	}

	//this is a tad silly
	min -= (0.1 + rand()*(1.0/RAND_MAX)*1e-3);
	max += (0.1 + rand()*(1.0/RAND_MAX)*1e-3);

	root = new QuadTreeNode();
	root->level = 0;
	//Root is the box [0,1]x[0,1]
	//BL
	root->corners[0] = min;
	root->corners[1] = min;
	//TL
	root->corners[2] = min;
	root->corners[3] = max;
	//TR
	root->corners[4] = max;
	root->corners[5] = max;
	//BR
	root->corners[6] = max;
	root->corners[7] = min;

	root->side = max-min;

	QuadTreeLevel* level_one = new QuadTreeLevel();
	level_one->nodes.push_back(root);
	levels.push_back(level_one);



	pts.resize(points.size());
	memcpy(&(pts[0]), &(points[0]), sizeof(double)*points.size());
	//all near and far ranges should be taken care of by this guy
	for(unsigned int i=0; i<points.size(); i+=2){
		recursive_add(this->root, points[i], points[i+1], i/2);
	}

	//This is a naive method and should probably be made smarter

	//After the quadtree is made, we go through the points and check if they should be in the near field of EVERY node


	// for(int i=0; i<points.size(); i+=2){
	// 	double x = points[i];
	// 	double y = points[i+1];
	// 	for(int j=0; j<levels.size(); j++){
	// 		QuadTreeLevel* current_level = levels[j];
	// 		for(int k=0; k< current_level->nodes.size(); k++){

	// 			QuadTreeNode* current_node = current_level->nodes[k];
				
	// 			int field = which_field(x,y,current_node);
	// 			if( field==1){
	// 				add_index(current_node->box.near_range, i/2);
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
	// 		if(is_stokes_){
	// 			assert(current_node->box.box_range.size()+
	// 			current_node->box.near_range.size()+
	// 			current_node->box.far_range.size() == points.size());
	// 		}else{
	// 			assert(current_node->box.box_range.size()+
	// 			current_node->box.near_range.size()+
	// 			current_node->box.far_range.size() == points.size()/2);
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

				double dist =  sqrt(pow(node_a->corners[0]-node_b->corners[0],2) 
								   +pow(node_a->corners[1]-node_b->corners[1],2));	
				//just need to check if the distance of the BL corners is <=2sqrt(2)
				if(  dist < node_a->side*sqrt(2)+1e-5){
					node_a->neighbors.push_back(node_b);
					node_b->neighbors.push_back(node_a);
				}
				else if(  dist < node_a->side*2*sqrt(2)+1e-5){
					node_a->far_neighbors.push_back(node_b);
					node_b->far_neighbors.push_back(node_a);
				}
			}		
		}
	}
	for(unsigned int j=0; j<levels.size(); j++){
		QuadTreeLevel* current_level = levels[j];
		for( QuadTreeNode* node_a : current_level->nodes ){
			for(QuadTreeNode* neighbor: node_a->neighbors){
				neighbor->box.near_range.insert(neighbor->box.near_range.end(),
					node_a->box.box_range.begin(), node_a->box.box_range.end());
			}		
		}
	}
}


void QuadTree::add_index(std::vector<unsigned int>& r, unsigned int ind){

	if(is_stokes_){
		r.push_back(2*ind);
		r.push_back(2*ind+1);
	}else{
		r.push_back(ind);
	}
}


int QuadTree::which_field(double x, double y, QuadTreeNode* node){
	double bl_x = node->corners[0];
	double bl_y = node->corners[1];
	double side = node->corners[3]-node->corners[1];
	if(x>= bl_x && x < bl_x+side && y>= bl_y && y< bl_y+side){
		return 0;//inside box
	} 
	if(x>= bl_x-side && x < bl_x+2*side && y>= bl_y-side && y< bl_y+2*side){
		return 1;//near field
	}
	 return 2;//far field
}


void QuadTree::recursive_add(QuadTreeNode* node, double x, double y, unsigned int mat_ind){
	
	add_index(node->box.box_range, mat_ind);
	
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
		if(node->is_leaf && node->box.box_range.size() > MAX_LEAF_DOFS){
			node_subdivide(node);
		}
	}
}

// 1) gives node its four children
// 2) puts these children in their proper level
// 3) gives these children their corners
void QuadTree::node_subdivide(QuadTreeNode* node){

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

	bl->side = node->side/2.0;
	br->side = node->side/2.0;
	tl->side = node->side/2.0;
	tr->side = node->side/2.0;



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

	for(unsigned int dof = 0; dof < node->box.box_range.size(); dof++){

		unsigned int ind = node->box.box_range[dof];
		//theres a non-trivial difference between the stokes
		// and laplace case here, so we split this up
		if(is_stokes_){
			if(ind%2==1) continue; //because in stokes case, two "dof"s correspond to one point
			double x = pts[ind];
			double y = pts[ind+1];
		
			if(x < midx && y< midy){
				bl->box.box_range.push_back(ind);
				bl->box.box_range.push_back(ind+1);
			}
			else if(x < midx && y>= midy){
				tl->box.box_range.push_back(ind);
				tl->box.box_range.push_back(ind+1);
			}else if(x >= midx && y< midy){
				br->box.box_range.push_back(ind);
				br->box.box_range.push_back(ind+1);
			}else{
				tr->box.box_range.push_back(ind);
				tr->box.box_range.push_back(ind+1);
			}
		}
		else{
			double x = pts[2*ind];
			double y = pts[2*ind+1];
	

			if(x < midx && y< midy){
				bl->box.box_range.push_back(ind);
			}
			else if(x < midx && y>= midy){
				tl->box.box_range.push_back(ind);
			}else if(x >= midx && y< midy){
				br->box.box_range.push_back(ind);
			}else{
				tr->box.box_range.push_back(ind);
			}
		}
	}

	if(bl->box.box_range.size() > MAX_LEAF_DOFS){
		node_subdivide(bl);
	}
	if(tl->box.box_range.size() > MAX_LEAF_DOFS){
		node_subdivide(tl);
	}
	if(tr->box.box_range.size() > MAX_LEAF_DOFS){
		node_subdivide(tr);
	}
	if(br->box.box_range.size() > MAX_LEAF_DOFS){
		node_subdivide(br);
	}
}

//TODO logging here
void QuadTree::print(){
	// std::cout<<"printing points first"<<std::endl;
	// for(int i=0; i<pts.size(); i+=2){

	// 	std::cout<<(i/2)<<": "<<pts[i]<<" "<<pts[i+1]<<std::endl;
	// }
	// std::cout<<"Done printing points"<<std::endl;
	// rec_print(root);


	// printf("Now printing levels\n");
	// for(int i=0; i<levels.size(); i++){
	// 	QuadTreeLevel* current_level = levels[i];
	// 	printf("\nLevel %d\n", i);
	// 	for(int j=0; j<current_level->nodes.size(); j++){
	// 		QuadTreeNode* current_node = current_level->nodes[j];
	// 		for(int k=0; k<current_node->box.box_range.size(); k++){
	// 			std::cout<<current_node->box.box_range[k]<<" ";
	// 		}
	// 	}
	// }
}


//TODO logging here
void QuadTree::rec_print(QuadTreeNode* n){

	// //printf("Current node has box with %d dofs\n", n->box.box_range.size());
	// printf("Corners are ");
	// for(int i=0; i<8; i+=2){
	// 	std::cout<<n->corners[i]<<" "<<n->corners[i+1]<<"   ";
	// }std::cout<<std::endl;


	// std::cout<<"Box"<<std::endl;
	// for(unsigned int i=0; i<n->box.box_range.size(); i++){
	// 	std::cout<<n->box.box_range[i]<<" ";
	// }
	// std::cout<<"\n"<<"Near"<<std::endl;
	// for(unsigned int i=0; i<n->box.near_range.size(); i++){
	// 	std::cout<<n->box.near_range[i]<<" ";
	// }std::cout<<std::endl;
	// if(n->bl){
	// 	printf("Node has children!\n");
	// 	printf("bl...\n");
	// 	rec_print(n->bl);
	// 	printf("tl...\n");
	// 	rec_print(n->tl);
	// 	printf("tr...\n");
	// 	rec_print(n->tr);
	// 	printf("br...\n");
	// 	rec_print(n->br);
	// }printf("Going up!\n");

}

} // namespace ie_solver