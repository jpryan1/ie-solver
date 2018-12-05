#ifndef _QUADTREE_H_
#define _QUADTREE_H_

#include "common.h"

#define MAX_LEAF_DOFS 128

// TODO all headers should have variable names explicit. 

namespace ie_solver{

struct QuadTreeNode{
	
	static unsigned int id_count;

	unsigned int id, level;	
	bool is_leaf, schur_updated = false;
	double side_length;
	
	QuadTreeNode *tl, *tr, *bl, *br;
	QuadTreeNode* children[4];
	std::vector<QuadTreeNode*> neighbors;

	InteractionLists interaction_lists;
	ie_Mat T, L, U, D_r, schur_update;
	
	//format is {BL, TL, TR, BR}
	double corners[8];

	QuadTreeNode(){
		id = id_count++;
		is_leaf = true;
		tl=NULL;
		tr=NULL;
		bl=NULL;
		br=NULL;
		for(int i=0; i<4; i++) children[i] = NULL;
	}
	~QuadTreeNode(){
		for(QuadTreeNode* child: children){
			if(child){
				delete child;
			}
		}
	}
};

struct QuadTreeLevel{
	std::vector<QuadTreeNode*> nodes;
};

class QuadTree {

public:

	bool is_stokes;
	double min, max;
	std::vector<double> pts;
	QuadTreeNode* root;
	std::vector<QuadTreeLevel*> levels;

	~QuadTree(){
		if(root){
			delete root;
		}
		for(QuadTreeLevel* level : levels){
			if(level){
				delete level;
			}
		}
		levels.clear();

	}
	
	void reset();

	void initialize_tree(const std::vector<double>& points, bool is_stokes);
	void recursive_add(QuadTreeNode* node, double x, double y,
		unsigned int mat_ind);
	void node_subdivide(QuadTreeNode* node);
	void add_index(std::vector<unsigned int>& r, unsigned int ind);

	int which_field(double x, double y, QuadTreeNode* node);

	void write_quadtree_to_file();

	// void print();
	// void rec_print(QuadTreeNode*);

};

} // namespace ie_solver

#endif