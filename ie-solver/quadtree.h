#ifndef _QUADTREE_H_
#define _QUADTREE_H_

#include "common.h"

#define MAX_LEAF_DOFS 128

// TODO all headers should have variable names explicit. 

namespace ie_solver{

struct QuadTreeNode{

	int id;
	static int id_count;

	QuadTreeNode* tl;
	QuadTreeNode* tr;
	QuadTreeNode* bl;
	QuadTreeNode* br;
	QuadTreeNode* children[4];

	std::vector<QuadTreeNode*> neighbors, far_neighbors;

	ie_Mat T, L, U, D_r, schur_update;
	Box box;
	
	double side;
	unsigned int level;
	bool is_leaf;
	bool schur_updated = false;
	
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
};

struct QuadTreeLevel{
	std::vector<QuadTreeNode*> nodes;
	
};

class QuadTree {

public:
	QuadTreeNode* root;
	std::vector<QuadTreeLevel*> levels;
	std::vector<double> pts;

	double min, max;

	int which_field(double, double, QuadTreeNode* );

	void initialize_tree(std::vector<double>, bool);
	void recursive_add(QuadTreeNode*, double, double, unsigned int);
	void node_subdivide(QuadTreeNode*);
	void add_index(std::vector<unsigned int>& r, unsigned int ind);


	void print();
	void rec_print(QuadTreeNode*);

private:
	bool is_stokes_;
};

} // namespace ie_solver

#endif