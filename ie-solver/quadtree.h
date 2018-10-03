#ifndef _QUADTREE_H_
#define _QUADTREE_H_

#include "common.h"

#define MAX_LEAF_DOFS 128

// TODO all headers should have variable names explicit. 

namespace ie_solver{

struct QuadTreeNode{
	int id;
	static int id_count;

	QuadTreeNode *tl, *tr, *bl, *br;
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
	~QuadTreeNode(){
		for(QuadTreeNode* child: children){
			if(child){
				delete child;
			}
		}
	}
};
//TODO replace this with a typedef?
struct QuadTreeLevel{
	std::vector<QuadTreeNode*> nodes;
};

class QuadTree {

public:
	~QuadTree(){
		if(root){
			delete root;
		}
		for(QuadTreeLevel* level : levels){
			if(level){
				delete level;
			}
		}
	}
	QuadTreeNode* root;
	std::vector<QuadTreeLevel*> levels;
	std::vector<double> pts;

	double min, max;

	void initialize_tree(const std::vector<double>& points, bool is_stokes);
	void recursive_add(QuadTreeNode* node, double x, double y,
		unsigned int mat_ind);
	void node_subdivide(QuadTreeNode* node);
	void add_index(std::vector<unsigned int>& r, unsigned int ind);

	int which_field(double x, double y, QuadTreeNode* node);

	// void print();
	// void rec_print(QuadTreeNode*);

private:
	bool is_stokes_;
};

} // namespace ie_solver

#endif