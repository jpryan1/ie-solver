#include "ie_solver_tools.h"
#include <cassert>

namespace ie_solver{

//TODO URGENTLY skelnear needs to be replaced as a variable name

//TODO here and elsewhere check on passing by reference

IeSolverTools::IeSolverTools(double id_tol, bool strong_admissibility_, 
	bool is_stokes_) {
	assert(id_tol > 0 && "id_tol must be greater than one to init tools.");
	this->id_tol = id_tol;
	strong_admissibility = strong_admissibility_;
	is_stokes = is_stokes_;
}


void IeSolverTools::populate_active_box(QuadTreeNode* node){

	// this function removes from the box any DoFs which have already been made 
	// redundant. It involves a bunch of annoying C++ functions and probably 
	// would look nicer in matlab.

	// populate active_box
	if(!node->is_leaf){
		for(QuadTreeNode* child: node->children){
			if(child->schur_updated){
				for(unsigned int i : child->interaction_lists.skel){
					node->interaction_lists.active_box.push_back(i);
				}
			}else{
				for(unsigned int i : child->interaction_lists.active_box){
					node->interaction_lists.active_box.push_back(i);
				}
			}
		}
	}else{
		node->interaction_lists.active_box = node->interaction_lists.original_box;
	}
}


void IeSolverTools::set_rs_ranges(InteractionLists& interaction_lists, 
	const std::vector<unsigned int>& prm, unsigned int sk, unsigned int rd) {

	assert(prm.size() == sk+rd);

	for(unsigned int i = 0; i<sk; i++){
    	interaction_lists.skel.push_back(interaction_lists.active_box[prm[i]]);
    	interaction_lists.permutation.push_back(prm[i]);	
    }
    for(unsigned int i = sk; i<sk+rd; i++){
        interaction_lists.redundant.push_back(interaction_lists.active_box[prm[i]]);
        interaction_lists.permutation.push_back(prm[i]);
    }
}

void IeSolverTools::set_skelnear_range(InteractionLists& interaction_lists){
	for(unsigned int i = 0; i<interaction_lists.skel.size(); i++){
        interaction_lists.skelnear.push_back(interaction_lists.skel[i]);
    }
    if(strong_admissibility){
	    for(unsigned int i=0; i<interaction_lists.near.size(); i++){
	        interaction_lists.skelnear.push_back(interaction_lists.near[i]);
	    }
	}
}

} // namespace ie_solver
