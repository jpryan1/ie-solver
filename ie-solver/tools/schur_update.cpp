#include "ie_solver_tools.h"
#include <cassert>

namespace ie_solver{

void IeSolverTools::get_all_schur_updates(ie_Mat& updates, 
	const std::vector<unsigned int>& BN, const QuadTreeNode* node){
	assert(node != nullptr && "get_all_schur_updates fails on null node.");
	assert(BN.size() > 0 && "get_all_schur_updates needs positive num of DOFs");

	if(!node->is_leaf) get_descendents_updates(updates, BN, node);
	
	if(strong_admissibility){
		for(QuadTreeNode* neighbor : node->neighbors){
			if(neighbor->schur_updated) get_update(updates, BN, neighbor);
			if(!neighbor->is_leaf) get_descendents_updates(updates, BN, neighbor);	
		}
	}
}


void IeSolverTools::get_descendents_updates(ie_Mat& updates, 
	const std::vector<unsigned int>& BN, const QuadTreeNode* node){
	assert(node != nullptr && "get_descendents_updates fails on null node.");
	assert(!node->is_leaf && 
		"get_descendents_updates must be called on non-leaf.");

	//by assumption, node is not a leaf	
	for(QuadTreeNode* child : node->children){
		if(child->schur_updated) get_update(updates, BN, child);
		if(!child->is_leaf) get_descendents_updates(updates, BN, child);
	}
}


void IeSolverTools::get_update(ie_Mat& update, 
	const std::vector<unsigned int>& BN, const QuadTreeNode* node){
	// node needs to check all its dofs against BN, enter interactions into 
	// corresponding locations
	// node only updated its own BN dofs, and the redundant ones are no longer 
	// relevant, so we only care about child's SN dofs
	std::vector<unsigned int> sn = node->interaction_lists.skelnear;
	// First create a list of Dofs that are also in node's skelnear, 
	// and with each one give the index in skelnear and the index in BN
	std::vector<unsigned int> BN_indices;
	std::vector<unsigned int> sn_indices;
	for(unsigned int i=0; i<sn.size(); i++){
		for(unsigned int j=0; j<BN.size(); j++){
			if(BN[j] == sn[i]){
				sn_indices.push_back(i);
				BN_indices.push_back(j);
			}
		}
	}	
	for(unsigned int i=0; i<BN_indices.size(); i++){
		for(unsigned int j=0; j<BN_indices.size(); j++){
			update.addset(BN_indices[i], BN_indices[j],
				node->schur_update.get(sn_indices[i], sn_indices[j]));
		}
	}
}

} // namespace ie_solver