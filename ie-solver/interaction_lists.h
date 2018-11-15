#ifndef _INTERACTION_LISTS_H_
#define _INTERACTION_LISTS_H_

namespace ie_solver{

struct InteractionLists {

	std::vector<unsigned int> box;
	std::vector<unsigned int> redundant;
	std::vector<unsigned int> skel;
	std::vector<unsigned int> near;
	std::vector<unsigned int> skelnear;

	std::vector<unsigned int> permutation;

};

} // namespace

#endif