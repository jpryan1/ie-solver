#ifndef _INTERACTION_LISTS_H_
#define _INTERACTION_LISTS_H_

namespace ie_solver{

struct InteractionLists {
// Note, maybe redundant and skel could be Range<int>s, that would make things faster maybe
	std::vector<unsigned int> box_range;
	std::vector<unsigned int> redundant_range;
	std::vector<unsigned int> skel_range;
	std::vector<unsigned int> near_range;
	std::vector<unsigned int> skelnear_range;
	// TODO this variable name is unhelpful, change it
	std::vector<unsigned int> p;

};

} // namespace

#endif