#ifndef _CIRCLE_H_
#define _CIRCLE_H_


#include "boundary.h"

namespace ie_solver{

class Circle : public Boundary {
	public:
		void initialize(int N);
		bool is_in_domain(Vec2& a);

};

} // namespace ie_solver

#endif