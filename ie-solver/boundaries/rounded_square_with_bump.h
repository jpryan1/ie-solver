#ifndef _ROUNDED_SQUARE_WITH_BUMP_H_
#define _ROUNDED_SQUARE_WITH_BUMP_H_


#include "boundary.h"

namespace ie_solver{

class RoundedSquareWithBump : public Boundary {
	public:
		
		enum BoundaryCondition{
			SINGLE_ELECTRON
		};
		
		void initialize(int N, int bc_enum);
		bool is_in_domain(Vec2& a);
};

} // namespace ie_solver

#endif