#ifndef _ROUNDED_SQUARE_H_
#define _ROUNDED_SQUARE_H_


#include "boundary.h"

namespace ie_solver{

class RoundedSquare : public Boundary {
	public:

		void initialize(int N);

		bool is_in_domain(Vec2& a);

};

} // namespace ie_solver

#endif