#ifndef _MOVING_CIRCLE_H_
#define _MOVING_CIRCLE_H_


#include "boundary.h"

namespace ie_solver{

class MovingCircle : public Boundary {
	public:
		// TOFO Do I actually need to define these things again, or is it 
		// implicit? 
		
		void initialize(int N);

		bool is_in_domain(Vec2& a);

		double circle_x, circle_y;

};

} // namespace ie_solver

#endif