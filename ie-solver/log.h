#ifndef _LOG_H_
#define _LOG_H_

#include <stdio.h>

namespace ie_solver{

struct LOG{
	enum LOG_LEVEL {
		INFO_, 
		WARNING_,
		ERROR_
	};
	static LOG_LEVEL log_level_;

	static void INFO(const std::string& message){
		if(log_level_ <= INFO_){
			std::cout << "INFO: " << message << std::endl;
		}
	}
	static void WARNING(const std::string& message){
		if(log_level_ <= WARNING_){
			std::cout << "WARNING: " << message << std::endl;
		}
	}
	static void ERROR(const std::string& message){
		if(log_level_ <= ERROR_){
			std::cout << "ERROR: " << message << std::endl;
		}
	}
};

} // namespace ie_solver

#endif 