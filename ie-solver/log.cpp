#include <iostream>
#include "log.h"

namespace ie_solver{

LOG::LOG_LEVEL LOG::log_level_ = LOG::LOG_LEVEL::INFO_;

void LOG::INFO(const std::string& message){
	if(log_level_ <= INFO_){
		std::cout << "INFO: " << message << std::endl;
	}
}


void LOG::WARNING(const std::string& message){
	if(log_level_ <= WARNING_){
		std::cout << "WARNING: " << message << std::endl;
	}
}


void LOG::ERROR(const std::string& message){
	if(log_level_ <= ERROR_){
		std::cout << "ERROR: " << message << std::endl;
	}
}

} // namespace ie_solver