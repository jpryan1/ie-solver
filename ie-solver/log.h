// Copyright 2019 John Paul Ryan
#ifndef _LOG_H_
#define _LOG_H_

#include <string>

namespace ie_solver {

struct LOG {
  enum LOG_LEVEL {
    INFO_,
    WARNING_,
    ERROR_
  };

  static LOG_LEVEL log_level_;

  static void INFO(const std::string& message);
  static void WARNING(const std::string& message);
  static void ERROR(const std::string& message);
};

}  // namespace ie_solver

#endif
