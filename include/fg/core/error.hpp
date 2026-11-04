// FourierGrating - error handling.
#pragma once

#include <stdexcept>
#include <string>

namespace fg {

/// Base exception for all FourierGrating errors (always carries context).
class Error : public std::runtime_error {
 public:
  explicit Error(const std::string& what) : std::runtime_error(what) {}
};

/// Thrown by functions that are specified in the plan but not implemented yet.
class NotImplemented : public Error {
 public:
  explicit NotImplemented(const std::string& what) : Error("not implemented: " + what) {}
};

}  // namespace fg

#define FG_REQUIRE(cond, msg)                                                      \
  do {                                                                             \
    if (!(cond))                                                                   \
      throw ::fg::Error(std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                        ": " + (msg));                                             \
  } while (0)
