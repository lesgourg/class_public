#include "exceptions.h"
#include <string>
#include <sstream>
#include <stdexcept>

void ThrowInvalidArgumentIf(bool condition, std::string string_for_printf, ...) {
  if (condition) {
    va_list args;
    va_start(args, string_for_printf);
    char error_message[100];
    sprintf(error_message, string_for_printf.c_str(), args);
    throw std::invalid_argument(error_message);
  }
}

void ThrowInvalidArgument(std::string string_for_printf, ...) {
  va_list args;
  va_start(args, string_for_printf);
  char error_message[100];
  sprintf(error_message, string_for_printf.c_str(), args);
  throw std::invalid_argument(error_message);
}
