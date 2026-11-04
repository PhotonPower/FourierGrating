# Common warning set for all FourierGrating targets.
add_library(fg_warnings INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(fg_warnings INTERFACE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion>)
  if(FG_WARNINGS_AS_ERRORS)
    target_compile_options(fg_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Werror>)
  endif()
endif()
