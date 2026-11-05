# Common warning set for all FourierGrating targets.
add_library(fg_warnings INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(fg_warnings INTERFACE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion>)
  if(FG_WARNINGS_AS_ERRORS)
    target_compile_options(fg_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Werror>)
  endif()
endif()
if(MSVC)
  # MSVC port (Plan section "M12: Build"): standard conformance, UTF-8 sources, M_PI etc.
  target_compile_options(fg_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/W4 /permissive- /utf-8 /bigobj
                                                   /external:anglebrackets /external:W0>)
  target_compile_definitions(fg_warnings INTERFACE _USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
  if(FG_WARNINGS_AS_ERRORS)
    target_compile_options(fg_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/WX>)
  endif()
endif()
