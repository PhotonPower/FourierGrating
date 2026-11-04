// NVIDIA's Windows CUDA Toolkit import libraries (cudart.lib, cusolver.lib, ...) are built with
// MSVC and embed objects that reference the MSVC CRT's /GS stack-cookie and SEH helpers. MinGW's
// runtime does not provide these symbols, so linking against the import libs with GCC/MinGW fails
// with "undefined reference to `__security_check_cookie'" (and `__GSHandlerCheck`/`__security_cookie`).
// They guard code paths MinGW never executes, so stubbing them out as no-ops is safe here.
#include <cstdint>

extern "C" {
std::uintptr_t __security_cookie = 0;
void __security_check_cookie(std::uintptr_t) {}
int __GSHandlerCheck() { return 0; }
}
