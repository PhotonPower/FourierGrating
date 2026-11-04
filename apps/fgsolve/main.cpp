// fgsolve - command line front end of FourierGrating.
//   fgsolve --version
//   fgsolve --info                 (backends, build options)
//   fgsolve <config.json>          (diffraction efficiencies and optional near fields)
#include <cstring>
#include <iostream>

#include "fg/core/error.hpp"
#include "fg/io/config.hpp"
#include "fg/linalg/backend.hpp"
#include "fg/version.hpp"

int main(int argc, char** argv) {
  if (argc < 2 || std::strcmp(argv[1], "--help") == 0) {
    std::cout << "usage: fgsolve [--version | --info | <config.json>]\n";
    return argc < 2 ? 1 : 0;
  }
  try {
    if (std::strcmp(argv[1], "--version") == 0) {
      std::cout << "FourierGrating " << FG_VERSION << "\n";
      return 0;
    }
    if (std::strcmp(argv[1], "--info") == 0) {
      std::cout << "FourierGrating " << FG_VERSION << "\n"
                << "cpu backend : " << fg::linalg::makeCpuBackend()->name() << "\n";
#ifdef FG_HAVE_CUDA
      std::cout << "cuda backend: " << fg::linalg::makeCudaBackend(0)->name() << "\n";
#else
      std::cout << "cuda backend: not built\n";
#endif
      return 0;
    }
    const auto cfg = fg::io::loadConfig(argv[1]);
    (void)cfg;
    return 0;
  } catch (const fg::NotImplemented& e) {
    std::cerr << "fgsolve: " << e.what() << "\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "fgsolve: error: " << e.what() << "\n";
    return 3;
  }
}
