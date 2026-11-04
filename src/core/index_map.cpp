#include "fg/core/index_map.hpp"

#include "fg/core/error.hpp"

namespace fg {

IndexMap::IndexMap(int Nx, int Nz) : Nx_(Nx), Nz_(Nz) {
  FG_REQUIRE(Nx >= 0 && Nz >= 0, "truncation orders must be non-negative");
}

int IndexMap::index(int n, int m) const {
  FG_REQUIRE(n >= -Nx_ && n <= Nx_ && m >= -Nz_ && m <= Nz_, "order out of truncation range");
  return (m + Nz_) + (2 * Nz_ + 1) * (n + Nx_);
}

std::pair<int, int> IndexMap::orders(int r) const {
  FG_REQUIRE(r >= 0 && r < size(), "single index out of range");
  const int n = -Nx_ + r / (2 * Nz_ + 1);
  const int m = r - (2 * Nz_ + 1) * (n + Nx_) - Nz_;
  return {n, m};
}

}  // namespace fg
