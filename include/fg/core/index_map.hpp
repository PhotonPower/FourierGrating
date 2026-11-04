// FourierGrating - mapping between Fourier orders (n, m) and the single index r.
// Plan eq:indexmap (corrects book eq. IX.14, which divides by N+1 instead of 2N+1).
#pragma once

#include <utility>

namespace fg {

/// Rectangular truncation |n| <= Nx, |m| <= Nz. Nz == 0 means a 1D grating.
/// r = (m + Nz) + (2 Nz + 1) (n + Nx), 0-based.
class IndexMap {
 public:
  IndexMap(int Nx, int Nz);

  int Nx() const { return Nx_; }
  int Nz() const { return Nz_; }
  bool is2D() const { return Nz_ > 0; }
  /// Number of retained orders n_f.
  int size() const { return (2 * Nx_ + 1) * (2 * Nz_ + 1); }

  /// Single index of order (n, m); throws fg::Error if out of range.
  int index(int n, int m) const;
  /// Orders (n, m) of single index r; throws fg::Error if out of range.
  std::pair<int, int> orders(int r) const;
  int n(int r) const { return orders(r).first; }
  int m(int r) const { return orders(r).second; }
  /// Index of the specular order (0, 0).
  int center() const { return index(0, 0); }

 private:
  int Nx_;
  int Nz_;
};

}  // namespace fg
