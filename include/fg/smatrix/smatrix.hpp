// FourierGrating - S-matrix algebra (Plan section "S-Matrix-Formalismus").
// Convention eq:Sdef:  (b_top_out ; a_bottom_out) = [S11 S12; S21 S22] (b_bottom_in ; a_top_in)
// a = down-going amplitudes, b = up-going amplitudes, each block n_h x n_h.
#pragma once

#include "fg/core/types.hpp"

namespace fg {

struct SMatrix {
  MatC S11, S12, S21, S22;

  int nh() const { return static_cast<int>(S11.rows()); }
  /// Neutral element of the star product (transparent layer of zero thickness).
  static SMatrix identity(int nh);
  /// Full 2nh x 2nh matrix [S11 S12; S21 S22].
  MatC full() const;
};

/// S-matrix from a transfer matrix T: (a_top; b_top) = T (a_bottom; b_bottom) (Plan eq:T2S).
/// Use only for thin slices (T11 must be well conditioned).
SMatrix fromTransfer(const MatC& T);

/// Redheffer star product C = lower * upper (Plan eq:star).
SMatrix star(const SMatrix& lower, const SMatrix& upper);

/// Interior amplitudes at the interface between a lower stack S_L and an upper stack S_U
/// (Plan eq:internal). Returns (a_q ; b_q).
VecC interfaceAmplitudes(const SMatrix& lower, const SMatrix& upper, const VecC& b_in_bottom,
                         const VecC& a_in_top);

}  // namespace fg
