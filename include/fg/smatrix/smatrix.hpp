// FourierGrating - S-matrix algebra (Plan section "S-Matrix-Formalismus").
// Convention eq:Sdef:  (b_top_out ; a_bottom_out) = [S11 S12; S21 S22] (b_bottom_in ; a_top_in)
// a = down-going amplitudes, b = up-going amplitudes, each block n_h x n_h.
#pragma once

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include "fg/core/error.hpp"
#include "fg/core/parallel.hpp"
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

/// Star product of S_0 * S_1 * ... * S_{K-1} as a balanced binary tree of fixed shape (Plan section
/// "M15", step 3): pairs (0,1), (2,3), ..., an unpaired last element moves up one level, then
/// recursively. The products of one level run on up to `threads` threads. Reference for StarChain.
/// S: SMatrix, or the device S-matrices of the resident chain (step 5) with their star product.
template <class S>
S starTree(std::vector<S> v, int threads, const std::function<S(const S&, const S&)>& product) {
  FG_REQUIRE(!v.empty(), "star product of an empty chain");
  while (v.size() > 1) {
    const int pairs = static_cast<int>(v.size() / 2);
    std::vector<S> next((v.size() + 1) / 2);
    parallelFor(pairs, threads, [&](int i) {
      const auto k = static_cast<std::size_t>(i);
      next[k] = product(v[2 * k], v[2 * k + 1]);
    });
    if (v.size() % 2) next.back() = std::move(v.back());
    v = std::move(next);
  }
  return std::move(v.front());
}
inline SMatrix starTree(std::vector<SMatrix> v, int threads = 1) {
  return starTree<SMatrix>(std::move(v), threads, [](const SMatrix& a, const SMatrix& b) { return star(a, b); });
}

/// Ordered star-product chain fed element by element (bottom to top).
/// Sequential: ((S_0 * S_1) * S_2) * ... (the v0.10 order, bit-identical).
/// Tree: the shape of starTree for any `threads` and memory O(threads + log K): aligned blocks of
/// 2^m >= threads elements are reduced in parallel, complete subtrees are merged like a binary
/// counter, the rest is folded from the right. The result is bit-identical for every thread count.
template <class S>
class BasicStarChain {
 public:
  using Product = std::function<S(const S&, const S&)>;
  BasicStarChain(bool tree, int threads, Product product)
      : tree_(tree), threads_(std::max(1, threads)), product_(std::move(product)) {
    while (block_ < threads_) {
      block_ *= 2;
      ++blockLevel_;
    }
  }
  void push(S x) {
    ++count_;
    if (!tree_) {
      total_ = count_ == 1 ? std::move(x) : product_(total_, x);
      return;
    }
    buf_.push_back(std::move(x));
    if (static_cast<int>(buf_.size()) == block_) {
      S b = starTree<S>(std::move(buf_), threads_, product_);
      buf_.clear();
      pushLevel(std::move(b), blockLevel_);
    }
  }
  /// The product of all pushed elements (at least one).
  S result() {
    FG_REQUIRE(count_ > 0, "star product of an empty chain");
    count_ = 0;
    if (!tree_) return std::move(total_);
    // the complete subtrees have strictly decreasing levels; the partial last block is the
    // rightmost subtree of the fixed-shape tree, so the remaining products fold from the right
    std::vector<S> parts;
    for (auto& [l, x] : stack_) parts.push_back(std::move(x));
    stack_.clear();
    if (!buf_.empty()) parts.push_back(starTree<S>(std::move(buf_), threads_, product_));
    buf_.clear();
    S acc = std::move(parts.back());
    for (std::size_t k = parts.size() - 1; k-- > 0;) acc = product_(parts[k], acc);
    return acc;
  }
  bool empty() const { return count_ == 0; }

 private:
  void pushLevel(S x, int level) {
    while (!stack_.empty() && stack_.back().first == level) {
      x = product_(stack_.back().second, x);
      stack_.pop_back();
      ++level;
    }
    stack_.emplace_back(level, std::move(x));
  }
  bool tree_;
  int threads_, block_ = 1, blockLevel_ = 0;
  Product product_;
  long count_ = 0;
  S total_{};                             ///< sequential
  std::vector<std::pair<int, S>> stack_;  ///< tree: complete subtrees (level, product)
  std::vector<S> buf_;                    ///< tree: current aligned block
};

class StarChain : public BasicStarChain<SMatrix> {
 public:
  explicit StarChain(bool tree = false, int threads = 1)
      : BasicStarChain<SMatrix>(tree, threads, [](const SMatrix& a, const SMatrix& b) { return star(a, b); }) {}
};

/// Interior amplitudes at the interface between a lower stack S_L and an upper stack S_U
/// (Plan eq:internal). Returns (a_q ; b_q).
VecC interfaceAmplitudes(const SMatrix& lower, const SMatrix& upper, const VecC& b_in_bottom,
                         const VecC& a_in_top);

/// Element with internal sources (Plan section "Streufeldformulierung", Verkettung):
/// (b_top_out ; a_bottom_out) = S (b_bottom_in ; a_top_in) + (u ; v), u, v: the amplitudes the sources
/// send out at the top and the bottom when nothing comes in.
struct SourcedSMatrix {
  SMatrix S;
  VecC u, v;
  static SourcedSMatrix identity(int nh);
};

/// Extended star product: the sources of both parts act like additional incoming amplitudes at the
/// common interface, u = u_U + U11 (1 - L12 U21)^-1 (u_L + L12 v_U), v = v_L + L22 (1 - U21 L12)^-1 (v_U + U21 u_L).
SourcedSMatrix star(const SourcedSMatrix& lower, const SourcedSMatrix& upper);

/// Interior amplitudes (a_q ; b_q) between two parts with sources.
VecC interfaceAmplitudes(const SourcedSMatrix& lower, const SourcedSMatrix& upper, const VecC& b_in_bottom,
                         const VecC& a_in_top);

}  // namespace fg
