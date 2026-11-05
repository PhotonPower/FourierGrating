#include "fg/fourier/sampling.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "fg/core/error.hpp"

namespace fg::fourier {

namespace {
bool isPow2(int n) { return n > 0 && (n & (n - 1)) == 0; }

/// In-place radix-2 FFT (forward, e^{-i ...}) of n = 2^k values with stride.
void fftPow2(cplx* a, int n, int stride) {
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i * stride], a[j * stride]);
  }
  for (int len = 2; len <= n; len <<= 1) {
    const double ang = -2.0 * PI / len;
    for (int k = 0; k < len / 2; ++k) {
      const cplx w(std::cos(ang * k), std::sin(ang * k));
      for (int i = k; i < n; i += len) {
        cplx& u = a[i * stride];
        cplx& v = a[(i + len / 2) * stride];
        const cplx t = v * w;
        v = u - t;
        u += t;
      }
    }
  }
}

/// Forward DFT of fixed length n: radix 2 for powers of two, otherwise Bluestein (chirp z) with a
/// power-of-two convolution of length L >= 2n - 1; O(n log n) for every n.
class FFTPlan {
 public:
  explicit FFTPlan(int n) : n_(n) {
    if (isPow2(n_)) return;
    L_ = 1;
    while (L_ < 2 * n_ - 1) L_ <<= 1;
    w_.resize(static_cast<std::size_t>(n_));
    for (int k = 0; k < n_; ++k) {  // w_k = exp(-i pi k^2 / n), k^2 reduced mod 2n (accuracy)
      const std::int64_t k2 = (static_cast<std::int64_t>(k) * k) % (2 * static_cast<std::int64_t>(n_));
      w_[static_cast<std::size_t>(k)] = std::exp(-I * (PI * static_cast<double>(k2) / n_));
    }
    b_.assign(static_cast<std::size_t>(L_), 0.0);
    b_[0] = std::conj(w_[0]);
    for (int k = 1; k < n_; ++k) b_[static_cast<std::size_t>(k)] = b_[static_cast<std::size_t>(L_ - k)] = std::conj(w_[static_cast<std::size_t>(k)]);
    fftPow2(b_.data(), L_, 1);
    buf_.resize(static_cast<std::size_t>(L_));
  }

  void forward(cplx* a, int stride) {
    if (L_ == 0) {
      fftPow2(a, n_, stride);
      return;
    }
    std::fill(buf_.begin(), buf_.end(), 0.0);
    for (int k = 0; k < n_; ++k) buf_[static_cast<std::size_t>(k)] = a[k * stride] * w_[static_cast<std::size_t>(k)];
    fftPow2(buf_.data(), L_, 1);
    for (int k = 0; k < L_; ++k) buf_[static_cast<std::size_t>(k)] = std::conj(buf_[static_cast<std::size_t>(k)] * b_[static_cast<std::size_t>(k)]);
    fftPow2(buf_.data(), L_, 1);  // inverse FFT via conjugation
    const double s = 1.0 / L_;
    for (int p = 0; p < n_; ++p) a[p * stride] = std::conj(buf_[static_cast<std::size_t>(p)]) * s * w_[static_cast<std::size_t>(p)];
  }

 private:
  int n_ = 0, L_ = 0;
  std::vector<cplx> w_, b_, buf_;
};
thread_local bool useFFT = false;

// ---- v0.10 paths (default, fft_sampling = false): direct DFT, radix 2 only for powers of two

std::vector<cplx> dftCoefficients1D(const std::function<cplx(double)>& f, double period,
                                        int pmax, int M) {
  FG_REQUIRE(M > 2 * pmax, "too few samples for the requested Fourier orders");
  std::vector<cplx> vals(static_cast<std::size_t>(M));
  for (int k = 0; k < M; ++k) vals[static_cast<std::size_t>(k)] = f(period * k / M);
  std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1));
  for (int p = -pmax; p <= pmax; ++p) {
    cplx s = 0.0;
    const cplx w = std::exp(-I * (2.0 * PI * p / M));
    cplx wk = 1.0;
    for (int k = 0; k < M; ++k) {
      s += vals[static_cast<std::size_t>(k)] * wk;
      wk *= w;
    }
    c[static_cast<std::size_t>(p + pmax)] = s / static_cast<double>(M);
  }
  return c;
}

CoefficientGrid dftCoefficients2D(const std::function<cplx(double, double)>& f, double dx,
                                      double dz, int Px, int Pz, int Mx, int Mz) {
  FG_REQUIRE(Mx > 2 * Px && Mz > 2 * Pz, "too few samples for the requested Fourier orders");
  // separable two-pass DFT: first along z for every x sample, then along x
  std::vector<std::vector<cplx>> rows(static_cast<std::size_t>(Mx));
  for (int i = 0; i < Mx; ++i) {
    const double x = dx * i / Mx;
    rows[static_cast<std::size_t>(i)] =
        dftCoefficients1D([&](double z) { return f(x, z); }, dz, Pz, Mz);
  }
  CoefficientGrid g(Px, Pz);
  for (int q = -Pz; q <= Pz; ++q) {
    const auto col = dftCoefficients1D(
        [&](double x) {
          const int i = static_cast<int>(std::lround(x / dx * Mx)) % Mx;
          return rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(q + Pz)];
        },
        dx, Px, Mx);
    for (int p = -Px; p <= Px; ++p) g(p, q) = col[static_cast<std::size_t>(p + Px)];
  }
  return g;
}

CoefficientGrid dftGridCoefficients2D(const std::vector<cplx>& f, int Mx, int Mz, int Px, int Pz) {
  FG_REQUIRE(Mx > 2 * Px && Mz > 2 * Pz, "too few samples for the requested Fourier orders");
  FG_REQUIRE(f.size() == static_cast<std::size_t>(Mx) * static_cast<std::size_t>(Mz), "sample array size");
  CoefficientGrid g(Px, Pz);
  const double norm = 1.0 / (static_cast<double>(Mx) * static_cast<double>(Mz));
  if (isPow2(Mx) && isPow2(Mz)) {
    std::vector<cplx> a(f);
    for (int i = 0; i < Mx; ++i) fftPow2(&a[static_cast<std::size_t>(i) * static_cast<std::size_t>(Mz)], Mz, 1);
    for (int k = 0; k < Mz; ++k) fftPow2(&a[static_cast<std::size_t>(k)], Mx, Mz);
    for (int p = -Px; p <= Px; ++p)
      for (int q = -Pz; q <= Pz; ++q)
        g(p, q) = a[static_cast<std::size_t>((p + Mx) % Mx) * static_cast<std::size_t>(Mz) +
                    static_cast<std::size_t>((q + Mz) % Mz)] * norm;
    return g;
  }
  const int nq = 2 * Pz + 1;
  std::vector<cplx> wz(static_cast<std::size_t>(Mz) * static_cast<std::size_t>(nq));
  for (int k = 0; k < Mz; ++k)
    for (int q = -Pz; q <= Pz; ++q)
      wz[static_cast<std::size_t>(k) * static_cast<std::size_t>(nq) + static_cast<std::size_t>(q + Pz)] =
          std::exp(-I * (2.0 * PI * static_cast<double>((static_cast<long>(q) * k) % Mz) / Mz));
  std::vector<cplx> rows(static_cast<std::size_t>(Mx) * static_cast<std::size_t>(nq), 0.0);
  for (int i = 0; i < Mx; ++i) {
    cplx* r = &rows[static_cast<std::size_t>(i) * static_cast<std::size_t>(nq)];
    const cplx* fi = &f[static_cast<std::size_t>(i) * static_cast<std::size_t>(Mz)];
    for (int k = 0; k < Mz; ++k) {
      const cplx* w = &wz[static_cast<std::size_t>(k) * static_cast<std::size_t>(nq)];
      for (int q = 0; q < nq; ++q) r[q] += fi[k] * w[q];
    }
  }
  for (int p = -Px; p <= Px; ++p)
    for (int q = 0; q < nq; ++q) {
      cplx s = 0.0;
      for (int i = 0; i < Mx; ++i)
        s += rows[static_cast<std::size_t>(i) * static_cast<std::size_t>(nq) + static_cast<std::size_t>(q)] *
             std::exp(-I * (2.0 * PI * static_cast<double>((static_cast<long>(p) * i) % Mx) / Mx));
      g(p, q - Pz) = s * norm;
    }
  return g;
}

}  // namespace

void setFFTSampling(bool on) { useFFT = on; }
bool fftSampling() { return useFFT; }

std::vector<cplx> fft(std::vector<cplx> a) {
  FFTPlan(static_cast<int>(a.size())).forward(a.data(), 1);
  return a;
}

std::vector<cplx> sampledCoefficients1D(const std::function<cplx(double)>& f, double period,
                                        int pmax, int M) {
  if (!useFFT) return dftCoefficients1D(f, period, pmax, M);
  FG_REQUIRE(M > 2 * pmax, "too few samples for the requested Fourier orders");
  std::vector<cplx> vals(static_cast<std::size_t>(M));
  for (int k = 0; k < M; ++k) vals[static_cast<std::size_t>(k)] = f(period * k / M);
  vals = fft(std::move(vals));
  std::vector<cplx> c(static_cast<std::size_t>(2 * pmax + 1));
  for (int p = -pmax; p <= pmax; ++p) c[static_cast<std::size_t>(p + pmax)] = vals[static_cast<std::size_t>((p + M) % M)] / static_cast<double>(M);
  return c;
}

std::vector<std::vector<cplx>> piecewiseSmoothCoefficients(
    const std::vector<double>& breaks, int K, const std::function<void(double, int, cplx*)>& f,
    double period, int pmax) {
  static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                               -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                               0.7966664774136267,  0.9602898564975363};
  static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                               0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  FG_REQUIRE(breaks.size() >= 2, "need at least one interval");
  std::vector<std::vector<cplx>> c(static_cast<std::size_t>(K), std::vector<cplx>(static_cast<std::size_t>(2 * pmax + 1), 0.0));
  const double Kg = 2.0 * PI / period;
  const double hmax = 1.0 / (Kg * std::max(pmax, 1));
  std::vector<cplx> val(static_cast<std::size_t>(K));
  for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
    const double a = breaks[k], b = breaks[k + 1];
    if (b <= a) continue;
    const int m = std::max(1, static_cast<int>(std::ceil((b - a) / hmax)));
    const double h = (b - a) / m;
    for (int s = 0; s < m; ++s) {
      const double xc = a + (s + 0.5) * h;
      for (int q = 0; q < 8; ++q) {
        const double x = xc + 0.5 * h * gx[q];
        const double w = 0.5 * h * gw[q] / period;
        f(x, static_cast<int>(k), val.data());
        const cplx e1 = std::exp(-I * (Kg * x));
        cplx ep = std::pow(e1, -pmax);
        for (int p = -pmax; p <= pmax; ++p) {
          for (int i = 0; i < K; ++i) c[static_cast<std::size_t>(i)][static_cast<std::size_t>(p + pmax)] += w * val[static_cast<std::size_t>(i)] * ep;
          ep *= e1;
        }
      }
    }
  }
  return c;
}

CoefficientGrid sampledCoefficients2D(const std::function<cplx(double, double)>& f, double dx,
                                      double dz, int Px, int Pz, int Mx, int Mz) {
  if (!useFFT) return dftCoefficients2D(f, dx, dz, Px, Pz, Mx, Mz);
  std::vector<cplx> v(static_cast<std::size_t>(Mx) * static_cast<std::size_t>(Mz));
  for (int i = 0; i < Mx; ++i)
    for (int k = 0; k < Mz; ++k) v[static_cast<std::size_t>(i) * static_cast<std::size_t>(Mz) + static_cast<std::size_t>(k)] = f(dx * i / Mx, dz * k / Mz);
  return gridCoefficients2D(v, Mx, Mz, Px, Pz);
}

CoefficientGrid gridCoefficients2D(const std::vector<double>& f, int Mx, int Mz, int Px, int Pz) {
  return gridCoefficients2D(std::vector<cplx>(f.begin(), f.end()), Mx, Mz, Px, Pz);
}

CoefficientGrid gridCoefficients2D(const std::vector<cplx>& f, int Mx, int Mz, int Px, int Pz) {
  if (!useFFT) return dftGridCoefficients2D(f, Mx, Mz, Px, Pz);
  FG_REQUIRE(Mx > 2 * Px && Mz > 2 * Pz, "too few samples for the requested Fourier orders");
  FG_REQUIRE(f.size() == static_cast<std::size_t>(Mx) * static_cast<std::size_t>(Mz), "sample array size");
  CoefficientGrid g(Px, Pz);
  const double norm = 1.0 / (static_cast<double>(Mx) * static_cast<double>(Mz));
  // separable FFT: all rows along z, then only the 2 Pz + 1 needed columns along x
  std::vector<cplx> a(f);
  FFTPlan pz(Mz), px(Mx);
  for (int i = 0; i < Mx; ++i) pz.forward(&a[static_cast<std::size_t>(i) * static_cast<std::size_t>(Mz)], 1);
  for (int q = -Pz; q <= Pz; ++q) {
    const int k = (q + Mz) % Mz;
    px.forward(&a[static_cast<std::size_t>(k)], Mz);
    for (int p = -Px; p <= Px; ++p)
      g(p, q) = a[static_cast<std::size_t>((p + Mx) % Mx) * static_cast<std::size_t>(Mz) + static_cast<std::size_t>(k)] * norm;
  }
  return g;
}

CoefficientGrid convolveCoefficients(const CoefficientGrid& a, const CoefficientGrid& b, int Px, int Pz) {
  FG_REQUIRE(b.Px() >= Px + a.Px() && b.Pz() >= Pz + a.Pz(), "convolveCoefficients: b too small");
  CoefficientGrid c(Px, Pz);
  for (int pp = -a.Px(); pp <= a.Px(); ++pp)
    for (int qq = -a.Pz(); qq <= a.Pz(); ++qq) {
      const cplx av = a(pp, qq);
      if (av == 0.0) continue;
      for (int p = -Px; p <= Px; ++p)
        for (int q = -Pz; q <= Pz; ++q) c(p, q) += av * b(p - pp, q - qq);
    }
  return c;
}

}  // namespace fg::fourier
