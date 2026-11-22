#include "cyl_basis.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/cyl/angular.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::detail {

namespace {

cplx isotropicEps(const Stack& stack, const std::string& name) {
  const Material& m = stack.material(name);
  FG_REQUIRE(m.isotropic, "cylindrical: material \"" + name + "\" must be isotropic (M20)");
  return m.eps;
}

/// Region scaled by k0 (normalized coordinates) and normalized.
geometry::Region scaledRegion(const geometry::Region& r, double k0) {
  geometry::Region out = r;
  auto scale = [k0](geometry::Contour2D& c) {
    for (auto& seg : c.segments)
      for (auto& p : seg.P) p = {p[0] * k0, p[1] * k0};
  };
  scale(out.outer);
  for (auto& h : out.holes) scale(h);
  out.normalize(1e-9 * k0);
  return out;
}

/// Block-diagonal factorization matrices of a class from those of its harmonics (frak A: (+, -) blocks of each).
cyl::CylQ blockDiagonal(const std::vector<cyl::CylQ>& q, Eigen::Index n) {
  const auto nb = static_cast<Eigen::Index>(q.size()), nf = nb * n;
  cyl::CylQ Q;
  for (MatC* M : {&Q.Qpp, &Q.Qpm, &Q.Qmp, &Q.Qmm, &Q.Qyy}) *M = MatC::Zero(nf, nf);
  const bool li = q.front().A.size() > 0;
  if (li) Q.A = MatC::Zero(2 * nf, 2 * nf);
  for (Eigen::Index b = 0; b < nb; ++b) {
    const cyl::CylQ& s = q[static_cast<std::size_t>(b)];
    Q.Qpp.block(b * n, b * n, n, n) = s.Qpp;
    Q.Qpm.block(b * n, b * n, n, n) = s.Qpm;
    Q.Qmp.block(b * n, b * n, n, n) = s.Qmp;
    Q.Qmm.block(b * n, b * n, n, n) = s.Qmm;
    Q.Qyy.block(b * n, b * n, n, n) = s.Qyy;
    if (li)
      for (Eigen::Index u = 0; u < 2; ++u)
        for (Eigen::Index v = 0; v < 2; ++v) Q.A.block(u * nf + b * n, v * nf + b * n, n, n) = s.A.block(u * n, v * n, n, n);
  }
  return Q;
}

/// The block of harmonic b of a block-diagonal class matrix (the inverse of blockDiagonal).
cyl::CylQ diagonalBlock(const cyl::CylQ& Q, Eigen::Index b, Eigen::Index n) {
  const Eigen::Index nf = Q.Qpp.rows();
  cyl::CylQ s;
  s.Qpp = Q.Qpp.block(b * n, b * n, n, n);
  s.Qpm = Q.Qpm.block(b * n, b * n, n, n);
  s.Qmp = Q.Qmp.block(b * n, b * n, n, n);
  s.Qmm = Q.Qmm.block(b * n, b * n, n, n);
  s.Qyy = Q.Qyy.block(b * n, b * n, n, n);
  if (Q.A.size()) {
    s.A.resize(2 * n, 2 * n);
    for (Eigen::Index u = 0; u < 2; ++u)
      for (Eigen::Index v = 0; v < 2; ++v) s.A.block(u * n, v * n, n, n) = Q.A.block(u * nf + b * n, v * nf + b * n, n, n);
  }
  return s;
}

}  // namespace

std::vector<std::unique_ptr<CylLayer>> buildCylLayers(const Stack& stack, double k0) {
  std::vector<std::unique_ptr<CylLayer>> out;
  for (const LayerEntry& le : stack.layers) {
    auto L = std::make_unique<CylLayer>();
    if (const auto* h = std::get_if<HomogeneousLayer>(&le.spec)) {
      L->thickness_ = h->thickness_um;
      L->eps = isotropicEps(stack, h->material);
      L->epsList = {L->eps};
    } else if (const auto* r = std::get_if<RingsLayer>(&le.spec)) {
      FG_REQUIRE(!r->radii_um.empty() && r->materials.size() == r->radii_um.size() + 1,
                 "rings: n radii and n + 1 materials (the last one outside)");
      L->thickness_ = r->thickness_um;
      L->homogeneous_ = false;
      for (std::size_t i = 0; i < r->radii_um.size(); ++i) {
        FG_REQUIRE(r->radii_um[i] > (i ? r->radii_um[i - 1] : 0.0),
                   "rings: radii increasing and > 0");
        L->rings.radii.push_back(r->radii_um[i] * k0);
        L->rings.eps.push_back(isotropicEps(stack, r->materials[i]));
      }
      L->rings.epsInf = isotropicEps(stack, r->materials.back());
      L->epsList = L->rings.eps;
      L->epsList.push_back(L->rings.epsInf);
    } else if (std::holds_alternative<Polygon2DLayer>(le.spec) || std::holds_alternative<Contour2DLayer>(le.spec)) {
      const int P = stack.cylindrical.rotation;
      FG_REQUIRE(P >= 1, "cylindrical: polygon2d and contour2d layers need \"symmetry\": {\"rotation\": P} (M21)");
      std::vector<geometry::Region> regions;
      std::string background;
      if (const auto* p = std::get_if<Polygon2DLayer>(&le.spec)) {
        L->thickness_ = p->thickness_um;
        background = p->background;
        for (const auto& inc : p->inclusions) {
          geometry::Region reg;
          reg.outer = geometry::Contour2D::polygon(inc.vertices);
          reg.material = inc.material;
          regions.push_back(std::move(reg));
        }
      } else {
        const auto& c = std::get<Contour2DLayer>(le.spec);
        L->thickness_ = c.thickness_um;
        background = c.background;
        regions = c.regions;
      }
      L->homogeneous_ = false;
      L->angular = true;
      double rmax = 0;
      for (const auto& reg : regions) {
        L->section.regions.push_back(scaledRegion(reg, k0));
        L->section.eps.push_back(isotropicEps(stack, reg.material));
        rmax = std::max(rmax, cyl::criticalRadii(L->section.regions.back()).back());
      }
      L->section.epsBg = isotropicEps(stack, background);
      // Plan sec:cylP: the cross section must go over into itself under a rotation by 2 pi / P
      if (P > 1 && !cyl::hasRotationSymmetry(L->section.regions, P, 1e-9 * rmax)) {
        const int found = cyl::rotationOrder(L->section.regions, std::max(P, 12), 1e-9 * rmax);
        throw Error("cylindrical: the cross section of a " + std::string(std::holds_alternative<Polygon2DLayer>(le.spec) ? "polygon2d" : "contour2d") +
                    " layer is not " + std::to_string(P) + "-fold rotationally symmetric about the axis (largest order found: " +
                    std::to_string(found) + ")");
      }
      L->epsList = L->section.eps;
      L->epsList.push_back(L->section.epsBg);
    } else {
      throw Error("cylindrical: only homogeneous, rings (M20), polygon2d and contour2d (M21) layers");
    }
    out.push_back(std::move(L));
  }
  return out;
}

Stack expandRingsProfiles(const Stack& stack, int defaultSlices) {
  Stack out = stack;
  out.layers.clear();
  for (const LayerEntry& le : stack.layers) {
    const auto* rp = std::get_if<RingsProfileLayer>(&le.spec);
    if (!rp) {
      out.layers.push_back(le);
      continue;
    }
    for (LayerEntry& e : ringsStaircase(*rp, defaultSlices)) {
      e.factorization = le.factorization;
      out.layers.push_back(std::move(e));
    }
  }
  return out;
}

CylBasis::CylBasis(const cyl::SpectralGrid& grid, int m, cplx epsSup, cplx epsSub,
                   const std::vector<std::unique_ptr<CylLayer>>& layers, cyl::CylFactorization f)
    : CylBasis(grid, cyl::HarmonicMap{m, 1, 0}, epsSup, epsSub, layers, f) {}

CylBasis::CylBasis(const cyl::SpectralGrid& grid, const cyl::HarmonicMap& map, cplx epsSup, cplx epsSub,
                   const std::vector<std::unique_ptr<CylLayer>>& layers, cyl::CylFactorization f)
    : grid_(grid),
      map_(map),
      k_(cyl::tiledK(grid.k, map.blocks())),
      ref_(cyl::homogeneousModes(epsSup, k_)),
      sub_(cyl::homogeneousModes(epsSub, k_)),
      f_(f) {
  const int nb = map.blocks();
  for (const auto& L : layers) {
    layers_.push_back(L.get());
    if (L->homogeneous()) {
      M_.emplace_back();
      Q_.emplace_back();
      Qyy_.emplace_back();
      continue;
    }
    if (L->angular) {
      Q_.push_back(cyl::coupledQ(grid, map, L->section));
    } else if (nb == 1) {
      Q_.push_back(cyl::cylQ(grid, map.m0, L->rings, f));
    } else {
      std::vector<cyl::CylQ> q;
      for (int b = 0; b < nb; ++b) q.push_back(cyl::cylQ(grid, map.m(b), L->rings, f));
      Q_.push_back(blockDiagonal(q, grid.k.size()));
    }
    M_.push_back(cyl::propagationMatrix(k_, Q_.back()));
    Qyy_.emplace_back(Q_.back().Qyy);
  }
}

std::size_t CylBasis::index(const StackLayer& layer) const {
  const auto it = std::find(layers_.begin(), layers_.end(), &layer);
  FG_REQUIRE(it != layers_.end(), "CylBasis: unknown layer");
  return static_cast<std::size_t>(it - layers_.begin());
}

void CylBasis::setPlaneWave(double kpar, double Psi) {
  src_.clear();
  const int nb = map_.blocks();
  const Eigen::Index n = grid_.k.size(), nf = nb * n;
  for (const StackLayer* l : layers_) {
    const auto& L = static_cast<const CylLayer&>(*l);
    const std::size_t i = src_.size();
    if (L.homogeneous()) {
      src_.emplace_back();
    } else if (L.angular) {
      src_.push_back(cyl::coupledSource(grid_, map_, L.section, kpar, Psi));
    } else if (nb == 1) {
      src_.push_back(cyl::planeWaveSource(grid_, map_.m0, L.rings, kpar, Psi, f_, Q_[i]));
    } else {  // ring layer in a class: the harmonics do not couple
      cyl::PlaneWaveSource B;
      B.By = VecC::Zero(nf);
      B.Bpm = MatC::Zero(2 * nf, 2);
      if (f_ == cyl::CylFactorization::LiR) B.Apm = MatC::Zero(2 * nf, 2);
      for (int b = 0; b < nb; ++b) {
        const cyl::PlaneWaveSource s =
            cyl::planeWaveSource(grid_, map_.m(b), L.rings, kpar, Psi, f_, diagonalBlock(Q_[i], b, n));
        B.By.segment(b * n, n) = s.By;
        for (Eigen::Index u = 0; u < 2; ++u) {
          B.Bpm.block(u * nf + b * n, 0, n, 2) = s.Bpm.middleRows(u * n, n);
          if (s.Apm.size()) B.Apm.block(u * nf + b * n, 0, n, 2) = s.Apm.middleRows(u * n, n);
        }
      }
      src_.push_back(std::move(B));
    }
  }
}

VecC CylBasis::source(const StackLayer& layer, cplx cp, cplx cm, cplx ey) const {
  const std::size_t i = index(layer);
  FG_REQUIRE(i < src_.size(), "CylBasis: setPlaneWave first");
  if (layer.homogeneous()) return VecC::Zero(4 * k_.size());
  return cyl::scatteredSource(k_, Qyy_[i], src_[i], cp, cm, ey);
}

VecC CylBasis::ey(std::size_t i, const VecC& hp, const VecC& hm, cplx eyBg) const {
  const VecC dy = -0.5 * k_.cwiseProduct(hp + hm);
  const auto& L = static_cast<const CylLayer&>(*layers_.at(i));
  if (L.homogeneous()) return dy / L.eps;
  if (eyBg == 0.0 || src_.empty()) return Qyy_[i].solve(dy);
  return Qyy_[i].solve(VecC(dy - eyBg * src_[i].By));
}

VecC CylBasis::g(std::size_t i, const VecC& epm, cplx cp, cplx cm) const {
  if (f_ != cyl::CylFactorization::LiR || layers_.at(i)->homogeneous() || Q_[i].A.size() == 0) return VecC();
  if (Alu_.size() != layers_.size()) Alu_.resize(layers_.size());
  if (!Alu_[i]) Alu_[i] = std::make_unique<Eigen::PartialPivLU<MatC>>(Q_[i].A);
  if ((cp == 0.0 && cm == 0.0) || src_.empty()) return Alu_[i]->solve(epm);
  return Alu_[i]->solve(VecC(epm - src_[i].Apm * Eigen::Vector2cd(cp, cm)));
}

MatC CylBasis::M(const StackLayer& layer, double /*y_um*/) const {
  const std::size_t i = index(layer);
  const auto& L = static_cast<const CylLayer&>(layer);
  // homogeneous layers: only for RK (the assembly uses the analytic modes)
  if (L.homogeneous()) {
    const Eigen::Index n = k_.size();
    const MatC Id = MatC::Identity(n, n);
    return cyl::propagationMatrix(
        k_, {L.eps * Id, MatC::Zero(n, n), MatC::Zero(n, n), L.eps * Id, L.eps * Id, MatC()});
  }
  return M_[i];
}

// M is block anti-diagonal: reduced path (Plan eq:rcwa_reduced)
LayerModes CylBasis::modes(const MatC& M) const {
  return layerModesAuto(M);
}

LayerModes CylBasis::homogeneousModes(const StackLayer& layer) const {
  return layerModesHomogeneous(cyl::homogeneousModes(static_cast<const CylLayer&>(layer).eps, k_));
}

double CylBasis::rhoEstimate(const StackLayer& layer) const {
  double rho = 0;
  for (Eigen::Index j = 0; j < grid_.k.size(); ++j)
    for (const cplx e : static_cast<const CylLayer&>(layer).epsList)
      rho = std::max(rho, std::abs(std::sqrt(grid_.k(j) * grid_.k(j) - e)));
  return rho;
}

}  // namespace fg::detail
