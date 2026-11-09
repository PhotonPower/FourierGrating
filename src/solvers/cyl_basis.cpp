#include "cyl_basis.hpp"

#include <algorithm>
#include <cmath>

#include "fg/core/error.hpp"
#include "fg/cyl/propagation.hpp"
#include "fg/solvers/rcwa.hpp"

namespace fg::detail {

namespace {

cplx isotropicEps(const Stack& stack, const std::string& name) {
  const Material& m = stack.material(name);
  FG_REQUIRE(m.isotropic, "cylindrical: material \"" + name + "\" must be isotropic (M20)");
  return m.eps;
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
    } else {
      throw Error("cylindrical: only homogeneous and rings layers (M20)");
    }
    out.push_back(std::move(L));
  }
  return out;
}

CylBasis::CylBasis(const cyl::SpectralGrid& grid, int m, cplx epsSup, cplx epsSub,
                   const std::vector<std::unique_ptr<CylLayer>>& layers, cyl::CylFactorization f)
    : grid_(grid),
      ref_(cyl::homogeneousModes(epsSup, grid.k)),
      sub_(cyl::homogeneousModes(epsSub, grid.k)) {
  for (const auto& L : layers) {
    layers_.push_back(L.get());
    M_.push_back(L->homogeneous()
                     ? MatC()
                     : cyl::propagationMatrix(grid.k, cyl::cylQ(grid, m, L->rings, f)));
  }
}

MatC CylBasis::M(const StackLayer& layer, double /*y_um*/) const {
  const auto it = std::find(layers_.begin(), layers_.end(), &layer);
  FG_REQUIRE(it != layers_.end(), "CylBasis: unknown layer");
  const auto& L = static_cast<const CylLayer&>(layer);
  // homogeneous layers: only for RK (the assembly uses the analytic modes)
  if (L.homogeneous()) {
    const Eigen::Index n = grid_.k.size();
    const MatC Id = MatC::Identity(n, n);
    return cyl::propagationMatrix(
        grid_.k, {L.eps * Id, MatC::Zero(n, n), MatC::Zero(n, n), L.eps * Id, L.eps * Id});
  }
  return M_[static_cast<std::size_t>(it - layers_.begin())];
}

// M is block anti-diagonal: reduced path (Plan eq:rcwa_reduced)
LayerModes CylBasis::modes(const MatC& M) const {
  return layerModesAuto(M);
}

LayerModes CylBasis::homogeneousModes(const StackLayer& layer) const {
  return layerModesHomogeneous(
      cyl::homogeneousModes(static_cast<const CylLayer&>(layer).eps, grid_.k));
}

double CylBasis::rhoEstimate(const StackLayer& layer) const {
  double rho = 0;
  for (Eigen::Index j = 0; j < grid_.k.size(); ++j)
    for (const cplx e : static_cast<const CylLayer&>(layer).epsList)
      rho = std::max(rho, std::abs(std::sqrt(grid_.k(j) * grid_.k(j) - e)));
  return rho;
}

}  // namespace fg::detail
