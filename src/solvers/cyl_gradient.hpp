// FourierGrating - gradients of cylindrical structures (M20, V21i).
#pragma once

#include <string>
#include <vector>

#include "fg/adjoint/gradient.hpp"
#include "fg/io/config.hpp"

namespace fg::detail {

/// Discrete adjoint of a Gaussian beam through a cylindrical stack (total field, RCWA): objective on R and T
/// (Objective::cylindrical); parameters "/layers/k/thickness_um", "/layers/k/radii_um/i" (rings layers),
/// "/materials/name/n|eps[/0|1]" (not the half spaces). The spectral grid is fixed.
GradientResult cylindricalGradient(const std::string& configJson, const io::RunConfig& cfg, const Wavelength& wl,
                                   const Objective& objective, const std::vector<std::string>& params);

}  // namespace fg::detail
