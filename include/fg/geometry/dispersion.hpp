// FourierGrating - wavelength-dependent material values (Plan section "Dispersive Materialien (M22)").
#pragma once

#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "fg/core/types.hpp"

namespace fg {

/// One material value as a function of the vacuum wavelength (micrometres). Every model has a natural
/// quantity (Plan sec:disp_models): Sellmeier and Drude-Lorentz give eps, Cauchy gives n; a table gives
/// what the key it stands under says (n, k under "n"; Re eps, Im eps under "eps" and "eps_principal").
class Dispersion {
 public:
  enum class Quantity { N, Eps };

  /// Plan eq:disp_sellmeier: eps = A + sum B_k lambda^2 / (lambda^2 - C_k), C_k in um^2.
  struct Sellmeier {
    double A = 1.0;
    std::vector<double> B, C_um2;
    /// optional validity range (um), e.g. of a formula imported from refractiveindex.info; outside: error
    double rangeLo_um = 0.0, rangeHi_um = std::numeric_limits<double>::infinity();
  };
  /// Plan eq:disp_cauchy: n = A + B / lambda^2 + C / lambda^4 + i k.
  struct Cauchy {
    double A = 1.0, B_um2 = 0.0, C_um4 = 0.0, k = 0.0;
  };
  /// Plan eq:disp_dl (Rakic form, energies in eV).
  struct DrudeLorentz {
    struct Oscillator {
      double f = 0, omega_eV = 0, gamma_eV = 0;
    };
    double epsInf = 1.0, omegaP_eV = 0.0, f0 = 0.0, gamma0_eV = 0.0;
    std::vector<Oscillator> oscillators;
  };
  /// Nodes (lambda_i, a_i, b_i), lambda strictly increasing after unit conversion; each column interpolated
  /// for itself (linear or not-a-knot cubic spline).
  struct Table {
    std::vector<double> lambda_um, a, b;
    std::vector<int> row;   ///< node i came from row[i] of the input (eV tables are reversed)
    bool cubic = false;
    bool extrapolateConstant = false;
    bool inlineData = false;  ///< "data" (values addressable by JSON pointer) instead of "file"
    std::vector<double> Ma, Mb;  ///< spline second derivatives (cubic)
  };

  Dispersion(std::string pointer, Quantity q, std::variant<Sellmeier, Cauchy, DrudeLorentz, Table> model);

  /// JSON pointer of the value in the configuration (e.g. "/materials/au/n"), used in every message.
  const std::string& pointer() const { return pointer_; }
  Quantity quantity() const { return quantity_; }
  const std::variant<Sellmeier, Cauchy, DrudeLorentz, Table>& model() const { return model_; }

  /// Natural quantity (n or eps) at lambda; throws fg::Error (with the pointer) outside the valid range or
  /// at a Sellmeier pole.
  cplx value(double lambda_um) const;
  /// eps = n^2 resp. the value.
  cplx eps(double lambda_um) const;
  cplx dValueDLambda(double lambda_um) const;
  cplx dEpsDLambda(double lambda_um) const;
  /// Derivative of the natural quantity with respect to a model number, given by the JSON pointer tokens
  /// below the value (e.g. {"sellmeier", "B", "0"}); throws for numbers that are not parameters (table
  /// wavelengths) and unknown pointers.
  cplx dValueDParam(const std::vector<std::string>& tokens, double lambda_um) const;

  /// Valid wavelength range [lo, hi] (tables; Sellmeier with "range_um"; otherwise (0, inf)).
  std::pair<double, double> validRange() const;
  /// Throws if lambda cannot be evaluated (range, pole), the message naming pointer, lambda and range.
  void check(double lambda_um) const;
  /// True if lambda lies outside a table that extrapolates constantly (warning in the diagnostics).
  bool extrapolated(double lambda_um) const;
  /// True if lambda is a node of a linear table (d eps / d lambda one-sided there).
  bool atLinearNode(double lambda_um) const;

 private:
  std::string pointer_;
  Quantity quantity_;
  std::variant<Sellmeier, Cauchy, DrudeLorentz, Table> model_;
};

/// h c / e in eV um (Plan sec:disp_models): hbar omega [eV] = 1.239841984 / lambda [um].
inline constexpr double HC_EV_UM = 1.239841984;

/// Not-a-knot cubic spline (second derivatives at the nodes) through (x_i, y_i), at least 4 nodes.
std::vector<double> splineSecondDerivatives(const std::vector<double>& x, const std::vector<double>& y);

}  // namespace fg
