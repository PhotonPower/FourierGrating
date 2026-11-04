#include "fg/solvers/solver.hpp"

#include "fg/core/error.hpp"

namespace fg {

namespace {
class PendingSolver final : public Solver {
 public:
  explicit PendingSolver(SolverOptions o) : opt_(std::move(o)) {}
  Result solve(const Stack&, const Wavelength&) override {
    throw NotImplemented(opt_.kind == SolverKind::RK
                             ? "RK solver (milestones M3/M4, Plan section 'Loeser A')"
                             : "RCWA solver (milestone M2, Plan section 'Loeser B')");
  }

 private:
  SolverOptions opt_;
};
}  // namespace

std::unique_ptr<Solver> makeSolver(const SolverOptions& opt) {
  return std::make_unique<PendingSolver>(opt);
}

}  // namespace fg
