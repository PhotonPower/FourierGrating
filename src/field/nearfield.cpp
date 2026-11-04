#include "fg/field/nearfield.hpp"

#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>

#include "fg/core/error.hpp"

namespace fg {

namespace {

/// Minimal NPY v1.0 writer (little endian, C order).
template <typename T>
void writeNpy(const std::string& file, const std::vector<T>& data, const std::vector<std::size_t>& shape,
              const char* descr) {
  std::ostringstream h;
  h << "{'descr': '" << descr << "', 'fortran_order': False, 'shape': (";
  for (std::size_t k = 0; k < shape.size(); ++k) h << shape[k] << (shape.size() == 1 ? "," : (k + 1 < shape.size() ? ", " : ""));
  h << "), }";
  std::string header = h.str();
  const std::size_t total = 10 + header.size() + 1;
  header.append((64 - total % 64) % 64, ' ');
  header.push_back('\n');
  std::ofstream out(file, std::ios::binary);
  FG_REQUIRE(out.good(), "cannot write " + file);
  out.write("\x93NUMPY", 6);
  const char ver[2] = {1, 0};
  out.write(ver, 2);
  const std::uint16_t len = static_cast<std::uint16_t>(header.size());
  out.write(reinterpret_cast<const char*>(&len), 2);
  out.write(header.data(), static_cast<std::streamsize>(header.size()));
  out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(T)));
}

}  // namespace

void NearFieldWriter::write(const FieldCut& cut, const std::string& dir) {
  std::filesystem::create_directories(dir);
  const std::size_t n1 = cut.axis1_um.size(), n2 = cut.axis2_um.size();
  const std::vector<std::size_t> shape = n2 ? std::vector<std::size_t>{n2, n1} : std::vector<std::size_t>{n1};
  const char* names[6] = {"Ex", "Ey", "Ez", "Hx", "Hy", "Hz"};
  std::vector<std::string> q;
  for (int c = 0; c < 3; ++c) {
    writeNpy(dir + "/" + names[c] + ".npy", cut.E[c], shape, "<c16");
    writeNpy(dir + "/" + names[3 + c] + ".npy", cut.H[c], shape, "<c16");
  }
  for (const char* n : names) q.push_back(n);
  const std::size_t N = cut.E[0].size();
  std::vector<double> e2(N), h2(N), sy(N);
  for (std::size_t k = 0; k < N; ++k) {
    e2[k] = std::norm(cut.E[0][k]) + std::norm(cut.E[1][k]) + std::norm(cut.E[2][k]);
    h2[k] = std::norm(cut.H[0][k]) + std::norm(cut.H[1][k]) + std::norm(cut.H[2][k]);
    sy[k] = 0.5 * std::real(cut.E[2][k] * std::conj(cut.H[0][k]) - cut.E[0][k] * std::conj(cut.H[2][k]));
  }
  writeNpy(dir + "/absE2.npy", e2, shape, "<f8");
  writeNpy(dir + "/absH2.npy", h2, shape, "<f8");
  writeNpy(dir + "/Sy.npy", sy, shape, "<f8");
  if (!cut.material.empty()) writeNpy(dir + "/material.npy", cut.material, shape, "<f8");
  q.insert(q.end(), {"absE2", "absH2", "Sy"});
  nlohmann::json meta;
  meta["cut"] = cut.cut;
  meta["dimension"] = "1D";
  meta["lambda_um"] = cut.lambda_um;
  meta["axes"]["axis1"] = {{"name", cut.axis1Name}, {"values_um", cut.axis1_um}};
  if (n2) meta["axes"]["axis2"] = {{"name", cut.axis2Name}, {"values_um", cut.axis2_um}};
  meta["quantities"] = q;
  meta["normalization"] = "E: |E_inc| = 1; H: Z0*H; Sy = 1/2 Re(E x H*)_y";
  std::ofstream(dir + "/meta.json") << meta.dump(2);
}

}  // namespace fg
