// Thin external CLASS mean -> existing full-C Irred Gaussian owner. SOURCE ONLY.
#include "irred/statistics.hpp"
#include <array>
#include <cfenv>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifndef IRRED_CLASS_BAO_BUILD_ID
#error "Bind the reviewed actual SDK build ID; no default/old artifact identity"
#endif
namespace {
using DS = irred::statistics::DensityStatus;
using NS = irred::numerics::Status;
constexpr std::size_t n = 13, maximum_models = 4, maximum_input_bytes = 16384;
constexpr std::size_t maximum_payload_bytes = 1048576;
constexpr double sensitivity = 1e-8;
constexpr const char *mean_sha = "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585";
constexpr const char *covariance_sha = "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509";
const std::array<std::string, n> ids = {
    "desi-dr2-00-DV_over_rs", "desi-dr2-01-DM_over_rs", "desi-dr2-02-DH_over_rs",
    "desi-dr2-03-DM_over_rs", "desi-dr2-04-DH_over_rs", "desi-dr2-05-DM_over_rs",
    "desi-dr2-06-DH_over_rs", "desi-dr2-07-DM_over_rs", "desi-dr2-08-DH_over_rs",
    "desi-dr2-09-DM_over_rs", "desi-dr2-10-DH_over_rs", "desi-dr2-11-DH_over_rs",
    "desi-dr2-12-DM_over_rs"};
const std::array<std::string, maximum_models> case_ids = {"anchor", "precision", "ns-minus", "ns-plus"};
void require(bool condition) { if (!condition) throw std::runtime_error("closed transport/domain refusal"); }
bool sha(const std::string &s) {
  if (s.size() != 64) return false;
  for (char c : s) if (!(c >= '0' && c <= '9') && !(c >= 'a' && c <= 'f')) return false;
  return true;
}
std::string word(std::istream &in, std::size_t limit = 64) {
  std::string s;
  require(bool(in >> s) && !s.empty() && s.size() <= limit);
  for (unsigned char c : s) require(c >= 0x21 && c <= 0x7e);
  return s;
}
double number(std::istream &in) {
  const auto text = word(in, 48);
  double value = 0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value,
                                     std::chars_format::general);
  require(result.ec == std::errc{} && result.ptr == text.data() + text.size() && std::isfinite(value));
  return value;
}
unsigned integer(std::istream &in) {
  const auto text = word(in, 8);
  unsigned value = 0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  require(result.ec == std::errc{} && result.ptr == text.data() + text.size());
  return value;
}
const char *density_status(DS s) {
  switch (s) {
  case DS::finite: return "finite";
  case DS::outside_support: return "outside_support";
  case DS::invalid_input: return "invalid_input";
  case DS::unsupported_domain: return "unsupported_domain";
  case DS::numerical_failure: return "numerical_failure";
  case DS::incompatible_metadata: return "incompatible_metadata";
  }
  return "unknown";
}
const char *numerical_status(NS s) {
  switch (s) {
  case NS::ok: return "ok";
  case NS::invalid_input: return "invalid_input";
  case NS::nonfinite_input: return "nonfinite_input";
  case NS::overflow: return "overflow";
  case NS::work_limit: return "work_limit";
  case NS::outside_domain: return "outside_domain";
  case NS::singular: return "singular";
  case NS::not_positive_definite: return "not_positive_definite";
  case NS::conditioning_budget_exceeded: return "conditioning_budget_exceeded";
  }
  return "unknown";
}
template <class Range> void values(const Range &range) {
  std::cout << '[';
  bool comma = false;
  for (const auto value : range) { if (comma) std::cout << ','; comma = true; std::cout << value; }
  std::cout << ']';
}
void strings() {
  std::cout << '[';
  for (std::size_t i = 0; i < n; ++i) { if (i) std::cout << ','; std::cout << std::quoted(ids[i]); }
  std::cout << ']';
}
struct Row {
  std::string id, state_sha;
  std::array<double, n> predictions{}, residuals{};
  bool input_complete = false, attempted = false, accepted = false;
  std::optional<irred::statistics::GaussianResult> result;
};
bool usable(const irred::statistics::GaussianResult &r) {
  return r.density.status == DS::finite && r.density.numerical_status == NS::ok &&
         std::isfinite(r.density.log_value) && std::isfinite(r.quadratic) && r.quadratic >= 0 &&
         std::isfinite(r.log_determinant) && std::isfinite(r.normalization) &&
         std::isfinite(r.backward_residual) && r.backward_residual >= 0 &&
         std::isfinite(r.estimated_forward_sensitivity) && r.estimated_forward_sensitivity >= 0 &&
         r.estimated_forward_sensitivity <= sensitivity;
}
} // namespace
int main(int argc, char **) {
  std::cout << std::setprecision(17);
  const char *stage = "transport";
  bool accepted = false, failed = false, data_complete = false;
  std::array<double, n> observed{};
  std::vector<Row> rows;
  std::optional<irred::statistics::Gaussian> gaussian;
  std::size_t preparation_attempts = 0, evaluation_attempts = 0;
  try {
    require(argc == 1 && sha(IRRED_CLASS_BAO_BUILD_ID));
    require(std::numeric_limits<double>::digits == 53 && std::fegetround() == FE_TONEAREST);
    std::string raw;
    char c = 0;
    while (std::cin.get(c)) { require(raw.size() < maximum_input_bytes); raw.push_back(c); }
    require(std::cin.eof());
    std::istringstream input(raw);
    require(word(input) == "CLASS_BAO_FIXED_MEAN_V1" && integer(input) == n);
    const auto models = integer(input);
    require(models > 0 && models <= maximum_models);
    require(word(input) == mean_sha && word(input) == covariance_sha);
    for (const auto &id : ids) require(word(input) == id);
    for (auto &value : observed) { value = number(input); require(value > 0); }
    std::vector<double> covariance(n * n);
    for (auto &value : covariance) value = number(input);
    data_complete = true;
    rows.reserve(models);
    for (unsigned i = 0; i < models; ++i) {
      rows.emplace_back();
      auto &row = rows.back();
      row.id = word(input);
      require(row.id == case_ids[i]);
      row.state_sha = word(input);
      require(sha(row.state_sha));
      for (std::size_t j = 0; j < n; ++j) {
        row.predictions[j] = number(input);
        require(row.predictions[j] > 0);
        row.residuals[j] = observed[j] - row.predictions[j];
        require(std::isfinite(row.residuals[j]));
      }
      row.input_complete = true;
    }
    std::string trailing;
    require(!(input >> trailing) && input.eof());
    irred::statistics::Metadata metadata;
    metadata.ordered_ids.assign(ids.begin(), ids.end());
    metadata.measure = "product-of-13-dimensionless-ratio-coordinates";
    metadata.table_identity = mean_sha;
    metadata.uncertainty_identity = covariance_sha;
    metadata.ordering_provenance = "exact DESI DR2 mean order and full covariance axes; final DH then DM";
    metadata.calibration_provenance = "released compression/reconstruction/fiducial calibration conditional";
    metadata.dependence_provenance = "full declared covariance; cross-probe dependence unresolved";
    metadata.source_semantics = "released fitted distance summaries; external fixed CLASS predictions";
    metadata.input_matrix_convention = "covariance";
    const auto arithmetic = irred::numerics::Arithmetic::binary64_legacy_v1;
    const auto payload = irred::statistics::gaussian_preparation_payload_bound(
        n, irred::statistics::MatrixKind::covariance, arithmetic, metadata);
    require(payload && *payload <= maximum_payload_bytes);
    stage = "preparation";
    ++preparation_attempts;
    gaussian.emplace(irred::statistics::prepare_gaussian(
        covariance, irred::statistics::MatrixKind::covariance, std::move(metadata), n * n,
        sensitivity, arithmetic));
    require(gaussian->status() == DS::finite && gaussian->numerical_status() == NS::ok);
    stage = "evaluation";
    const auto scratch = gaussian->evaluation_payload_bound(n, false);
    require(scratch && *scratch <= maximum_payload_bytes);
    accepted = true;
    for (auto &row : rows) {
      row.attempted = true;
      ++evaluation_attempts;
      row.result = gaussian->evaluate(row.residuals, ids, sensitivity);
      row.accepted = usable(*row.result);
      accepted = accepted && row.accepted;
    }
    stage = accepted ? "complete" : "evaluation";
  } catch (const std::exception &) { failed = true; accepted = false; }
    catch (...) { failed = true; accepted = false; }
  std::cout << "{\"schema_version\":1,\"interface_id\":\"external-class-bao-gaussian-native/v1\","
               "\"status\":" << std::quoted(accepted ? "accepted" : "refused")
            << ",\"stage\":" << std::quoted(stage)
            << ",\"sdk_build_id\":" << std::quoted(IRRED_CLASS_BAO_BUILD_ID)
            << ",\"preparation_status\":";
  if (gaussian) std::cout << std::quoted(density_status(gaussian->status())); else std::cout << "null";
  std::cout << ",\"preparation_numerical_status\":";
  if (gaussian) std::cout << std::quoted(numerical_status(gaussian->numerical_status())); else std::cout << "null";
  std::cout << ",\"arithmetic_id\":";
  if (gaussian) std::cout << std::quoted(gaussian->metadata().arithmetic_id); else std::cout << "null";
  std::cout << ",\"producer_policy\":{\"maximum_forward_sensitivity\":" << sensitivity
            << ",\"maximum_matrix_elements\":169,\"maximum_models\":4,\"maximum_input_bytes\":16384,"
               "\"maximum_library_payload_bytes\":1048576},\"ordered_ids\":";
  strings();
  std::cout << ",\"observed\":";
  if (data_complete) values(observed); else std::cout << "null";
  std::cout << ",\"work\":{\"preparation_attempts\":" << preparation_attempts
            << ",\"evaluation_attempts\":" << evaluation_attempts << "},\"rows\":[";
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const auto &row = rows[i];
    if (i) std::cout << ',';
    std::cout << "{\"index\":" << i << ",\"id\":" << std::quoted(row.id)
              << ",\"common_state_sha256\":" << std::quoted(row.state_sha)
              << ",\"input_complete\":" << (row.input_complete ? "true" : "false")
              << ",\"attempted\":" << (row.attempted ? "true" : "false")
              << ",\"accepted\":" << (row.accepted ? "true" : "false")
              << ",\"density_status\":";
    if (row.result) std::cout << std::quoted(density_status(row.result->density.status)); else std::cout << "null";
    std::cout << ",\"numerical_status\":";
    if (row.result) std::cout << std::quoted(numerical_status(row.result->density.numerical_status)); else std::cout << "null";
    std::cout << ",\"payload\":";
    if (row.accepted) {
      const auto &r = *row.result;
      std::cout << "{\"predictions\":"; values(row.predictions);
      std::cout << ",\"residuals\":"; values(row.residuals);
      std::cout << ",\"quadratic\":" << r.quadratic << ",\"log_determinant\":" << r.log_determinant
                << ",\"normalization\":" << r.normalization << ",\"log_density\":" << r.density.log_value
                << ",\"backward_residual\":" << r.backward_residual
                << ",\"estimated_forward_sensitivity\":" << r.estimated_forward_sensitivity << '}';
    } else std::cout << "null";
    std::cout << '}';
  }
  std::cout << "],\"error\":";
  if (failed) std::cout << "{\"kind\":\"native_refusal\",\"stage\":" << std::quoted(stage) << '}';
  else std::cout << "null";
  std::cout << "}\n";
  return accepted ? 0 : 2;
}
