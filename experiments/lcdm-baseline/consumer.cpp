// Experiment-specific typed transport; physical equations are installed Irred.
#include "irred/bao_conditional.hpp"
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
using namespace irred;
namespace {
void require(bool b) {
  if (!b)
    throw std::runtime_error("invalid typed input or native status");
}
double number() {
  double x;
  require(bool(std::cin >> x) && std::isfinite(x));
  return x;
}
std::string token() {
  std::string s;
  require(bool(std::cin >> std::quoted(s)) && s.size() <= 4096 && !s.empty());
  return s;
}
void values(const std::vector<double> &v) {
  std::cout << '[';
  for (size_t j = 0; j < v.size(); ++j) {
    if (j)
      std::cout << ',';
    std::cout << v[j];
  }
  std::cout << ']';
}
} // namespace
int main(int argc, char **argv) {
  try {
    std::cout << std::setprecision(17);
    cosmology::SoundHorizonRequest point{
        {number(), number(), number(), number(), number()}, number(), token()};
    unsigned n = 0, nz = 0;
    require(bool(std::cin >> n >> nz) && n == 13 && nz > 0 && nz <= 16);
    std::vector<double> z;
    for (unsigned j = 0; j < nz; ++j)
      z.push_back(number());
    bao::DensityInput d;
    for (unsigned j = 0; j < n; ++j) {
      d.ordered_ids.push_back(token());
      double rz = number();
      unsigned kind = 0;
      require(bool(std::cin >> kind) && kind <= 2);
      d.queries.push_back({rz, static_cast<bao::Observable>(kind)});
      d.observed.push_back(number());
    }
    for (unsigned j = 0; j < n * n; ++j)
      d.covariance.push_back(number());
    std::string trailing;
    require(!(std::cin >> trailing));
    d.role = bao::RowRole::released_fitted_distance_summary;
    d.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
    d.table_identity =
        "CobayaSampler/bao_data@bb0c1c9009dc76d1391300e169e8df38fd1096db/"
        "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt "
        "SHA2569ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a64058"
        "5";
    d.covariance_identity = "same release full13 covariance "
                            "SHA256252a143274c8a07c78694c119617d36594f6d7965d00"
                            "319ca611c6ffb886e509";
    d.ordering_provenance = "exact13 mean row order, bound to full covariance "
                            "axis order by reviewed release";
    d.calibration_provenance =
        "released fitted ratio compression; release "
        "reconstruction/fiducial/ruler conventions remain conditional";
    d.dependence_provenance = "full declared13 covariance; cross-probe overlap "
                              "unknown, no independent joint probe claim";
    auto prepared = bao::prepare_density(
        std::move(d), {13, 169, 8192, 1048576, 1e-8,
                       numerics::Arithmetic::longdouble_cpu_v1});
    require(prepared.status() == statistics::DensityStatus::finite);
    cosmology::EarlyLatePolicy early{
        1e-12,   2e-14,
        1e-14,   2e-13,
        1000000, 4000000,
        48,      64,
        4194304, {1e-12, 2e-14, 1000000, 48, 64, 4000000, 4194304}};
    require(argc == 1 || (argc == 2 && std::string(argv[1]) ==
                                           "--underbudget-producer-control"));
    if (argc == 2)
      early.relative_tolerance_ratio = 2e-14;
    bao::ConditionalDensityPolicy policy;
    policy.predictions = early;
    policy.maximum_models = 2;
    policy.maximum_queries = 13;
    policy.maximum_string_bytes = 8192;
    policy.maximum_native_bytes = 8388608;
    policy.maximum_total_callbacks = 4000000;
    policy.maximum_forward_sensitivity = 1e-8;
    policy.requested = 7;
    auto twice = point;
    twice.model.h0_km_s_mpc *= 2;
    cosmology::SoundHorizonRequest points[] = {point, twice};
    auto batch = prepared.evaluate_conditional(points, policy);
    require(batch.status == statistics::DensityStatus::finite &&
            batch.numerical_status == numerics::Status::ok &&
            batch.slots.size() == 2);
    std::cout << "{\"schema_version\":1,\"scientific_ids\":{\"early_late_"
                 "equation\":\""
              << cosmology::early_late_equation_id
              << "\",\"conditional_density\":\"" << bao::conditional_density_id
              << "\",\"physical_model\":\"" << cosmology::sound_horizon_model_id
              << "\"},\"arithmetic\":{\"density_arithmetic_id\":\""
              << prepared.metadata().arithmetic_id
              << "\",\"double_mantissa_bits\":"
              << std::numeric_limits<double>::digits
              << ",\"long_double_mantissa_bits\":"
              << std::numeric_limits<long double>::digits
              << ",\"long_double_max_exponent\":"
              << std::numeric_limits<long double>::max_exponent
              << ",\"round_to_nearest\":"
              << (std::fegetround() == FE_TONEAREST ? "true" : "false")
              << "},\"producer_policy\":{\"distance_"

                 "absolute_mpc\":"
              << early.absolute_tolerance_mpc
              << ",\"distance_relative\":" << early.relative_tolerance
              << ",\"ratio_absolute\":" << early.absolute_tolerance_ratio
              << ",\"ratio_relative\":" << early.relative_tolerance_ratio
              << ",\"sound_absolute_mpc\":"
              << early.sound.absolute_tolerance_mpc
              << ",\"sound_relative\":" << early.sound.relative_tolerance
              << "},\"density\":[";
    for (unsigned i = 0; i < 2; ++i) {
      if (i)
        std::cout << ',';
      const auto &s = batch.slots[i];
      std::cerr << "density point " << i << " numerical_status "
                << unsigned(s.numerical_status) << " prediction "
                << unsigned(s.predictions_state.numerical_status) << " density "
                << unsigned(s.density_state.numerical_status) << " callbacks "
                << s.callbacks << "\n";
      if (s.numerical_status != numerics::Status::ok) {
        std::vector<double> zs;
        unsigned diagnostic_mask = 0;
        for (const auto &q : prepared.source().queries) {
          zs.push_back(q.z);
          diagnostic_mask |= cosmology::early_late_mask(
              q.observable == bao::Observable::transverse_over_ruler
                  ? cosmology::EarlyLateOutput::dm_over_rs
              : q.observable == bao::Observable::hubble_over_ruler
                  ? cosmology::EarlyLateOutput::dh_over_rs
                  : cosmology::EarlyLateOutput::dv_over_rs);
        }
        const auto diagnostic = cosmology::evaluate_early_late(
            points[i], zs, diagnostic_mask, early);
        std::cerr << "diagnostic provider batch " << unsigned(diagnostic.status)
                  << " ruler "
                  << (diagnostic.ruler ? int(diagnostic.ruler->status) : -1)
                  << " callbacks " << diagnostic.callbacks << "\n";
        for (size_t j = 0; j < diagnostic.rows.size(); ++j)
          for (unsigned k = 6; k < 9; ++k)
            if (diagnostic_mask & (1u << k)) {
              const auto &v = diagnostic.rows[j].outputs[k];
              std::cerr << "row " << j << " output " << k << " cause "
                        << unsigned(v.status) << " empirical_error "
                        << v.error_estimate << "\n";
            }
      }
      require(s.projection_log_density_error_estimate.has_value() && s.result &&
              s.result->density.status == statistics::DensityStatus::finite &&
              s.predictions.size() == 13 && s.residuals.size() == 13 &&
              s.predictions_state.availability ==
                  cosmology::Availability::available &&
              s.density_state.availability ==
                  cosmology::Availability::available);
      const auto &r = *s.result;
      std::cout << "{\"status\":\"ok\",\"predictions\":";
      values(s.predictions);
      std::cout << ",\"quadratic\":" << r.quadratic
                << ",\"log_determinant\":" << r.log_determinant
                << ",\"normalization\":" << r.normalization
                << ",\"log_density\":" << r.density.log_value
                << ",\"projection_estimate\":"
                << *s.projection_log_density_error_estimate
                << ",\"callbacks\":" << s.callbacks << '}';
    }
    std::cout << "],\"background\":[";
    size_t callbacks = batch.callbacks;
    unsigned mask =
        cosmology::early_late_mask(cosmology::EarlyLateOutput::e) |
        cosmology::early_late_mask(cosmology::EarlyLateOutput::dm_mpc) |
        cosmology::early_late_mask(cosmology::EarlyLateOutput::dl_mpc);
    for (unsigned i = 0; i < 2; ++i) {
      if (i)
        std::cout << ',';
      early.maximum_total_callbacks =
          std::min(size_t(4000000) - callbacks, early.maximum_total_callbacks);
      auto b = cosmology::evaluate_early_late(points[i], z, mask, early);
      require(b.status == numerics::Status::ok && b.rows.size() == nz);
      callbacks += b.callbacks;
      std::cout << '[';
      for (unsigned j = 0; j < nz; ++j) {
        if (j)
          std::cout << ',';
        std::cout << "{\"status\":\"ok\",\"z\":" << z[j];
        for (auto o :
             {cosmology::EarlyLateOutput::e, cosmology::EarlyLateOutput::dm_mpc,
              cosmology::EarlyLateOutput::dl_mpc}) {
          const auto &v = b.rows[j].outputs[unsigned(o)];
          require(v.status == numerics::Status::ok && v.value);
          std::cout << ",\""
                    << (o == cosmology::EarlyLateOutput::e        ? "E"
                        : o == cosmology::EarlyLateOutput::dm_mpc ? "DM"
                                                                  : "DL")
                    << "\":" << *v.value;
        }
        std::cout << '}';
      }
      std::cout << ']';
    }
    std::cout << "],\"callbacks\":" << callbacks << "}\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 2;
  }
}
