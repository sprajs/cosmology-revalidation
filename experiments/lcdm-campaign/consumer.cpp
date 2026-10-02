// Bounded experiment transport. Production equations belong to Irreducible.
#include "irred/bao_thermal.hpp"
#include <cmath>
#include <cfenv>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
void require(bool value) { if (!value) throw std::runtime_error("typed input/native admission failed"); }
double number() { double x; require(bool(std::cin >> x) && std::isfinite(x)); return x; }
unsigned count(unsigned maximum) { unsigned n; require(bool(std::cin >> n) && n <= maximum); return n; }
void vector(const std::vector<double>& v) {
  std::cout << '[';
  for (size_t i=0;i<v.size();++i) std::cout << (i?",":"") << v[i];
  std::cout << ']';
}
void state(const bao::OutputState& s) {
  std::cout<<"{\"availability\":"<<unsigned(s.availability)
    <<",\"status\":"<<unsigned(s.status)<<",\"numerical_status\":"<<unsigned(s.numerical_status)<<'}';
}
void model(const cosmology::ThermalObservableRequest& r, const char* id) {
  const auto& m=r.model;
  std::cout<<"{\"id\":\""<<id<<"\",\"H0\":"<<m.h0_km_s_mpc
    <<",\"omega_b\":"<<m.physical_baryon_density<<",\"omega_c\":"<<m.physical_cdm_density
    <<",\"Tcmb\":"<<m.tcmb_kelvin<<",\"omega_other\":"<<m.physical_massless_nonphoton_density
    <<",\"z_drag\":"<<r.z_drag<<",\"species\":[";
  for(size_t j=0;j<m.species.size();++j) {
    const auto& sp=m.species[j];
    std::cout<<(j?",":"")<<"{\"mass_eV\":"<<sp.mass_ev<<",\"temperature_K\":"
      <<sp.temperature_today_kelvin<<",\"g\":"<<sp.statistical_weight<<'}';
  }
  std::cout<<"]}";
}
}
int main() { try {
  std::cout << std::setprecision(17);
  const auto nm=count(2); require(nm==2);
  std::vector<cosmology::ThermalObservableRequest> models;
  for(unsigned i=0;i<nm;++i) {
    cosmology::ThermalObservableRequest r;
    auto& m=r.model;
    m.h0_km_s_mpc=number(); m.physical_baryon_density=number();
    m.physical_cdm_density=number(); m.tcmb_kelvin=number();
    m.physical_massless_nonphoton_density=number(); r.z_drag=number();
    const auto ns=count(1);
    for(unsigned j=0;j<ns;++j) m.species.push_back({number(),number(),number()});
    r.drag_origin="chosen supplied endpoint; physical drag prediction blocked";
    r.source_origin=i?"campaign massive FD control/v1":"campaign massless control/v1";
    models.push_back(r);
  }
  const auto n=count(13); require(n==13);
  bao::DensityInput input;
  for(unsigned i=0;i<n;++i) {
    auto z=number(); auto kind=count(2); auto y=number();
    input.queries.push_back({z,static_cast<bao::Observable>(kind)});
    input.observed.push_back(y); input.ordered_ids.push_back("DESI-DR2-row-"+std::to_string(i));
  }
  for(unsigned i=0;i<n*n;++i) input.covariance.push_back(number());
  std::string extra; require(!(std::cin>>extra));
  input.role=bao::RowRole::released_fitted_distance_summary;
  input.covariance_unit=bao::CovarianceUnit::dimensionless_ratio_squared;
  input.table_identity="bao_data@bb0c1c9009dc76d1391300e169e8df38fd1096db; hash-admitted original13";
  input.covariance_identity="hash-admitted full13 covariance in original order";
  input.ordering_provenance="unchanged mean row order; full covariance uses identical axes";
  input.calibration_provenance="released fitted compression; source template/reconstruction applicability conditional";
  input.dependence_provenance="full supplied covariance; other-probe covariance unknown; no joint inference";
  const bao::PreparationPolicy prep{13,169,8192,1048576,1e-8,numerics::Arithmetic::longdouble_cpu_v1};
  auto prepared=bao::prepare_density(std::move(input),prep);
  if(prepared.status()!=statistics::DensityStatus::finite) {
    std::cout<<"{\"schema_version\":1,\"accepted\":false,\"failure_stage\":\"observation_preparation\",\"preparation_status\":"
      <<unsigned(prepared.status())<<",\"numerical_status\":"<<unsigned(prepared.numerical_status())
      <<",\"slots\":[],\"callbacks\":0}\n";
    return 3;
  }
  bao::ThermalDensityPolicy p;
  p.maximum_models=2; p.maximum_queries=13; p.maximum_string_bytes=8192;
  p.maximum_native_bytes=16*1024*1024; p.maximum_total_callbacks=1000000000;
  p.maximum_forward_sensitivity=1e-8; p.maximum_projection_log_density_error=1e-8;
  p.requested=7;
  auto& q=p.predictions;
  q.absolute_tolerance_mpc=1e-10; q.relative_tolerance=2e-12;
  q.absolute_tolerance_ratio=1e-12; q.relative_tolerance_ratio=5e-12;
  q.maximum_callbacks_per_point=100000000; q.maximum_total_callbacks=500000000;
  q.thermal.absolute_tolerance=1e-14; q.thermal.relative_tolerance=2e-14;
  q.thermal.maximum_total_callbacks=500000000;
  auto b=prepared.evaluate_thermal(models,p);
  const char* ids[]={"thermal-massless-supplied-drag","thermal-massive-FD-supplied-drag"};
  std::cout<<"{\"schema_version\":1,\"method\":\""<<bao::thermal_density_id<<"\",\"mapping\":\""<<cosmology::thermal_physical_mapping_id
    <<"\",\"requested\":"<<p.requested<<",\"batch_status\":"<<unsigned(b.status)
    <<",\"numerical_status\":"<<unsigned(b.numerical_status)<<",\"arithmetic\":{\"id\":\""<<prepared.metadata().arithmetic_id
    <<"\",\"double_mantissa_bits\":"<<std::numeric_limits<double>::digits
    <<",\"long_double_mantissa_bits\":"<<std::numeric_limits<long double>::digits
    <<",\"long_double_max_exponent\":"<<std::numeric_limits<long double>::max_exponent
    <<",\"round_to_nearest\":"<<(std::fegetround()==FE_TONEAREST?"true":"false")
    <<"},\"producer_policy\":{\"absolute_mpc\":"<<q.absolute_tolerance_mpc<<",\"relative\":"<<q.relative_tolerance
    <<",\"ratio_absolute\":"<<q.absolute_tolerance_ratio<<",\"ratio_relative\":"<<q.relative_tolerance_ratio
    <<",\"moment_absolute\":"<<q.thermal.absolute_tolerance<<",\"moment_relative\":"<<q.thermal.relative_tolerance
    <<",\"maximum_callbacks_per_point\":"<<q.maximum_callbacks_per_point
    <<",\"maximum_callbacks_per_model\":"<<q.maximum_total_callbacks<<",\"maximum_callbacks_total\":"<<p.maximum_total_callbacks
    <<"},\"policy_metadata\":{\"momentum_method\":\""<<cosmology::thermal_momentum_method_id(q.thermal.momentum_method)
    <<"\",\"provider_maximum_depth\":"<<q.maximum_depth<<",\"provider_maximum_points\":"<<q.maximum_points
    <<",\"provider_maximum_native_bytes\":"<<q.maximum_native_bytes
    <<",\"momentum_maximum_callbacks_per_evaluation\":"<<q.thermal.maximum_callbacks_per_evaluation
    <<",\"momentum_maximum_total_callbacks\":"<<q.thermal.maximum_total_callbacks
    <<",\"momentum_maximum_depth\":"<<q.thermal.maximum_depth<<",\"momentum_maximum_points\":"<<q.thermal.maximum_points
    <<",\"momentum_maximum_species\":"<<q.thermal.maximum_species<<",\"momentum_maximum_native_bytes\":"<<q.thermal.maximum_native_bytes
    <<",\"density_maximum_models\":"<<p.maximum_models<<",\"density_maximum_queries\":"<<p.maximum_queries
    <<",\"density_maximum_string_bytes\":"<<p.maximum_string_bytes<<",\"density_maximum_native_bytes\":"<<p.maximum_native_bytes
    <<",\"density_maximum_forward_sensitivity\":"<<p.maximum_forward_sensitivity
    <<",\"density_maximum_projection_log_density_error\":"<<p.maximum_projection_log_density_error
    <<",\"preparation_maximum_queries\":"<<prep.maximum_queries<<",\"preparation_maximum_matrix_elements\":"<<prep.maximum_matrix_elements
    <<",\"preparation_maximum_string_bytes\":"<<prep.maximum_string_bytes<<",\"preparation_maximum_native_bytes\":"<<prep.maximum_native_bytes
    <<",\"preparation_maximum_forward_sensitivity\":"<<prep.maximum_forward_sensitivity
    <<"},\"model_order\":[\""<<ids[0]<<"\",\""<<ids[1]<<"\"],\"model_sources\":[";
  for(size_t i=0;i<models.size();++i){if(i)std::cout<<',';model(models[i],ids[i]);}
  std::cout<<"],\"queries\":[";
  const auto& source=prepared.source();
  for(size_t i=0;i<source.queries.size();++i) {
    const auto& r=source.queries[i];
    std::cout<<(i?",":"")<<"{\"z\":"<<r.z<<",\"kind\":"<<unsigned(r.observable)
      <<",\"observed\":"<<source.observed[i]<<",\"id\":\""<<source.ordered_ids[i]<<"\"}";
  }
  std::cout<<"],\"slots\":[";
  bool accepted=b.status==statistics::DensityStatus::finite && b.numerical_status==numerics::Status::ok && b.slots.size()==2;
  for(size_t i=0;i<b.slots.size();++i) {
    const auto& s=b.slots[i]; if(i)std::cout<<',';
    std::cout<<"{\"source_index\":"<<i<<",\"model_source\":";model(s.source,ids[i]);
    std::cout<<",\"numerical_status\":"<<unsigned(s.numerical_status)<<",\"preparation_status\":"<<unsigned(s.preparation_status)
      <<",\"predictions_state\":";state(s.predictions_state);
    std::cout<<",\"residuals_state\":";state(s.residuals_state);
    std::cout<<",\"density_state\":";state(s.density_state);
    std::cout<<",\"density_status\":";if(s.result)std::cout<<unsigned(s.result->density.status);else std::cout<<"null";
    std::cout<<",\"predictions\":";
    vector(s.predictions); std::cout<<",\"residuals\":";vector(s.residuals);
    std::cout<<",\"callbacks\":"<<s.callbacks<<",\"preparation_callbacks\":"<<s.preparation_callbacks<<",\"outer_callbacks\":"<<s.outer_callbacks<<",\"momentum_callbacks\":"<<s.momentum_callbacks;
    std::cout<<",\"projection_estimate\":"; if(s.projection_log_density_error_estimate)std::cout<<*s.projection_log_density_error_estimate;else std::cout<<"null";
    const bool ok=s.numerical_status==numerics::Status::ok && s.result && s.result->density.status==statistics::DensityStatus::finite && s.predictions.size()==13;
    accepted=accepted&&ok;
    std::cout<<",\"density\":";
    if(ok) { const auto& d=*s.result; std::cout<<"{\"quadratic\":"<<d.quadratic<<",\"log_determinant\":"<<d.log_determinant<<",\"normalization\":"<<d.normalization<<",\"log_density\":"<<d.density.log_value<<'}'; }
    else std::cout<<"null";
    std::cout<<'}';
  }
  std::cout<<"],\"callbacks\":"<<b.callbacks<<",\"preparation_callbacks\":"<<b.preparation_callbacks
    <<",\"outer_callbacks\":"<<b.outer_callbacks<<",\"momentum_callbacks\":"<<b.momentum_callbacks
    <<",\"accepted\":"<<(accepted?"true":"false")<<"}\n";
  // Output keeps partial failures even when one model fails numerical admission.
  return accepted?0:3;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 2; } }
