// Source-only experiment transport draft; production mathematics is Irreducible's.
// Frozen distinct6f GaussianBox SDK; source-only until an explicit compute lease.
#include "irred/gaussian_box.hpp"
#include "bounds.hpp"
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <initializer_list>
#include <limits>
#include <locale>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/resource.h>
#include <utility>
#include <vector>

#if !defined(RELEASED_BOX_ENGINE_REVISION) || !defined(RELEASED_BOX_BUILD_ID) || \
    !defined(RELEASED_BOX_MANIFEST_SHA256) || !defined(RELEASED_BOX_ARCHIVE_SHA256) || \
    !defined(RELEASED_BOX_CLI_SHA256) || !defined(RELEASED_BOX_HEADER_SHA256)
#error Supply the frozen6f GaussianBox SDK identity declarations; historical SDKs are not admitted.
#endif

namespace {
using namespace irred;
using namespace irred::statistics;
constexpr std::size_t n=3492, original_p=47, p=46;
constexpr std::size_t native_payload=std::size_t(2)*1024*1024*1024;
constexpr std::size_t address_limit=std::size_t(3)*1024*1024*1024;
constexpr std::size_t output_limit=1024*1024;
constexpr unsigned native_seconds=900;
constexpr std::string_view interface_id="released-fixed44-box-native/v3";
constexpr std::string_view contract_sha="ee885213f48cdb04d3b6d0bada41fea3d65524feb9f33b0d3bda985645aa0741";
constexpr std::string_view request_sha="664b7054f1d2c3d23ae7af170c163e8d717213e2ccf05e5945e1971b56635de8";
// These declarations are frozen source identities. The controller verifies the
// actual source/archive/header/toolchain and consumer executable before/after.
static_assert(std::string_view(RELEASED_BOX_ENGINE_REVISION)=="6f869532c1951ed1afd9f2506b5d05c6cfd03c82");
static_assert(std::string_view(RELEASED_BOX_BUILD_ID)=="f22c25423cfb9cbac3c2b91a4e514b13ce604e92e7010f55a9aa0bdd42f40f59");
static_assert(std::string_view(RELEASED_BOX_MANIFEST_SHA256)=="f2e4d6a22257c13f8ab46de51cdfc92654d9857a82a7bc3ea5b8c0edf234f062");
static_assert(std::string_view(RELEASED_BOX_ARCHIVE_SHA256)=="1cb2b85ad292334f3b5041d669187b04a6fca0916bca8c41070b0ff878499dff");
static_assert(std::string_view(RELEASED_BOX_CLI_SHA256)=="09b5bfd05ce6c057423b1d41f6f81f8e4d0db44dfd35424fdccd0b34fd9b22d9");
static_assert(std::string_view(RELEASED_BOX_HEADER_SHA256)=="8d54de01dffc17a11c1cf625afe19df3655b375b7bcdb9d977f8b97baa86b3ce");
constexpr std::string_view guide_sha="7db2a06a5ef95c07729d0dbec4acd7d06e3b0eeb0678137e3d5c06e1d1a35571";
constexpr std::string_view compiler_sha="f04191f6a7b2cd7d9a62e1745872b8a6088791e5af6955488c69c9b2c4668bc9";
constexpr std::string_view standard_library_sha="f5fc7380f2ae46fa4053a64be04e7b98109f1066a4bbfff3c37042488aa0be0e";
constexpr std::string_view measure="46-dimensional Lebesgue measure in original ordered coordinates0..43,45,46";
constexpr std::string_view fixed_provenance="original44=literal0 point mass outside active46 measure; physical identity unresolved";
constexpr std::string_view prior_identity="lstsq_results.txt@c447f0fea703fcd0fff57de5000947b5ca81286b;SHA256=37d2d423d06b6a2100c47578eb9b1c566a575caf6ade0ed61f6b9586d2a35c95;literal closed endpoints";
static_assert(released_box_transport::active_original.size()==p);
static_assert(released_box_transport::active_original[43]==43&&
    released_box_transport::active_original[44]==45&&released_box_transport::active_original[45]==46);

void require(bool yes,const char* message) { if(!yes) throw std::runtime_error(message); }
bool normal(double x) { return std::isfinite(x)&&(x==0||std::fpclassify(x)==FP_NORMAL); }
void text(std::string_view s) {
  std::cout<<'"';
  for(unsigned char c:s) {
    if(c=='"'||c=='\\') std::cout<<'\\'<<char(c);
    else if(c<32||c>=127) {
      const char* h="0123456789abcdef";
      std::cout<<"\\u00"<<h[c/16]<<h[c%16];
    } else std::cout<<char(c);
  }
  std::cout<<'"';
}
void scalar(double x) { if(normal(x)) std::cout<<x;else std::cout<<"null"; }
void interval(const BoxInterval& x) {
  std::cout<<"{\"lower\":";scalar(x.lower);std::cout<<",\"upper\":";scalar(x.upper);std::cout<<'}';
}
bool valid_interval(const BoxInterval& x) { return normal(x.lower)&&normal(x.upper)&&x.lower<=x.upper; }
void values(std::span<const double> v) {
  std::cout<<'[';for(std::size_t i=0;i<v.size();++i) {if(i)std::cout<<',';scalar(v[i]);}std::cout<<']';
}
void intervals(std::span<const BoxInterval> v) {
  std::cout<<'[';for(std::size_t i=0;i<v.size();++i) {if(i)std::cout<<',';interval(v[i]);}std::cout<<']';
}
template<class T> void status(std::optional<T> x) {if(x)std::cout<<unsigned(*x);else std::cout<<"null";}
void bound(std::optional<std::size_t> x) {if(x)std::cout<<*x;else std::cout<<"null";}
bool hex_identity(std::string_view s,std::size_t width) {
  if(s.size()!=width)return false;
  for(char c:s)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false;
  return true;
}
void limits() {
  for(const auto& item:std::initializer_list<std::pair<int,rlim_t>>{
      {RLIMIT_AS,address_limit},{RLIMIT_CPU,native_seconds},{RLIMIT_FSIZE,output_limit}}) {
    const rlimit limit{item.second,item.second};
    require(setrlimit(item.first,&limit)==0,"process resource limit refused");
  }
  for(const char* name:{"OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"})
    require(setenv(name,"1",1)==0,"single-thread environment refused");
}
std::vector<double> read(const std::filesystem::path& directory,const char* name,std::size_t count) {
  const auto path=directory/name;
  require(std::filesystem::is_regular_file(std::filesystem::symlink_status(path)),"canonical file must be regular and not symlink");
  const auto bytes=count*sizeof(double);
  require(std::filesystem::file_size(path)==bytes,"exact canonical byte size mismatch");
  std::vector<double> data(count);
  std::ifstream stream(path,std::ios::binary);
  stream.read(reinterpret_cast<char*>(data.data()),static_cast<std::streamsize>(bytes));
  require(bool(stream)&&stream.peek()==std::char_traits<char>::eof(),"bounded canonical read mismatch");
  for(double value:data)require(normal(value),"canonical finite-normal binary64 domain mismatch");
  return data;
}
BoxPolicy frozen_policy() {
  BoxPolicy policy;
  policy.design={n*n,native_payload,1e-10};
  policy.maximum_log_probability_width=1e-9;
  policy.maximum_quantile_width=1e-8;
  policy.cdf_absolute_radius=2e-11;
  policy.maximum_cdf_nodes=8192;
  policy.maximum_bisections=96;
  policy.maximum_cdf_evaluations=256;
  return policy;
}
struct Attempt {
  const char* stage="input_admission";
  std::string error;
  std::optional<DensityStatus> gaussian_status,design_status,box_status;
  std::optional<numerics::Status> gaussian_numerical,design_numerical;
  std::optional<DesignRank> rank;
  std::optional<double> triangular_condition,transpose_condition;
  std::optional<std::size_t> gaussian_bound,design_bound,box_bound,evaluation_bound;
  std::optional<GaussianBoxResult> result;
  bool input_complete=false,design_consumed=false,output_complete=false;
};
bool complete(const GaussianBoxResult& r,const BoxPolicy& policy) {
  if(r.status!=DensityStatus::finite||r.numerical_status!=numerics::Status::ok||
      r.stage!=BoxStage::complete||!r.gaussian_completion_available||!r.endpoint_margins_available||
      !r.rectangle_enclosure_available||!r.normalization_enclosures_available||
      !r.quantile_enclosure_available||!r.endpoint_cdf_enclosures_available||
      r.unboxed_mean.size()!=p||r.unboxed_variance.size()!=p||
      r.standardized_lower_margins.size()!=p||r.standardized_upper_margins.size()!=p)return false;
  for(double v:r.unboxed_mean)if(!normal(v))return false;
  for(double v:r.unboxed_variance)if(!normal(v)||v<=0)return false;
  for(const auto& v:r.standardized_lower_margins)if(!valid_interval(v))return false;
  for(const auto& v:r.standardized_upper_margins)if(!valid_interval(v))return false;
  for(double v:{r.minimum_quadratic,r.reported_profile_quadratic,r.log_design_precision_determinant,
      r.source_covariance_log_determinant,r.maximum_variance_sensitivity,r.profile_stationarity,
      r.excluded_mass_upper})if(!normal(v))return false;
  for(const auto& v:{r.box_probability,r.log_box_probability,r.log_prior_volume,r.log_relative_box_integral,
      r.log_prior_normalized_relative_evidence,r.log_observation_normalized_evidence,r.requested_quantile,
      r.lower_endpoint_cdf,r.upper_endpoint_cdf})if(!valid_interval(v))return false;
  return r.minimum_quadratic>=0&&r.reported_profile_quadratic>=0&&
      r.maximum_variance_sensitivity>=0&&r.maximum_variance_sensitivity<=policy.design.maximum_forward_sensitivity&&
      r.profile_stationarity>=0&&r.profile_stationarity<=policy.design.maximum_forward_sensitivity&&
      r.excluded_mass_upper>0&&r.excluded_mass_upper<1&&r.box_probability.lower>0&&r.box_probability.upper<=1&&
      r.log_box_probability.upper<=0&&
      r.log_box_probability.upper-r.log_box_probability.lower<=policy.maximum_log_probability_width&&
      r.requested_quantile.lower>=released_box_transport::lower.back()&&
      r.requested_quantile.upper<=released_box_transport::upper.back()&&
      r.requested_quantile.upper-r.requested_quantile.lower<=policy.maximum_quantile_width&&
      r.lower_endpoint_cdf.lower>=0&&r.lower_endpoint_cdf.upper<=1&&
      r.upper_endpoint_cdf.lower>=0&&r.upper_endpoint_cdf.upper<=1&&
      r.cdf_evaluations<=policy.maximum_cdf_evaluations&&
      r.cdf_node_evaluations<=policy.maximum_cdf_nodes*policy.maximum_cdf_evaluations&&
      r.bisections<=2*policy.maximum_bisections;
}
void emit(const Attempt& a,const BoxPolicy& policy) {
  const bool accepted=a.output_complete&&a.result&&a.gaussian_status==DensityStatus::finite&&
      a.design_status==DensityStatus::finite&&a.box_status==DensityStatus::finite;
  std::cout<<std::setprecision(17)<<"{\"schema_version\":1,\"interface_id\":";text(interface_id);
  std::cout<<",\"contract_sha256\":";text(contract_sha);
  std::cout<<",\"request_sha256\":";text(request_sha);
  std::cout<<",\"sdk_identity\":{\"engine_revision\":";text(RELEASED_BOX_ENGINE_REVISION);
  std::cout<<",\"build_id\":";text(RELEASED_BOX_BUILD_ID);
  std::cout<<",\"manifest_sha256\":";text(RELEASED_BOX_MANIFEST_SHA256);
  std::cout<<",\"archive_sha256\":";text(RELEASED_BOX_ARCHIVE_SHA256);
  std::cout<<",\"cli_sha256\":";text(RELEASED_BOX_CLI_SHA256);
  std::cout<<",\"gaussian_box_header_sha256\":";text(RELEASED_BOX_HEADER_SHA256);
  std::cout<<",\"gaussian_box_guide_sha256\":";text(guide_sha);
  std::cout<<",\"compiler_executable_sha256\":";text(compiler_sha);
  std::cout<<",\"standard_library_sha256\":";text(standard_library_sha);
  std::cout<<",\"source_inventory_count\":311";
  std::cout<<",\"verification_owner\":\"distinct bounded controller; echoed compile-time declarations are not self-certified\"}";
  std::cout<<",\"target\":{\"observations\":3492,\"original_columns\":47,\"active_columns\":46,\"active_original_indices\":[";
  for(std::size_t j=0;j<p;++j)std::cout<<(j?",":"")<<released_box_transport::active_original[j];
  std::cout<<"],\"fixed_coordinates\":[{\"original_index\":44,\"value\":0,\"measure\":\"point mass outside active46 Lebesgue\"}],\"parameter_measure\":";text(measure);
  std::cout<<",\"prior_identity\":";text(prior_identity);
  std::cout<<",\"requested_marginal\":{\"original_index\":46,\"active_index\":45,\"cumulative_probability\":0.5,\"unit\":\"5log10(H0/[1km/s/Mpc])\"},\"support_lower\":";values(released_box_transport::lower);
  std::cout<<",\"support_upper\":";values(released_box_transport::upper);
  std::cout<<",\"endpoints\":\"closed literal binary64\",\"row_order\":\"released-row-0 through released-row-3491; unchanged original FITS order\"}";
  std::cout<<",\"source_identities\":{\"source_revision\":\"c447f0fea703fcd0fff57de5000947b5ca81286b\",\"packet_source_revision\":\"e881efd16cd6782435cfd3af4fa8fc2ee4d7ba0d\",\"lineage_sha256\":\"be07b74a1550dcc3501f332246213cea0a827e7ff1a7bb8ab3ac3f404401d7b6\",\"constrained_sha256\":\"2b8fab46b097dec96162130f85b9e1e675cdeb52af0032a631c666b834e5d499\",\"canonical_C_sha256\":\"ef4c2703047e3f8f3b74a77889d9aa74d66df99a249021f2ea9b3757b0c46471\",\"canonical_X_sha256\":\"7ac0bcf92b3658ff4af1002cfc75f49200fc1c45bf88dc1f91d19bb39a3ad288\",\"canonical_y_sha256\":\"e24516870e16f1fc4362f6c8837965163695ce21f425a7ffee7cd549cd43b61b\",\"hash_admission_owner\":\"controller before/after source and canonical byte verification\"}";
  std::cout<<",\"arithmetic\":{\"covariance_id\":\"F02/longdouble-cpu/v1\",\"qr_id\":\"longdouble-cpu/v1\",\"double_mantissa_bits\":"<<std::numeric_limits<double>::digits
    <<",\"long_double_mantissa_bits\":"<<std::numeric_limits<long double>::digits
    <<",\"long_double_max_exponent\":"<<std::numeric_limits<long double>::max_exponent
    <<",\"round_to_nearest\":"<<(std::fegetround()==FE_TONEAREST?"true":"false")
    <<",\"binary64_ieee\":"<<(std::numeric_limits<double>::is_iec559?"true":"false")
    <<",\"little_endian\":"<<(std::endian::native==std::endian::little?"true":"false")<<"}";
  std::cout<<",\"policy\":{\"maximum_elements\":"<<policy.design.maximum_elements<<",\"maximum_payload_bytes\":"<<policy.design.maximum_payload_bytes
    <<",\"maximum_forward_sensitivity\":"<<policy.design.maximum_forward_sensitivity
    <<",\"maximum_log_probability_width\":"<<policy.maximum_log_probability_width<<",\"maximum_quantile_width\":"<<policy.maximum_quantile_width
    <<",\"cdf_absolute_radius\":"<<policy.cdf_absolute_radius<<",\"maximum_cdf_nodes\":"<<policy.maximum_cdf_nodes
    <<",\"maximum_bisections_per_inverse\":"<<policy.maximum_bisections<<",\"maximum_cdf_evaluations\":"<<policy.maximum_cdf_evaluations
    <<",\"maximum_total_cdf_nodes\":"<<policy.maximum_cdf_nodes*policy.maximum_cdf_evaluations<<"}";
  std::cout<<",\"resources\":{\"jobs\":1,\"threads\":1,\"address_limit_bytes\":"<<address_limit<<",\"output_limit_bytes\":"<<output_limit<<",\"native_cpu_seconds\":"<<native_seconds<<",\"wall_timeout_owner\":\"controller\"}";
  std::cout<<",\"stages\":{\"last_stage\":";text(a.stage);std::cout<<",\"error\":";if(a.error.empty())std::cout<<"null";else text(a.error);
  std::cout<<",\"input_complete\":"<<(a.input_complete?"true":"false")<<",\"gaussian_status\":";status(a.gaussian_status);
  std::cout<<",\"gaussian_numerical_status\":";status(a.gaussian_numerical);
  std::cout<<",\"design_status\":";status(a.design_status);std::cout<<",\"design_numerical_status\":";status(a.design_numerical);
  std::cout<<",\"design_rank\":";status(a.rank);std::cout<<",\"equilibrated_triangular_condition_inf\":";if(a.triangular_condition)scalar(*a.triangular_condition);else std::cout<<"null";
  std::cout<<",\"equilibrated_transpose_triangular_condition_inf\":";if(a.transpose_condition)scalar(*a.transpose_condition);else std::cout<<"null";
  std::cout<<",\"design_method\":\"retained-whitened-pivoted-householder-qr/v1\",\"box_preparation_status\":";status(a.box_status);
  std::cout<<",\"box_preparation_numerical_status\":null,\"box_preparation_numerical_status_scope\":\"not exposed by current API\",\"design_consumed_by_box\":"<<(a.design_consumed?"true":"false")<<"}";
  std::cout<<",\"payload_bounds\":{\"gaussian_preparation\":";bound(a.gaussian_bound);std::cout<<",\"design_preparation\":";bound(a.design_bound);
  std::cout<<",\"box_preparation\":";bound(a.box_bound);std::cout<<",\"box_evaluation\":";bound(a.evaluation_bound);std::cout<<"}";
  std::cout<<",\"result_status\":";if(a.result)std::cout<<unsigned(a.result->status);else std::cout<<"null";
  std::cout<<",\"result_numerical_status\":";if(a.result)std::cout<<unsigned(a.result->numerical_status);else std::cout<<"null";
  std::cout<<",\"method_id\":";if(a.result)text(a.result->method_id);else std::cout<<"null";
  std::cout<<",\"enclosure_scope\":";if(a.result)text(a.result->enclosure_scope);else std::cout<<"null";
  std::cout<<",\"result_stage\":";if(a.result)std::cout<<unsigned(a.result->stage);else std::cout<<"null";
  std::cout<<",\"availability\":";
  if(a.result) {
    const auto& r=*a.result;
    std::cout<<"{\"gaussian_completion\":"<<(r.gaussian_completion_available?"true":"false")
      <<",\"endpoint_margins\":"<<(r.endpoint_margins_available?"true":"false")
      <<",\"rectangle_enclosure\":"<<(r.rectangle_enclosure_available?"true":"false")
      <<",\"normalization_enclosures\":"<<(r.normalization_enclosures_available?"true":"false")
      <<",\"quantile_enclosure\":"<<(r.quantile_enclosure_available?"true":"false")
      <<",\"endpoint_cdf_enclosures\":"<<(r.endpoint_cdf_enclosures_available?"true":"false")<<'}';
  } else std::cout<<"null";
  // Operation-owned availability preserves earned partial output without zeros
  // from unadmitted intervals being represented as physical predictions.
  std::cout<<",\"completion\":";
  if(a.result&&a.result->gaussian_completion_available) {
    const auto& r=*a.result;
    std::cout<<"{\"mean\":";values(r.unboxed_mean);std::cout<<",\"variance\":";values(r.unboxed_variance);
    std::cout<<",\"named_variance_original46\":";if(r.unboxed_variance.size()==p)scalar(r.unboxed_variance.back());else std::cout<<"null";
    std::cout<<",\"minimum_quadratic\":";scalar(r.minimum_quadratic);std::cout<<",\"reported_postcast_profile_quadratic\":";scalar(r.reported_profile_quadratic);
    std::cout<<",\"log_design_precision_determinant\":";scalar(r.log_design_precision_determinant);std::cout<<",\"source_covariance_log_determinant\":";scalar(r.source_covariance_log_determinant);
    std::cout<<",\"maximum_variance_sensitivity\":";scalar(r.maximum_variance_sensitivity);std::cout<<",\"profile_stationarity\":";scalar(r.profile_stationarity);std::cout<<'}';
  } else std::cout<<"null";
  std::cout<<",\"box_diagnostics\":";
  if(a.result&&a.result->endpoint_margins_available) {
    const auto& r=*a.result;
    std::cout<<"{\"standardized_lower_margins\":";intervals(r.standardized_lower_margins);std::cout<<",\"standardized_upper_margins\":";intervals(r.standardized_upper_margins);
    std::cout<<",\"excluded_mass_upper\":";scalar(r.excluded_mass_upper);std::cout<<",\"box_probability\":";
    if(r.rectangle_enclosure_available)interval(r.box_probability);else std::cout<<"null";
    std::cout<<",\"log_box_probability\":";if(r.rectangle_enclosure_available)interval(r.log_box_probability);else std::cout<<"null";
    std::cout<<",\"log_prior_volume\":";interval(r.log_prior_volume);std::cout<<'}';
  } else std::cout<<"null";
  std::cout<<",\"normalizations\":";
  if(a.result&&a.result->normalization_enclosures_available) {
    const auto& r=*a.result;std::cout<<"{\"log_relative_box_integral\":";interval(r.log_relative_box_integral);
    std::cout<<",\"log_prior_normalized_relative_evidence\":";interval(r.log_prior_normalized_relative_evidence);
    std::cout<<",\"log_observation_normalized_evidence\":";interval(r.log_observation_normalized_evidence);std::cout<<'}';
  } else std::cout<<"null";
  std::cout<<",\"median\":";
  if(a.result&&a.result->quantile_enclosure_available) {
    const auto& r=*a.result;std::cout<<"{\"original_index\":46,\"active_index\":45,\"cumulative_probability\":0.5,\"quantile\":";interval(r.requested_quantile);
    std::cout<<",\"lower_endpoint_cdf\":";if(r.endpoint_cdf_enclosures_available)interval(r.lower_endpoint_cdf);else std::cout<<"null";
    std::cout<<",\"upper_endpoint_cdf\":";if(r.endpoint_cdf_enclosures_available)interval(r.upper_endpoint_cdf);else std::cout<<"null";std::cout<<'}';
  } else std::cout<<"null";
  std::cout<<",\"work\":{\"cdf_node_evaluations\":";if(a.result)std::cout<<a.result->cdf_node_evaluations;else std::cout<<"null";
  std::cout<<",\"cdf_evaluations\":";if(a.result)std::cout<<a.result->cdf_evaluations;else std::cout<<"null";
  std::cout<<",\"bisections\":";if(a.result)std::cout<<a.result->bisections;else std::cout<<"null";
  std::cout<<"},\"output_complete\":"<<(a.output_complete?"true":"false")<<",\"accepted\":"<<(accepted?"true":"false")
    <<",\"qualification\":\"native conditional completion/enclosure only; independent original-input comparison unassessed; observational qualification blocked\"}\n";
}
}
int main(int argc,char** argv) {
  std::cout.imbue(std::locale::classic());
  Attempt attempt;const auto policy=frozen_policy();
  try {
    limits();
    require(argc==2&&std::string_view(argv[1]).size()<=4096,"exact canonical input directory required");
    require(std::endian::native==std::endian::little&&sizeof(double)==8&&std::numeric_limits<double>::is_iec559&&
        std::numeric_limits<double>::digits==53&&std::numeric_limits<long double>::digits>=64&&
        std::numeric_limits<long double>::max_exponent>=16384&&std::fegetround()==FE_TONEAREST,"arithmetic profile refused");
    require(hex_identity(RELEASED_BOX_ENGINE_REVISION,40)&&hex_identity(RELEASED_BOX_BUILD_ID,64)&&
        hex_identity(RELEASED_BOX_MANIFEST_SHA256,64)&&hex_identity(RELEASED_BOX_ARCHIVE_SHA256,64)&&
        hex_identity(RELEASED_BOX_CLI_SHA256,64)&&hex_identity(RELEASED_BOX_HEADER_SHA256,64),"new SDK identity declarations malformed");
    const std::filesystem::path directory(argv[1]);
    auto covariance=read(directory,"C.f64",n*n),original_design=read(directory,"X.f64",n*original_p),y=read(directory,"y.f64",n);
    std::vector<double> active_design;active_design.reserve(n*p);
    for(std::size_t i=0;i<n;++i)for(auto j:released_box_transport::active_original)
      active_design.push_back(original_design[i*original_p+j]);
    std::vector<double>().swap(original_design);
    attempt.input_complete=true;
    Metadata metadata;metadata.measure="product of unchanged3492 released observation coordinates";
    metadata.table_identity="SH0ES released compact y@c447f0fea703fcd0fff57de5000947b5ca81286b;canonical y SHA pinned by controller";
    metadata.uncertainty_identity="unchanged full3492x3492 C;canonical C SHA pinned by controller";
    metadata.ordering_provenance="original FITS row order; no masks, removal or permutation";
    metadata.calibration_provenance="released compact constraints unchanged; primary reduction discrepancies retained";
    metadata.dependence_provenance="full supplied covariance; unknown external/event/cross-probe covariance remains unknown";
    metadata.source_semantics="released fitted/calibrated compact products; conditional Gaussian target; not raw independent observations";
    for(std::size_t i=0;i<n;++i)metadata.ordered_ids.push_back("released-row-"+std::to_string(i));
    DesignMetadata design_metadata;design_metadata.residual_unit="released magnitude-like coordinate";
    design_metadata.design_identity="unchanged X=L.T;select only original44=literal0 out of3492x47;active0..43,45,46";
    design_metadata.dependence_identity=metadata.dependence_provenance;
    for(auto j:released_box_transport::active_original) {
      design_metadata.ordered_parameter_ids.push_back("released-parameter-"+std::to_string(j));
      design_metadata.parameter_units.push_back(j==41||j==43?"mag/dex":j==46?"mag:5log10(H0/[1km/s/Mpc])":"mag");
    }
    attempt.stage="gaussian_preparation";
    attempt.gaussian_bound=gaussian_preparation_payload_bound(n,MatrixKind::covariance,numerics::Arithmetic::longdouble_cpu_v1,metadata);
    require(attempt.gaussian_bound&&*attempt.gaussian_bound<=native_payload,"Gaussian payload quota refused");
    auto gaussian=prepare_gaussian(covariance,MatrixKind::covariance,metadata,n*n,1e-10,numerics::Arithmetic::longdouble_cpu_v1);
    attempt.gaussian_status=gaussian.status();attempt.gaussian_numerical=gaussian.numerical_status();
    require(gaussian.status()==DensityStatus::finite,"Gaussian preparation refused");
    std::vector<double>().swap(covariance);
    attempt.stage="design_preparation";
    attempt.design_bound=DesignProfile::preparation_payload_bound(gaussian,p,design_metadata);
    require(attempt.design_bound&&*attempt.design_bound<=native_payload,"DesignProfile payload quota refused");
    auto design=DesignProfile::prepare(std::move(gaussian),active_design,metadata.ordered_ids,design_metadata,policy.design);
    attempt.design_status=design.status();attempt.design_numerical=design.numerical_status();attempt.rank=design.rank();
    if(design.status()==DensityStatus::finite) {
      attempt.triangular_condition=design.equilibrated_triangular_condition_inf();
      attempt.transpose_condition=design.equilibrated_transpose_triangular_condition_inf();
    }
    require(design.status()==DensityStatus::finite&&design.rank()==DesignRank::full_within_conditioning_contract,"DesignProfile preparation refused");
    std::vector<double>().swap(active_design);
    BoxSupport support;support.ordered_parameter_ids=design_metadata.ordered_parameter_ids;
    support.lower.assign(released_box_transport::lower.begin(),released_box_transport::lower.end());
    support.upper.assign(released_box_transport::upper.begin(),released_box_transport::upper.end());
    support.parameter_measure=std::string(measure);support.prior_identity=std::string(prior_identity);
    support.fixed_coordinate_provenance=std::string(fixed_provenance);
    attempt.stage="box_preparation";
    attempt.box_bound=GaussianBox::preparation_payload_bound(design,support);
    require(attempt.box_bound&&*attempt.box_bound<=native_payload,"GaussianBox preparation payload quota refused");
    auto box=GaussianBox::prepare(std::move(design),std::move(support),policy);
    attempt.box_status=box.status();attempt.design_consumed=box.status()==DensityStatus::finite;
    require(box.status()==DensityStatus::finite,"GaussianBox preparation refused");
    attempt.stage="box_evaluation";
    attempt.evaluation_bound=box.evaluation_payload_bound();
    require(attempt.evaluation_bound&&*attempt.evaluation_bound<=native_payload,"GaussianBox evaluation payload quota refused");
    attempt.result=box.evaluate(y,metadata.ordered_ids,BoxMarginalRequest{45,.5},policy);
    attempt.output_complete=complete(*attempt.result,policy);
    attempt.stage=attempt.output_complete?"native_complete":"native_refusal_or_output_postcondition";
  } catch(const std::exception& error) {attempt.error=std::string(error.what()).substr(0,256);}
  emit(attempt,policy);
  return attempt.output_complete?0:3;
}
