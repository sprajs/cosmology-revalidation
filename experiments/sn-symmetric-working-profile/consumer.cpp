// NEW stored S_hat target. Existing Irred owns every factor/profile operation.
#include "irred/statistics.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#ifndef IRRED_SN_BUILD_ID
#error "Exact reviewed f3b19/6d495 SDK build ID required"
#endif
namespace {
using DS = irred::statistics::DensityStatus;
using NS = irred::numerics::Status;
constexpr size_t selected_n = 1657, max_n = 1701, max_input = 128 * 1024 * 1024;
constexpr size_t max_library_payload = 256 * 1024 * 1024;
constexpr double sensitivity = 1e-8;
constexpr const char *table_sha = "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8";
constexpr const char *cov_sha = "abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc";
const std::array<std::string,4> cases = {"anchor","precision","ns-minus","ns-plus"};
void need(bool p) { if (!p) throw std::runtime_error("SN closed transport/domain refusal"); }
bool sha(const std::string &s) {
  if (s.size()!=64) return false;
  for(char c:s) if (!(c>='0'&&c<='9')&&!(c>='a'&&c<='f')) return false;
  return true;
}
struct Input {
  size_t bytes = 0;
  int get() { const int c=std::cin.get(); if(c!=EOF) need(++bytes<=max_input); return c; }
  std::string word(size_t cap=80) {
    std::string out; int c;
    do { c=get(); } while(c==' '||c=='\n'||c=='\t'||c=='\r');
    need(c!=EOF);
    do { need(c>=33&&c<=126&&out.size()<cap); out.push_back(static_cast<char>(c)); c=get(); }
    while(c!=EOF&&c!=' '&&c!='\n'&&c!='\t'&&c!='\r');
    return out;
  }
  size_t integer() {
    auto s=word(8); size_t v=0;
    auto r=std::from_chars(s.data(),s.data()+s.size(),v);
    need(r.ec==std::errc{}&&r.ptr==s.data()+s.size()); return v;
  }
  double number() {
    auto s=word(32); double v=0;
    auto r=std::from_chars(s.data(),s.data()+s.size(),v,std::chars_format::general);
    need(r.ec==std::errc{}&&r.ptr==s.data()+s.size()&&std::isfinite(v)); return v;
  }
  void end() { int c; while((c=get())!=EOF) need(c==' '||c=='\n'||c=='\t'||c=='\r'); need(std::cin.eof()); }
};
const char* status(DS x) {
  switch(x) {case DS::finite:return "finite";case DS::outside_support:return "outside_support";
    case DS::invalid_input:return "invalid_input";case DS::unsupported_domain:return "unsupported_domain";
    case DS::numerical_failure:return "numerical_failure";case DS::incompatible_metadata:return "incompatible_metadata";}
  return "unknown";
}
const char* status(NS x) {
  switch(x) {case NS::ok:return "ok";case NS::invalid_input:return "invalid_input";
    case NS::nonfinite_input:return "nonfinite_input";case NS::overflow:return "overflow";
    case NS::work_limit:return "work_limit";case NS::outside_domain:return "outside_domain";
    case NS::singular:return "singular";case NS::not_positive_definite:return "not_positive_definite";
    case NS::conditioning_budget_exceeded:return "conditioning_budget_exceeded";}
  return "unknown";
}
template<class R> void array(const R& r) {std::cout<<'[';bool comma=false;for(const auto& x:r){if(comma)std::cout<<',';comma=true;std::cout<<x;}std::cout<<']';}
struct Row {
  std::string state_sha;
  std::vector<double> mu, residual;
  bool complete=false,attempted=false,accepted=false;
  std::optional<irred::statistics::ProfileResult> result;
};
bool finite_payload(const irred::statistics::ProfileResult& r,size_t n) {
  if(r.status!=DS::finite||r.numerical_status!=NS::ok||r.adjusted_residuals.size()!=n||
     !std::isfinite(r.coefficient)||!std::isfinite(r.quadratic)||r.quadratic<0) return false;
  for(double x:r.adjusted_residuals) if(!std::isfinite(x))return false;
  for(double x:{r.backward_residual,r.estimated_forward_sensitivity,r.coefficient_solve_backward_residual,
                r.coefficient_solve_forward_sensitivity,r.residual_l1,r.solution_norm_inf,
                r.adjusted_residual_l1,r.adjusted_solution_norm_inf}) if(!std::isfinite(x)||x<0)return false;
  return true;
}
bool valid(const irred::statistics::ProfileResult& r,size_t n) {
  return finite_payload(r,n)&&r.estimated_forward_sensitivity<=sensitivity&&r.coefficient_solve_forward_sensitivity<=sensitivity;
}
void payload(const irred::statistics::ProfileResult& r) {
  std::cout<<"{\"M_coefficient\":"<<r.coefficient<<",\"quadratic\":"<<r.quadratic
    <<",\"relative_profile_score\":"<<-0.5*r.quadratic<<",\"adjusted_residuals\":";
  array(r.adjusted_residuals);
  std::cout<<",\"backward_residual\":"<<r.backward_residual
    <<",\"estimated_forward_sensitivity\":"<<r.estimated_forward_sensitivity
    <<",\"coefficient_solve_backward_residual\":"<<r.coefficient_solve_backward_residual
    <<",\"coefficient_solve_forward_sensitivity\":"<<r.coefficient_solve_forward_sensitivity
    <<",\"residual_l1\":"<<r.residual_l1<<",\"solution_norm_inf\":"<<r.solution_norm_inf
    <<",\"adjusted_residual_l1\":"<<r.adjusted_residual_l1
    <<",\"adjusted_solution_norm_inf\":"<<r.adjusted_solution_norm_inf<<'}';
}
}
int main(int argc,char** argv) {
  std::cout<<std::setprecision(17);
  size_t n=0,prepare_attempts=0,profile_prepare_attempts=0,evaluation_attempts=0;
  Input input; std::vector<size_t> indices;std::vector<std::string> ids;
  std::vector<double> observed;std::array<Row,4> rows;
  std::string selected_sha;const char* stage="transport";bool accepted=false,error=false;
  std::string control="none";bool raw_only=false,profile_payload=false;
  std::optional<irred::statistics::Gaussian> gaussian;
  std::optional<DS> preparation_status;
  std::optional<NS> preparation_numerical_status;
  std::optional<irred::statistics::ProfileOperator> profile;
  try {
    need((argc==1||argc==3)&&std::string(IRRED_SN_BUILD_ID)=="6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e");
    need(std::numeric_limits<double>::digits==53&&std::fegetround()==FE_TONEAREST);
    std::vector<double> covariance;
    if(argc==3&&std::string(argv[1])=="--control") {
      control=argv[2];need(control=="analytic"||control=="conditioning"||control=="indefinite"||control=="nonsymmetric"||
        control=="nonfinite"||control=="wrong-order"||control=="missing-row");
      n=control=="missing-row"?selected_n-1:selected_n;need(n==selected_n);
      indices.resize(n);ids.reserve(n);observed.resize(n);covariance.resize(n*n);
      for(size_t i=0;i<n;++i){indices[i]=i;ids.push_back("full-n-analytic:row:"+std::to_string(i));observed[i]=(double(int(i%17))-8)/16;
        for(size_t j=0;j<n;++j)covariance[i*n+j]=(i==j?2.0:0.0)+((i%2==j%2)?0.0625:-0.0625);}
      if(control=="wrong-order"){std::swap(indices[0],indices[1]);need(indices[0]<indices[1]);}
      if(control=="conditioning"||control=="indefinite") {
        std::fill(covariance.begin(),covariance.end(),0.0);
        for(size_t i=0;i<n;++i)covariance[i*n+i]=1.0;
        if(control=="conditioning")covariance.back()=std::ldexp(1.0,-48);
        else covariance[1]=covariance[n]=2.0; // exact projection of B01=3,B10=1; positive diagonal, second pivot fails
      }
      if(control=="nonsymmetric")covariance[1]+=0.015625;
      if(control=="nonfinite")covariance[0]=std::numeric_limits<double>::quiet_NaN();
      for(size_t k=0;k<4;++k){auto& row=rows[k];row.state_sha=std::string(64,'0');row.mu.resize(n);row.residual.resize(n);
        const double shift=k==1?1.0:k==2?-0.5:0.0;
        for(size_t i=0;i<n;++i){row.mu[i]=-shift;row.residual[i]=observed[i]+shift;}row.complete=true;}
    } else {
      raw_only=argc==3&&std::string(argv[1])=="--raw-B-refusal";need(argc==1||raw_only);
      need(input.word()==(raw_only?"SN_ORIGINAL_B_REFUSAL_ONLY_V1":"SN_SHAT_COMMON_M_V1"));
      n=input.integer();need(n==selected_n&&input.integer()==4);
      need(input.word()==table_sha&&input.word()==cov_sha);selected_sha=input.word();need(sha(selected_sha));
    indices.resize(n);ids.reserve(n);
    for(size_t i=0;i<n;++i) {indices[i]=input.integer();need(indices[i]<1701&&(i==0||indices[i]>indices[i-1]));
      ids.push_back(std::string(table_sha)+":row:"+std::to_string(indices[i]));}
    covariance.resize(n*n);
    if(raw_only) {
      std::ifstream binary(argv[2],std::ios::binary);need(bool(binary));
      for(auto& x:covariance){std::array<unsigned char,8> bytes{};binary.read(reinterpret_cast<char*>(bytes.data()),8);need(bool(binary));
        std::uint64_t bits=0;for(unsigned char b:bytes)bits=(bits<<8)|b;x=std::bit_cast<double>(bits);need(std::isfinite(x));}
      need(binary.get()==EOF);input.end();
    } else {
      observed.resize(n);for(auto& x:observed)x=input.number();
      for(auto& x:covariance)x=input.number();
    for(size_t i=0;i<4;++i) {
      need(input.word()==cases[i]);auto& row=rows[i];row.state_sha=input.word();need(sha(row.state_sha));
      row.mu.resize(n);row.residual.resize(n);
      for(size_t j=0;j<n;++j){row.mu[j]=input.number();row.residual[j]=observed[j]-row.mu[j];need(std::isfinite(row.residual[j]));}
      row.complete=true;
    }
    input.end();
    }
    }
    irred::statistics::Metadata m; m.ordered_ids=ids;m.measure="product d(magnitude)";
    m.table_identity=table_sha;m.uncertainty_identity=std::string(cov_sha)+(raw_only?";original-B-refusal-only-sha=":";NEW-selected-S_hat-f64be-sha=")+selected_sha;
    m.ordering_provenance="original table occurrence order; zHD>.01 OR calibrator; exact selected indices";
    m.calibration_provenance="released m_b_corr and CEPH_DIST; full STAT+SYS includes Cepheid-host covariance";
    m.dependence_provenance="SN+SH0ES shared events/host/systematics; cross-probe covariance unknown";
    m.source_semantics="released fitted summary; external CLASS means and source calibrator host moduli";
    m.input_matrix_convention="covariance";m.matrix_validation_scope=irred::statistics::MatrixValidationScope::selected_covariance_only;
    m.treatment="NEW unconstrained common-M relative profile; no prior/density/evidence";
    if(control!="none") {
      m.table_identity="synthetic-full-n-SN-route-control/v1";m.uncertainty_identity="synthetic-"+control;
      m.ordering_provenance="synthetic occurrence indices0..1656; no observational rows";
      m.calibration_provenance="synthetic control; calibration not applicable";
      m.dependence_provenance="synthetic same-target engineering check";
      m.source_semantics="synthetic full-n rank-one/diagonal/refusal fixture; not released input";
    }
    const auto arithmetic=irred::numerics::Arithmetic::binary64_legacy_v1;
    auto bound=irred::statistics::gaussian_preparation_payload_bound(n,irred::statistics::MatrixKind::covariance,arithmetic,m);
    need(bound&&*bound<=max_library_payload);
    stage="covariance-preparation";++prepare_attempts;
    gaussian.emplace(irred::statistics::prepare_gaussian(covariance,irred::statistics::MatrixKind::covariance,std::move(m),n*n,sensitivity,arithmetic));
    preparation_status=gaussian->status();preparation_numerical_status=gaussian->numerical_status();
    if(raw_only) {stage="original-B-preparation-only-refusal";need(gaussian->status()!=DS::finite);throw std::runtime_error("expected original B refusal; no profile or score attempted");}
    need(gaussian->status()==DS::finite&&gaussian->numerical_status()==NS::ok);
    const auto retained=gaussian->retained_payload_bound(),scratch=gaussian->evaluation_payload_bound(n,true);
    need(retained&&scratch&&*retained<=max_library_payload&&*scratch<=max_library_payload-*retained);
    stage="profile-preparation";++profile_prepare_attempts;
    std::vector<double> ones(n,1.0);
    profile.emplace(std::move(*gaussian).prepare_offset_profile(ones,ids,sensitivity));
    need(profile->status()==DS::finite&&profile->numerical_status()==NS::ok);
    need(std::isfinite(profile->gram())&&profile->gram()>0&&std::isfinite(profile->cached_response_backward_residual())&&
      profile->cached_response_backward_residual()>=0&&std::isfinite(profile->cached_response_forward_sensitivity())&&
      profile->cached_response_forward_sensitivity()>=0&&profile->cached_response_solution().size()==n);
    for(double x:profile->cached_response_solution())need(std::isfinite(x));
    profile_payload=true;
    need(profile->cached_response_forward_sensitivity()<=sensitivity);
    covariance.clear();covariance.shrink_to_fit();
    stage="evaluation";accepted=true;
    for(auto& row:rows){row.attempted=true;++evaluation_attempts;row.result=profile->evaluate(row.residual,ids,sensitivity);
      row.accepted=valid(*row.result,n);accepted=accepted&&row.accepted;if(!row.accepted)break;}
    stage=accepted?"complete":"evaluation";
  }catch(...){error=true;accepted=false;}
  std::cout<<"{\"interface\":\"SN-S_hat-common-M-native/v1\",\"target_id\":\"released-sn-symmetric-working-profile/v1\",\"control\":"<<std::quoted(control)
    <<",\"original_B_refusal_only\":"<<(raw_only?"true":"false")<<",\"status\":"
    <<std::quoted(accepted?"accepted":"refused")<<",\"stage\":"<<std::quoted(stage)
    <<",\"sdk_build_id\":"<<std::quoted(IRRED_SN_BUILD_ID)<<",\"n\":"<<n
    <<",\"selected_covariance_sha256\":"<<std::quoted(selected_sha)<<",\"source_indices\":";array(indices);
  std::cout<<",\"input_bytes\":"<<input.bytes<<",\"maximum_forward_sensitivity\":"<<sensitivity
    <<",\"arithmetic_policy\":\"F02/binary64-legacy/v1\",\"work\":{\"gaussian_prepare_attempts\":"<<prepare_attempts
    <<",\"profile_prepare_attempts\":"<<profile_prepare_attempts<<",\"evaluate_attempts\":"<<evaluation_attempts
    <<",\"condition_inverse_column_upper_bound\":"<<n*prepare_attempts<<"},\"preparation\":";
  if(preparation_status)std::cout<<"{\"status\":"<<std::quoted(status(*preparation_status))<<",\"numerical_status\":"<<std::quoted(status(*preparation_numerical_status))<<'}';else std::cout<<"null";
  std::cout<<",\"profile_preparation\":";
  if(profile){std::cout<<"{\"status\":"<<std::quoted(status(profile->status()))<<",\"numerical_status\":"<<std::quoted(status(profile->numerical_status()))
    <<",\"payload_available\":"<<(profile_payload?"true":"false")
    <<",\"cached_response_backward_residual\":";if(profile_payload)std::cout<<profile->cached_response_backward_residual();else std::cout<<"null";
    std::cout<<",\"cached_response_forward_sensitivity\":";if(profile_payload)std::cout<<profile->cached_response_forward_sensitivity();else std::cout<<"null";
    std::cout<<",\"gram\":";if(profile_payload)std::cout<<profile->gram();else std::cout<<"null";
    std::cout<<",\"cached_response_solution\":";if(profile_payload)array(profile->cached_response_solution());else std::cout<<"null";
    std::cout<<",\"adjusted_solution_vector\":null,\"native_stationarity\":null}";}else std::cout<<"null";
  std::cout<<",\"rows\":[";
  for(size_t i=0;i<4;++i){if(i)std::cout<<',';const auto& row=rows[i];std::cout<<"{\"case_id\":"<<std::quoted(cases[i])
    <<",\"common_state_sha256\":"<<std::quoted(row.state_sha)<<",\"input_complete\":"<<(row.complete?"true":"false")
    <<",\"attempted\":"<<(row.attempted?"true":"false")<<",\"accepted\":"<<(row.accepted?"true":"false")<<",\"status\":";
    if(row.result)std::cout<<std::quoted(status(row.result->status));else std::cout<<"null";
    std::cout<<",\"numerical_status\":";if(row.result)std::cout<<std::quoted(status(row.result->numerical_status));else std::cout<<"null";
    std::cout<<",\"payload\":";if(row.result&&finite_payload(*row.result,n))payload(*row.result);else std::cout<<"null";std::cout<<'}';}
  std::cout<<"],\"exception\":"<<(error?"true":"false")<<",\"normalized_density\":null,\"prior\":null,\"joint_target\":null}\n";
  return accepted?0:2;
}
