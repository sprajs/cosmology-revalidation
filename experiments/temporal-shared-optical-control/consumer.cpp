// Ignored external SDK control. No production implementation is changed.
#include "irred/temporal_photometry.hpp"
#include "irred/detector_selection.hpp"
#include "irred/photometry_calibration.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace p = irred::photometry;
namespace d = irred::detector;
using S = irred::numerics::Status;
using W = long double;
using Clock = std::chrono::steady_clock;
constexpr W weights[] = {0.25L, 0.75L};
constexpr W sensitivity = 3.2e-12L;
constexpr W wide_eps = std::numeric_limits<W>::epsilon();
constexpr W binary_eps = std::numeric_limits<double>::epsilon();
const char *status(S s) {
  switch (s) {
  case S::ok: return "ok";
  case S::invalid_input: return "invalid_input";
  case S::nonfinite_input: return "nonfinite_input";
  case S::overflow: return "overflow";
  case S::work_limit: return "work_limit";
  case S::outside_domain: return "outside_domain";
  case S::singular: return "singular";
  case S::not_positive_definite: return "not_positive_definite";
  case S::conditioning_budget_exceeded: return "conditioning_budget_exceeded";
  }
  return "unknown";
}
const char *coverage(p::TemporalCoverage c) {
  switch (c) {
  case p::TemporalCoverage::unassessed: return "unassessed";
  case p::TemporalCoverage::no_overlap: return "no_overlap";
  case p::TemporalCoverage::partial: return "partial";
  case p::TemporalCoverage::full: return "full";
  }
  return "unknown";
}
template<class T> void optional(const std::optional<T> &v) {
  if (v) std::cout << *v; else std::cout << "null";
}
void outcome(const p::Outcome &o) {
  std::cout << "{\"availability\":" << static_cast<unsigned>(o.availability)
            << ",\"status\":\"" << status(o.numerical_status) << "\",\"value\":";
  optional(o.value); std::cout << '}';
}
bool admitted(const p::TemporalBatch &b) {
  if (b.status != S::ok || b.rows.size() != 8) return false;
  for (const auto &r : b.rows)
    if (r.admission_status != S::ok || r.expected_photons.availability != p::Availability::available ||
        r.expected_photons.numerical_status != S::ok || !r.expected_photons.value ||
        !std::isfinite(*r.expected_photons.value) || *r.expected_photons.value < 0) return false;
  return true;
}
void temporal(const p::TemporalBatch &b) {
  std::cout << "{\"status\":\"" << status(b.status) << "\",\"segment_work\":" << b.segment_work
            << ",\"required_rows_admitted\":" << (admitted(b) ? "true" : "false") << ",\"rows\":[";
  for (std::size_t i=0; i<b.rows.size(); ++i) {
    const auto &r=b.rows[i]; if (i) std::cout << ',';
    std::cout << "{\"index\":" << i << ",\"grid_index\":" << r.source.grid_index
              << ",\"band_index\":" << r.source.band_index << ",\"source_epoch_second\":"
              << r.source.source_epoch_observer_second << ",\"observer_lower_second\":" << r.source.observer_lower_second
              << ",\"observer_upper_second\":" << r.source.observer_upper_second
              << ",\"admission_status\":\"" << status(r.admission_status) << "\",\"coverage\":\"" << coverage(r.coverage)
              << "\",\"observer_duration_second\":" << r.observer_duration_second
              << ",\"covered_observer_second\":" << r.covered_observer_second << ",\"covered_fraction\":" << r.covered_fraction
              << ",\"segment_work\":" << r.segment_work << ",\"photons\":";
    outcome(r.expected_photons); std::cout << ",\"mean_flux\":"; outcome(r.mean_flux_watt_per_square_metre);
    std::cout << ",\"energy\":"; outcome(r.energy_joule); std::cout << '}';
  }
  std::cout << "]}";
}
void batch(const d::LikelihoodBatch &b) {
  std::cout << "{\"status\":\"" << status(b.status) << "\",\"photons\":" << b.source.expected_transmitted_photons
            << ",\"QE\":" << b.source.quantum_efficiency << ",\"full_exposure_second\":" << b.source.observer_exposure_second
            << ",\"sigma_electrons\":" << b.source.read_noise_rms_electrons << ",\"poisson_terms\":" << b.poisson_terms
            << ",\"omitted_tail\":" << b.omitted_poisson_tail_estimate << ",\"rows\":[";
  for (std::size_t i=0;i<b.rows.size();++i) {
    const auto &r=b.rows[i]; if(i)std::cout<<',';
    std::cout << "{\"status\":\"" << status(r.status) << "\",\"detected\":" << (r.source.detected ? "true":"false")
              << ",\"threshold_adu\":" << r.source.detection_threshold_adu << ",\"electron_count\":";
    optional(r.source.electron_count); std::cout << ",\"measured_adu\":";optional(r.source.measured_adu);
    std::cout << ",\"zero_probability\":" << (r.zero_probability?"true":"false") << ",\"log_value\":";
    optional(r.log_value);std::cout << ",\"detection_probability\":";optional(r.detection_probability);
    std::cout << ",\"log_error\":" << r.numerical_error_estimate << '}';
  }
  std::cout << "]}";
}
struct Law { S state=S::invalid_input; W record=0,event=0,record_error=0,event_error=0; bool zero_record=false,zero_event=false; };
Law law(const d::LikelihoodBatch &b,std::size_t r) {
  Law v;
  if(b.status!=S::ok){v.state=b.status;return v;}
  if(b.rows.size()!=6)return v;
  const auto &data=b.rows[r], &probe=b.rows[3+r];
  if(data.status!=S::ok){v.state=data.status;return v;}
  if(probe.status!=S::ok){v.state=probe.status;return v;}
  if(!probe.detection_probability)return v;
  v.zero_record=data.zero_probability;
  if(!v.zero_record){if(!data.log_value)return v;v.record=*data.log_value;v.record_error=data.numerical_error_estimate+4*binary_eps*(1+std::abs(v.record));}
  const W q=*probe.detection_probability;v.zero_event=q==0;
  if(!v.zero_event){
    v.event=std::log(q);
    W error=b.omitted_poisson_tail_estimate+128*wide_eps*b.poisson_terms;
    if(!probe.zero_probability){if(!probe.log_value)return v;error=std::max(error,std::exp(W(*probe.log_value))*probe.numerical_error_estimate);}
    v.event_error=error/q+4*binary_eps*(1+std::abs(v.event));
  }
  v.state=S::ok;return v;
}
struct Channel { d::Moments moments; std::array<d::LikelihoodBatch,3> attempts; std::array<bool,3> present{true,false,false}; };
Law sensitive(const Channel &c,std::size_t r) {
  auto v=law(c.attempts[0],r);if(v.state!=S::ok)return v;
  W er=0,ee=0;
  for(unsigned a=1;a<3;++a){if(!c.present[a])continue;const auto e=law(c.attempts[a],r);if(e.state!=S::ok){v.state=e.state;return v;}
    if(e.zero_record!=v.zero_record||e.zero_event!=v.zero_event){v.state=S::conditioning_budget_exceeded;return v;}
    if(!v.zero_record)er=std::max(er,std::abs(e.record-v.record)+e.record_error);
    if(!v.zero_event)ee=std::max(ee,std::abs(e.event-v.event)+e.event_error);
  }
  v.record_error+=er;v.event_error+=ee;return v;
}
struct Mix { W log=0,error=0;bool zero=true; };
Mix mixture(const std::array<Law,2>&v,bool event) {
  Mix m;W total=0,relative_error=0,peak=-INFINITY;
  for(unsigned s=0;s<2;++s)if(!(event?v[s].zero_event:v[s].zero_record))peak=std::max(peak,std::log(weights[s])+(event?v[s].event:v[s].record));
  if(!std::isfinite(peak))return m;
  for(unsigned s=0;s<2;++s)if(!(event?v[s].zero_event:v[s].zero_record)){
    const W q=std::exp(std::log(weights[s])+(event?v[s].event:v[s].record)-peak);
    total+=q;relative_error+=q*std::expm1(event?v[s].event_error:v[s].record_error);
  }
  m.zero=false;m.log=peak+std::log(total);
  const W relative=relative_error/total+18*wide_eps;
  m.error=relative<1?-std::log1p(-relative):INFINITY;return m;
}
bool log_gate(W value,W error){return std::isfinite(value)&&std::isfinite(error)&&error<=5e-12L+1e-8L*(1+std::abs(value));}
bool emit_mix(const std::array<Channel,8>&channels,unsigned r,bool selected_only,
              const Channel* replacement=nullptr,unsigned replacement_index=7) {
  const auto channel=[&](unsigned i)->const Channel&{return replacement&&i==replacement_index?*replacement:channels[i];};
  std::array<Law,2> state;bool ok=true;S result_status=S::ok;
  std::cout << "{\"record_index\":"<<r<<",\"conditional\":[";
  for(unsigned s=0;s<2;++s){state[s].state=S::ok;
    for(unsigned c=0;c<4;++c){auto v=sensitive(channel(4*s+c),r);if(v.state!=S::ok){state[s].state=v.state;result_status=v.state;ok=false;continue;}
      state[s].record+=v.record;state[s].event+=v.event;state[s].record_error+=v.record_error;state[s].event_error+=v.event_error;
      state[s].zero_record|=v.zero_record;state[s].zero_event|=v.zero_event;
    }
    if(s)std::cout<<',';
    std::cout<<"{\"state_id\":\"S"<<s<<"\",\"mass\":"<<weights[s]<<",\"status\":\""<<status(state[s].state)<<"\",\"log_record\":"<<state[s].record
             <<",\"log_event\":"<<state[s].event<<",\"record_error\":"<<state[s].record_error<<",\"event_error\":"<<state[s].event_error<<'}';
  }
  const auto record=mixture(state,false),event=mixture(state,true);
  const W value=record.log-(selected_only?event.log:0),error=record.error+(selected_only?event.error:0);
  ok=ok&&!record.zero&&!event.zero&&log_gate(record.log,record.error)&&log_gate(event.log,event.error)&&log_gate(value,error);
  if(!ok&&result_status==S::ok)result_status=S::conditioning_budget_exceeded;
  if(selected_only&&r==2){ok=false;result_status=S::invalid_input;}
  std::cout<<"],\"selected_only\":"<<(selected_only?"true":"false")<<",\"status\":\""<<status(result_status)<<"\",\"log_joint\":";
  if(ok)std::cout<<record.log;else std::cout<<"null";
  std::cout<<",\"log_event\":";if(ok)std::cout<<event.log;else std::cout<<"null";
  std::cout<<",\"log_value\":";if(ok)std::cout<<value;else std::cout<<"null";
  std::cout<<",\"joint_log_error\":"<<record.error<<",\"event_log_error\":"<<event.error<<",\"value_log_error\":"<<error;
  if(!ok){std::cout<<",\"product_per_exposure_log_record\":null,\"product_per_exposure_log_event\":null,\"record_relative_difference\":null,\"event_relative_difference\":null}";return false;}
  W independent_record=0,independent_event=0;
  for(unsigned e=0;e<2;++e){std::array<Law,2> local;for(unsigned s=0;s<2;++s)for(unsigned b=0;b<2;++b){const auto v=sensitive(channel(4*s+2*e+b),r);local[s].record+=v.record;local[s].event+=v.event;local[s].zero_record|=v.zero_record;local[s].zero_event|=v.zero_event;}
    independent_record+=mixture(local,false).log;independent_event+=mixture(local,true).log;
  }
  std::cout<<",\"product_per_exposure_log_record\":"<<independent_record<<",\"product_per_exposure_log_event\":"<<independent_event
           <<",\"record_relative_difference\":"<<std::expm1(record.log-independent_record)<<",\"event_relative_difference\":"<<std::expm1(event.log-independent_event)<<'}';return true;
}
int main() {
  std::cout<<std::setprecision(21);
  const std::array<double,3> time{-1,0,2};const std::array<double,2>wavelength{1,2};
  const double C=0x1p-80;
  const std::array<double,6>luminosity{.75*C,1.1875*C,C,1.5*C,1.5*C,2.125*C};
  const p::TemporalGrid grid{"synthetic-bilinear-ux16-v1","original analytic synthetic control",time,wavelength,luminosity};
  const std::array<double,2>blue{2,3},red{3,4},low{.25,.25},high{.75,.75};
  const std::array<p::TemporalBand,4>bands{{{"S0/blue","supplied synthetic optical law",{blue,low}},
    {"S0/red","supplied synthetic optical law",{red,high}},{"S1/blue","supplied synthetic optical law",{blue,high}},
    {"S1/red","supplied synthetic optical law",{red,low}}}};
  const auto start=Clock::now();
  auto owner=p::prepare_temporal(std::span(&grid,1),bands,{1,4,64,6,1048576});
  const auto prepared=Clock::now();
  std::array<p::TemporalExposure,8>rows;
  for(unsigned s=0;s<2;++s)for(unsigned e=0;e<2;++e)for(unsigned b=0;b<2;++b)
    rows[4*s+2*e+b]={0,2*s+b,1,1,1,10,e?12.:9.,e?15.:11.};
  const p::TemporalPolicy policy{p::transmitted_photons,8,1048576,512};
  const auto native=p::evaluate_temporal(owner,rows,policy);const auto evaluated=Clock::now();
  std::cout<<"{\"schema\":\"external-temporal-optical-sdk-control/v1\",\"model_id\":\""<<p::temporal_model_id
           <<"\",\"detector_model_id\":\""<<d::model_id<<"\",\"selection_id\":\""<<d::selection_id
           <<"\",\"sampled_empirical_sensitivity\":"<<p::calibration_sampled_relative_sensitivity
           <<",\"temporal_photon_endpoint_sensitivity\":"<<sensitivity<<",\"prepared_status\":\""<<status(owner.status())
           <<"\",\"grid_count\":"<<owner.grid_count()<<",\"band_count\":"<<owner.band_count()<<",\"retained_temporal_payload_bytes\":";
  optional(owner.retained_payload_bytes());std::cout<<",\"prepare_seconds\":"<<std::chrono::duration<double>(prepared-start).count()
    <<",\"temporal_evaluate_seconds\":"<<std::chrono::duration<double>(evaluated-prepared).count()<<",\"primary\":";
  temporal(native);std::cout<<",\"detector_controls\":[";
  std::size_t all_terms=0;
  std::array<std::array<Channel,8>,2>retained;
  if(admitted(native))for(unsigned noise=0;noise<2;++noise){if(noise)std::cout<<',';
    auto &channels=retained[noise];
    std::cout<<"{\"sigma_electrons\":"<<(noise?.75:0)<<",\"attempts\":[";
    for(unsigned i=0;i<8;++i){const unsigned c=i%4,b=c%2;const double n=*native.rows[i].expected_photons.value;
      d::Input input{d::PhotonLaw::poisson_arrivals,n,b?.5:.75,.25,.125,c<2?2.:3.,noise?.75:0.,2,1};
      channels[i].moments=d::moments(input);
      std::array<d::Observation,6>obs;
      for(unsigned r=0;r<3;++r){obs[r].detection_threshold_adu=1.5;obs[r].detected=!(r==2&&c==3);
        if(obs[r].detected){if(noise)obs[r].measured_adu=(c==0||c==3)?3.:3.5;else obs[r].electron_count=(c==0||c==3)?2:3;}
        obs[3+r].detection_threshold_adu=1.5;}
      if(i)std::cout<<',';
      std::cout<<"{\"state_index\":"<<i/4<<",\"channel_index\":"<<c<<",\"lambda_electrons\":";
      optional(channels[i].moments.poisson_mean_electrons);std::cout<<",\"batches\":[";
      for(unsigned a=0;a<3;++a){if(a)std::cout<<',';if(a&&n==0){std::cout<<"null";continue;}channels[i].present[a]=true;
        if(a)input.expected_transmitted_photons=std::nextafter(double(W(n)*(1+(a==1?-sensitivity:sensitivity))),a==1?0.:INFINITY);
        if(a&&(!std::isnormal(input.expected_transmitted_photons)||input.expected_transmitted_photons<=0)){
          channels[i].attempts[a].source=input;channels[i].attempts[a].status=S::outside_domain;
        }else channels[i].attempts[a]=d::likelihood(input,obs,{6,1048576,256});
        all_terms+=channels[i].attempts[a].poisson_terms;batch(channels[i].attempts[a]);}
      std::cout<<"]}";
    }
    std::cout<<"],\"records\":[";for(unsigned r=0;r<3;++r){if(r)std::cout<<',';emit_mix(channels,r,r==1);}std::cout<<"]}";
  }
  std::cout<<"],\"primary_poisson_terms\":"<<all_terms<<",\"actual_temporal_controls\":[";
  std::size_t all_temporal_work=native.segment_work;
  const char*ids[]={"no-time-overlap","zero-duration","invalid-redshift","required-state-row-refusal","temporal-work-refusal","temporal-row-admission-refusal"};
  for(unsigned k=0;k<6;++k){if(k)std::cout<<',';auto changed=rows;auto cp=policy;
    if(k==0){changed[0].observer_lower_second=20;changed[0].observer_upper_second=22;}
    if(k==1)changed[0].observer_upper_second=changed[0].observer_lower_second;
    if(k==2)changed[0].redshift=-1;
    if(k==3)changed[7].band_index=4;
    if(k==4)cp.maximum_segment_work=0;
    if(k==5)cp.maximum_rows=7;
    const auto result=p::evaluate_temporal(owner,changed,cp);
    all_temporal_work+=result.segment_work;
    std::cout<<"{\"id\":\""<<ids[k]<<"\",\"joint_withheld\":"<<(!admitted(result)?"true":"false")<<",\"native\":";temporal(result);std::cout<<'}';
  }
  std::cout<<"],\"actual_detector_controls\":[";
  const char*did[]={"detector-QE-domain","detector-work-refusal","selected-nondetection","zero-noise-threshold-boundary"};
  for(unsigned k=0;k<4;++k){if(k)std::cout<<',';
    d::Input input{d::PhotonLaw::poisson_arrivals,1,.75,.25,.125,2,.75,2,1};
    d::Observation obs{1.5,true,{},3.,d::SelectionMeasure::joint_detection_record};d::SelectionPolicy dp{1,1048576,256};
    if(k==0)input.quantum_efficiency=1.01;
    if(k==1)dp.maximum_poisson_terms=1;
    if(k==2){obs.detected=false;obs.measured_adu.reset();obs.measure=d::SelectionMeasure::selected_only;}
    if(k==3){input.read_noise_rms_electrons=0;obs.measured_adu.reset();obs.electron_count=0;}
    const auto result=d::likelihood(input,std::span(&obs,1),dp);
    all_terms+=result.poisson_terms;
    std::cout<<"{\"id\":\""<<did[k]<<"\",\"native\":";batch(result);std::cout<<'}';
  }
  std::cout<<"],\"actual_external_reducer_controls\":[";
  if(admitted(native))for(unsigned k=0;k<3;++k){if(k)std::cout<<',';
    const auto &channels=retained[1];
    std::cout<<"{\"id\":\""<<(k==0?"required-QE-refusal":k==1?"required-work-refusal":"selected-censored-refusal")<<"\",\"retained_state_masses\":[0.25,0.75],\"replacement_batch\":";
    Channel replacement;
    if(k<2){auto input=channels[7].attempts[0].source;std::array<d::Observation,6>obs;
      for(unsigned j=0;j<6;++j)obs[j]=channels[7].attempts[0].rows[j].source;
      if(k==0)input.quantum_efficiency=1.01;
      replacement.attempts[0]=d::likelihood(input,obs,{6,1048576,k==1?1u:256u});
      all_terms+=replacement.attempts[0].poisson_terms;batch(replacement.attempts[0]);
    }else std::cout<<"null";
    std::cout<<",\"reducer\":";emit_mix(channels,k==2?2:0,k==2,k<2?&replacement:nullptr);std::cout<<'}';
  }
  std::cout<<"],\"all_poisson_terms\":"<<all_terms<<",\"all_temporal_segment_work\":"<<all_temporal_work
           <<",\"detector_six_row_payload_bound_bytes\":";optional(d::likelihood_payload_bound(6));
  std::cout<<",\"total_native_seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<"}\n";
  return owner.status()==S::ok&&admitted(native)?0:1;
}
