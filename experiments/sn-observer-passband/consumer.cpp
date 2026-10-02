#include "irred/background.hpp"
#include "irred/sampled_photometry.hpp"
#include "irred/quantities.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace irred;
int main(int argc,char**argv) {
  try {
    if(argc!=4) throw std::runtime_error("expected row, passband and metre-output paths");
    std::ifstream rows(argv[1]), band(argv[2]);
    std::vector<cosmology::Request> requests;
    std::vector<unsigned> row_ids;
    unsigned id; double zhd,zhel;
    constexpr auto mask=cosmology::Observable::radial|cosmology::Observable::luminosity_shape;
    constexpr auto physical=static_cast<unsigned>(cosmology::Observable::flat_distances_volume);
    while(rows>>id>>zhd>>zhel) {
      row_ids.push_back(id);
      requests.emplace_back(zhd,mask|physical,cosmology::Observer{zhel,cosmology::Convention::released_zhd_zhel},cosmology::PhysicalScale{70});
      requests.emplace_back(zhd,mask|physical,cosmology::Observer{zhd,cosmology::Convention::geometric_same_redshift},cosmology::PhysicalScale{70});
    }
    if(!rows.eof()||row_ids.size()!=321) throw std::runtime_error("source rows");
    for(double z:{0.,1e-10,1e-6})requests.emplace_back(z,mask|physical,cosmology::Observer{z,cosmology::Convention::geometric_same_redshift},cosmology::PhysicalScale{70});
    cosmology::EvaluationPolicy policy;
    policy.integration=numerics::IntegrationPolicy{1e-16,1e-15,100000,30};
    policy.maximum_queries=700;policy.maximum_callbacks=1000000;policy.maximum_segment_visits=100000;policy.maximum_native_bytes=1048576;
    const auto model=cosmology::prepare(cosmology::ConstantQ{0},cosmology::FlatFLRW{});
    const auto result=model.evaluate(requests,policy);
    if(result.status!=cosmology::Status::ok||result.slots.size()!=requests.size())throw std::runtime_error("background batch");
    std::cout<<std::setprecision(17)<<"{\"variant\":\"conditional-sdss-observer-historical-optical/v1\",\"background\":[";
    for(std::size_t i=0;i<result.slots.size();++i){
      const auto&s=result.slots[i];
      if(s.admission_status!=cosmology::Status::ok||!s.radial.value||!s.luminosity_shape.value||!s.physical.value)throw std::runtime_error("background output");
      if(i)std::cout<<',';
      std::cout<<"{\"index\":"<<i<<",\"zHD\":"<<requests[i].z_expansion<<",\"zHEL\":"<<requests[i].observer->redshift<<",\"radial\":"<<s.radial.value->integral<<",\"shape\":"<<*s.luminosity_shape.value<<",\"DL_mpc\":"<<s.physical.value->luminosity_mpc<<'}';
    }
    std::vector<double>w,t;double angstrom,transmission;
    while(band>>angstrom>>transmission){w.push_back(static_cast<double>(static_cast<long double>(angstrom)*1e-10L));t.push_back(transmission);}
    if(!band.eof()||w.size()!=910)throw std::runtime_error("source passband");
    std::ofstream metres(argv[3]);metres<<std::setprecision(17);
    for(std::size_t i=0;i<w.size();++i)metres<<w[i]<<' '<<t[i]<<'\n';
    metres.close();if(!metres)throw std::runtime_error("metre output");
    std::vector<double>sw{1e-7,1.2e-6},sl{1,1};
    photometry::SampledPolicy pp{7,1000,1000};
    std::cout<<"],\"photometry\":[";
    unsigned n=0;
    for(double z:{0.,.1,.06707}){
      const auto pr=photometry::evaluate_sampled({{sw,sl},{w,t},1e20,z,1,1},pp);
      if(pr.admission_status!=numerics::Status::ok||!pr.flux_watt_per_square_metre.value||!pr.energy_joule.value||!pr.expected_photons.value)throw std::runtime_error("photometry output");
      if(n++)std::cout<<',';
      std::cout<<"{\"z\":"<<z<<",\"flux\":"<<*pr.flux_watt_per_square_metre.value<<",\"energy\":"<<*pr.energy_joule.value<<",\"photons\":"<<*pr.expected_photons.value<<'}';
    }
    const std::vector<cosmology::Request> bad{
      {0.1,mask,cosmology::Observer{0.2,cosmology::Convention::geometric_same_redshift}},
      {0.1,mask,cosmology::Observer{-1,cosmology::Convention::released_zhd_zhel}},
      {-0.1,mask,cosmology::Observer{0.1,cosmology::Convention::released_zhd_zhel}}};
    const auto invalid=model.evaluate(bad,policy);
    bool bg=invalid.status==cosmology::Status::ok&&invalid.slots.size()==3;
    if(bg){
      for(unsigned i=0;i<3;++i){
        const auto&s=invalid.slots[i];
        bg&=!s.luminosity_shape.value;
        bg&=s.luminosity_shape.status==(i<2?cosmology::Status::incompatible_convention:cosmology::Status::unsupported_domain);
      }
    }
    auto invalid_t=t;invalid_t[100]=1.1;
    const auto bad_t=photometry::evaluate_sampled({{sw,sl},{w,invalid_t},1e20,0,1,1},pp);
    auto invalid_w=w;invalid_w[100]=invalid_w[99];
    const auto bad_w=photometry::evaluate_sampled({{sw,sl},{invalid_w,t},1e20,0,1,1},pp);
    const bool pt=bad_t.admission_status==numerics::Status::outside_domain&&!bad_t.energy_joule.value;
    const bool pw=bad_w.admission_status==numerics::Status::outside_domain&&!bad_w.energy_joule.value;
    std::cout<<"],\"invalid_controls\":{\"background_no_payload\":"<<(bg?"true":"false")<<",\"transmission_rejected\":"<<(pt?"true":"false")<<",\"wavelength_order_rejected\":"<<(pw?"true":"false")<<"}}\n";
    if(!bg||!pt||!pw)throw std::runtime_error("invalid control admitted");
    return 0;
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
