// Original BSD-3-Clause experiment adapter. Physics/numerics stay in Irreducible.
#include "irred/gaussian_design.hpp"
#include <bit>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace irred::statistics;
namespace {
constexpr std::size_t n=3492, p=47, bytes=std::size_t(2)*1024*1024*1024;
std::vector<double> read(const std::string& path,std::size_t count) {
  std::ifstream s(path,std::ios::binary); std::vector<double> v(count);
  s.read(reinterpret_cast<char*>(v.data()),count*sizeof(double));
  if(!s || s.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("exact binary size mismatch");
  return v;
}
void require(bool admitted,const char* message) { if(!admitted) throw std::runtime_error(message); }
}
int main(int argc,char** argv) { try {
  require(argc==2 && std::endian::native==std::endian::little,"one canonical input directory; little-endian host required");
  const std::string dir=argv[1]; auto c=read(dir+"/C.f64",n*n), x=read(dir+"/X.f64",n*p), y=read(dir+"/y.f64",n);
  Metadata m; m.measure="relative full compact Gaussian quadratic; no parameter measure";
  m.table_identity="SH0ES compact 3492 original rows"; m.uncertainty_identity="pinned supplied full C";
  m.ordering_provenance="original FITS row order; no masks"; m.calibration_provenance="released compact constraints retained";
  m.dependence_provenance="supplied full C; no joint-probe independence"; m.source_semantics="released high-level compact products";
  DesignMetadata md; md.residual_unit="released magnitude-like coordinate"; md.design_identity="pinned full X=L.T,3492x47"; md.dependence_identity=m.dependence_provenance;
  for(std::size_t i=0;i<n;++i) m.ordered_ids.push_back("released-row-"+std::to_string(i));
  for(std::size_t j=0;j<p;++j) {
    md.ordered_parameter_ids.push_back("released-parameter-"+std::to_string(j));
    md.parameter_units.push_back(j<37?"mag":j==41||j==43?"mag/dex":j==46?"mag:5log10(H0/[1km/s/Mpc])":"unresolved original coordinate; see lineage");
  }
  auto bound=gaussian_preparation_payload_bound(n,MatrixKind::covariance,irred::numerics::Arithmetic::longdouble_cpu_v1,m);
  require(bound && *bound<=bytes,"Gaussian payload bound");
  auto g=prepare_gaussian(c,MatrixKind::covariance,m,n*n,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  require(g.status()==DensityStatus::finite,"Gaussian admission failed");
  std::vector<double>().swap(c); // Source image buffer no longer needed.
  bound=DesignProfile::preparation_payload_bound(g,p,md); require(bound && *bound<=bytes,"design payload bound");
  DesignPolicy policy{n*n,bytes,1e-10}; auto design=DesignProfile::prepare(std::move(g),x,m.ordered_ids,md,policy);
  require(design.status()==DensityStatus::finite,"full47 design admission failed");
  std::vector<double>().swap(x); // The retained profile owns its design.
  auto fit=design.evaluate(y,m.ordered_ids,policy); require(fit.status==DensityStatus::finite,"fit admission failed");
  std::vector<double> w(p,0); w[46]=1; LinearFunctionalMetadata fm;
  fm.functional_identity="SH0ES released original e46/v1"; fm.output_unit="mag";
  for(std::size_t j=0;j<p;++j) fm.weight_units.push_back(j==46?"mag/mag":"zero weight on original coordinate");
  auto variance=design.estimator_variance(w,md.ordered_parameter_ids,std::move(fm),policy);
  require(variance.status==DensityStatus::finite,"variance admission failed");
  std::cout<<std::setprecision(17)<<"{\"coefficients\":[";
  for(std::size_t j=0;j<p;++j) std::cout<<(j?",":"")<<fit.coefficients[j];
  std::cout<<"],\"quadratic\":"<<fit.quadratic<<",\"relative_log_score\":"<<fit.relative_log_score
    <<",\"variance46\":"<<variance.variance<<",\"variance_sensitivity\":"<<variance.estimated_forward_sensitivity
    <<",\"stationarity\":"<<fit.normalized_normal_equation_residual<<",\"method\":\""<<design.method_id()
    <<"\",\"variance_method\":\""<<variance.method_id<<"\",\"sensitivities\":[";
  // Declared synthetic perturbations of compact constraints, with X,C unchanged.
  // Retain covariance/QR factors; no physical calibration name is inferred here.
  for(std::size_t row=3207;row<=3214;++row) {
    auto yp=y; yp[row]+=0.01; auto plus=design.evaluate(yp,m.ordered_ids,policy);
    yp[row]=y[row]-0.01; auto minus=design.evaluate(yp,m.ordered_ids,policy);
    require(plus.status==DensityStatus::finite && minus.status==DensityStatus::finite,"perturbation admission failed");
    std::cout<<(row==3207?"":",")<<"{\"row\":"<<row<<",\"delta_y\":0.01,\"beta46_plus\":"<<plus.coefficients[46]
      <<",\"beta46_minus\":"<<minus.coefficients[46]<<",\"q_plus\":"<<plus.quadratic<<",\"q_minus\":"<<minus.quadratic<<"}";
  }
  std::cout<<"]}\n"; return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; } }
