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
std::vector<double> select_zero(std::span<const double> x,std::size_t rows,std::size_t columns,std::size_t fixed) {
  require(fixed<columns && x.size()==rows*columns,"fixed-coordinate selection dimensions");
  std::vector<double> selected; selected.reserve(rows*(columns-1));
  for(std::size_t i=0;i<rows;++i) for(std::size_t j=0;j<columns;++j) if(j!=fixed) selected.push_back(x[i*columns+j]);
  return selected;
}
void control() {
  Metadata m; m.measure="synthetic exact fixed-zero control";m.ordered_ids={"a","b","c"};
  m.table_identity="synthetic";m.uncertainty_identity="identity";m.ordering_provenance="declared";
  m.calibration_provenance="none synthetic";m.dependence_provenance="identity";m.source_semantics="synthetic";
  const std::vector<double> identity={1,0,0,0,1,0,0,0,1};
  auto g=prepare_gaussian(identity,MatrixKind::covariance,m,9,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  DesignMetadata md;md.ordered_parameter_ids={"first","second"};md.parameter_units={"u","u"};
  md.residual_unit="u";md.design_identity="exact zero-coordinate control";md.dependence_identity="identity";
  auto x=select_zero(identity,3,3,2);auto design=DesignProfile::prepare(std::move(g),x,m.ordered_ids,md);
  const std::vector<double> y={2,3,7};auto fit=design.evaluate(y,m.ordered_ids);
  require(fit.status==DensityStatus::finite && fit.coefficients==std::vector<double>{2,3} && fit.quadratic==49,"exact fixed-zero linear control");
  auto wrong=m.ordered_ids;std::swap(wrong[0],wrong[1]);
  require(design.evaluate(y,wrong).status!=DensityStatus::finite,"invalid row order admitted");
  bool refused=false;try {select_zero(identity,3,3,3);}catch(const std::exception&) {refused=true;}
  require(refused,"invalid fixed support admitted");
  std::cout<<"{\"exact_fixed_zero\":true,\"invalid_order_refused\":true,\"invalid_support_refused\":true}\n";
}
}
int main(int argc,char** argv) { try {
  if(argc==2 && std::string(argv[1])=="--self-test") {control();return 0;}
  require((argc==2 || argc==3) && std::endian::native==std::endian::little,"canonical input directory and optional released-fixed44/v1; little-endian host required");
  const bool fixed=argc==3;
  require(!fixed || std::string(argv[2])=="released-fixed44/v1","unknown target");
  const std::string dir=argv[1]; auto c=read(dir+"/C.f64",n*n), x=read(dir+"/X.f64",n*p), y=read(dir+"/y.f64",n);
  std::vector<std::size_t> active;
  for(std::size_t j=0;j<p;++j) if(!fixed || j!=44) active.push_back(j);
  if(fixed) {
    // Literal-zero coordinate contributes exactly zero; preserve every row/C.
    x=select_zero(x,n,p,44);
  }
  Metadata m; m.measure="relative full compact Gaussian quadratic; no parameter measure";
  if(fixed) m.measure="relative fixed44=0 compact Gaussian quadratic; original active46 coordinates";
  m.table_identity="SH0ES compact 3492 original rows"; m.uncertainty_identity="pinned supplied full C";
  m.ordering_provenance="original FITS row order; no masks"; m.calibration_provenance="released compact constraints retained";
  m.dependence_provenance="supplied full C; no joint-probe independence"; m.source_semantics="released high-level compact products";
  DesignMetadata md; md.residual_unit="released magnitude-like coordinate"; md.design_identity="pinned full X=L.T,3492x47"; md.dependence_identity=m.dependence_provenance;
  if(fixed) md.design_identity="pinned X=L.T active original axes0..43,45,46;3492x46; beta44=literal0";
  for(std::size_t i=0;i<n;++i) m.ordered_ids.push_back("released-row-"+std::to_string(i));
  for(auto j:active) {
    md.ordered_parameter_ids.push_back("released-parameter-"+std::to_string(j));
    md.parameter_units.push_back(j<37 || (fixed && (j==37||j==38||j==39||j==40||j==42||j==45))?"mag":j==41||j==43?"mag/dex":j==46?"mag:5log10(H0/[1km/s/Mpc])":"unresolved original coordinate; see lineage");
  }
  auto bound=gaussian_preparation_payload_bound(n,MatrixKind::covariance,irred::numerics::Arithmetic::longdouble_cpu_v1,m);
  require(bound && *bound<=bytes,"Gaussian payload bound");
  auto g=prepare_gaussian(c,MatrixKind::covariance,m,n*n,1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  require(g.status()==DensityStatus::finite,"Gaussian admission failed");
  std::vector<double>().swap(c); // Source image buffer no longer needed.
  bound=DesignProfile::preparation_payload_bound(g,active.size(),md); require(bound && *bound<=bytes,"design payload bound");
  DesignPolicy policy{n*n,bytes,1e-10}; auto design=DesignProfile::prepare(std::move(g),x,m.ordered_ids,md,policy);
  require(design.status()==DensityStatus::finite,"selected design admission failed");
  std::vector<double>().swap(x); // The retained profile owns its design.
  auto fit=design.evaluate(y,m.ordered_ids,policy); require(fit.status==DensityStatus::finite,"fit admission failed");
  const auto contrast=active.size()-1;
  std::vector<double> w(active.size(),0); w[contrast]=1; LinearFunctionalMetadata fm;
  fm.functional_identity="SH0ES released original e46/v1"; fm.output_unit="mag";
  for(auto j:active) fm.weight_units.push_back(j==46?"mag/mag":"zero weight on original coordinate");
  auto variance=design.estimator_variance(w,md.ordered_parameter_ids,std::move(fm),policy);
  require(variance.status==DensityStatus::finite,"variance admission failed");
  std::cout<<std::setprecision(17)<<"{\"coefficients\":[";
  std::size_t k=0;
  for(std::size_t j=0;j<p;++j) std::cout<<(j?",":"")<<(fixed&&j==44?0:fit.coefficients[k++]);
  std::cout<<"],\"quadratic\":"<<fit.quadratic<<",\"relative_log_score\":"<<fit.relative_log_score
    <<",\"variance46\":"<<variance.variance<<",\"variance_sensitivity\":"<<variance.estimated_forward_sensitivity
    <<",\"stationarity\":"<<fit.normalized_normal_equation_residual<<",\"rank\":"<<active.size()<<",\"method\":\""<<design.method_id()
    <<"\",\"variance_method\":\""<<variance.method_id<<"\",\"sensitivities\":[";
  // Historical controls use synthetic .01; fixed46 uses pinned printed source
  // calibration-constraint mean shifts. Both retain X,C and prepared factors.
  const std::vector<std::size_t> rows=fixed?std::vector<std::size_t>{3211,3213,3214}:std::vector<std::size_t>{3207,3208,3209,3210,3211,3212,3213,3214};
  for(std::size_t ri=0;ri<rows.size();++ri) {
    const auto row=rows[ri]; const double delta=fixed?(row==3211?.10:row==3213?.032:.0263):.01;
    auto yp=y; yp[row]+=delta; auto plus=design.evaluate(yp,m.ordered_ids,policy);
    yp[row]=y[row]-delta; auto minus=design.evaluate(yp,m.ordered_ids,policy);
    require(plus.status==DensityStatus::finite && minus.status==DensityStatus::finite,"perturbation admission failed");
    std::cout<<(ri?",":"")<<"{\"row\":"<<row<<",\"delta_y\":"<<delta<<",\"beta46_plus\":"<<plus.coefficients[contrast]
      <<",\"beta46_minus\":"<<minus.coefficients[contrast]<<",\"q_plus\":"<<plus.quadratic<<",\"q_minus\":"<<minus.quadratic<<"}";
  }
  std::cout<<"]}\n"; return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; } }
