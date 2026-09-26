/* Private conditional-estimator change; disabled unless exact env value 1.
   No flux, fitted residual, or truth value enters the domain construction. */
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PROSP_DOMAIN_NEPOCH 117
#define PROSP_DOMAIN_NCID 64

typedef struct {
 char cid[64]; int configured,n,iterations; double mjd[PROSP_DOMAIN_NEPOCH];
 double z,off,tlo,thi,klo,khi,rawlo,rawhi,lo,hi,slo,shi,zmax;
 unsigned long long means,fcns;
 double phase_lo,phase_hi,edge_phase,edge_peak,peak_lo,peak_hi;
} PROSP_DOMAIN_STATE;
static PROSP_DOMAIN_STATE prosp_domain_states[PROSP_DOMAIN_NCID];
static int prosp_domain_count=0,prosp_domain_mode=-1;
static void prosp_domain_finish(void) {
 int i; for(i=0;i<prosp_domain_count;i++) {
  PROSP_DOMAIN_STATE *s=&prosp_domain_states[i];
  printf("PROSP_DOMAIN_SUMMARY %s %d %d %llu %llu %.17g %.17g %.17g %.17g %.17g %.17g\n",
   s->cid,s->configured,s->iterations,s->means,s->fcns,
   s->phase_lo,s->phase_hi,s->edge_phase,s->peak_lo,s->peak_hi,s->edge_peak);
 } fflush(stdout);
}
static void prosp_domain_abort(const char *reason,const char *cid) {
 fprintf(stderr,"PROSP_DOMAIN_ABORT %s %s\n",cid,reason);
 fflush(stdout);fflush(stderr);exit(86);
}
static int prosp_domain_on(void) {
 if(prosp_domain_mode<0) {
  const char *s=getenv("PROSP_HARD_PEAK_DOMAIN");
  if(!s || strcmp(s,"0")==0) prosp_domain_mode=0;
  else if(strcmp(s,"1")==0) {prosp_domain_mode=1;atexit(prosp_domain_finish);}
  else prosp_domain_abort("environment_must_be_absent_0_or_1","GLOBAL");
 }
 return prosp_domain_mode;
}
static PROSP_DOMAIN_STATE *prosp_domain_state(const char *cid) {
 int i;
 if(strlen(cid)==0 || strlen(cid)>=64) prosp_domain_abort("invalid_CID",cid);
 for(i=0;i<prosp_domain_count;i++)
  if(strcmp(prosp_domain_states[i].cid,cid)==0) return &prosp_domain_states[i];
 if(prosp_domain_count==PROSP_DOMAIN_NCID) prosp_domain_abort("CID_capacity",cid);
 PROSP_DOMAIN_STATE *s=&prosp_domain_states[prosp_domain_count++];
 strcpy(s->cid,cid);s->phase_lo=s->edge_phase=s->edge_peak=s->peak_lo=DBL_MAX;
 s->phase_hi=s->peak_hi=-DBL_MAX;return s;
}
static void prosp_domain_tables(const char *cid,double klo,double khi,
 double *tlo,double *thi,double *slo,double *shi) {
 int nt=SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_TREST];
 int ns=SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_SHAPEPAR];
 if(nt<2 || ns<2) prosp_domain_abort("SNooPy_table_not_loaded",cid);
 *tlo=fmax(SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][1],klo);
 *thi=fmin(SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][nt],khi);
 *slo=SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][1];
 *shi=SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][ns];
 if(!isfinite(klo)||!isfinite(khi)||klo>=khi||!isfinite(*tlo)||
    !isfinite(*thi)||*tlo>=*thi||!isfinite(*slo)||!isfinite(*shi)||*slo>=*shi)
  prosp_domain_abort("invalid_table_bounds",cid);
}
/* Matches R8 FCNSNLC ordering, then separately R4 LOAD_EPALL cut arithmetic. */
static int prosp_domain_endpoint_ok(int n,const double *mjd,double z,double off,
 double p,double tlo,double thi) {
 int i; for(i=0;i<n;i++) {
  double phase=((mjd[i]-off)-p)/(1.0+z);
  float tobs=(float)(mjd[i]-(p+off));
  float z1=1.0f+(float)z;
  float mask_phase=tobs/z1;
  if(!isfinite(phase)||phase<tlo||phase>thi||
     !isfinite(mask_phase)||mask_phase<tlo||mask_phase>thi) return 0;
 }return 1;
}
void prosp_peak_bounds__(char *cid,int *iter,int *n,double *mjd,
 double *z,double *zstep,double *off,double *initial,
 double *klo,double *khi,double *zmax,int *model,int *snoopy,
 double *lower,double *upper) {
 if(!prosp_domain_on()) return;
 int i,j; double tlo,thi,slo,shi,rawlo=-DBL_MAX,rawhi=DBL_MAX,lo,hi;
 PROSP_DOMAIN_STATE *s=prosp_domain_state(cid);
 if(*model!=*snoopy || *n!=PROSP_DOMAIN_NEPOCH || *zstep!=0.0 ||
    !isfinite(*z)||*z<0||!isfinite(*zmax)||*z>*zmax||
    !isfinite(*off)||!isfinite(*initial)||*iter<1)
  prosp_domain_abort("unsupported_model_metadata_or_free_redshift",cid);
 prosp_domain_tables(cid,*klo,*khi,&tlo,&thi,&slo,&shi);
 for(i=0;i<*n;i++) {
  if(!isfinite(mjd[i])) prosp_domain_abort("nonfinite_MJD",cid);
  rawlo=fmax(rawlo,mjd[i]-(1.0+*z)*thi);
  rawhi=fmin(rawhi,mjd[i]-(1.0+*z)*tlo);
 }
 lo=rawlo-*off;hi=rawhi-*off;
 if(!isfinite(lo)||!isfinite(hi)||lo>=hi) prosp_domain_abort("empty_domain",cid);
 for(j=0;j<1024 && !prosp_domain_endpoint_ok(*n,mjd,*z,*off,lo,tlo,thi);j++)
  lo=nextafter(lo,hi);
 if(j==1024) prosp_domain_abort("lower_roundoff_support_failure",cid);
 for(j=0;j<1024 && !prosp_domain_endpoint_ok(*n,mjd,*z,*off,hi,tlo,thi);j++)
  hi=nextafter(hi,lo);
 if(j==1024 || lo>=hi) prosp_domain_abort("upper_roundoff_support_failure",cid);
 if(*initial<lo||*initial>hi) prosp_domain_abort("initial_peak_outside_domain",cid);
 if(s->configured) {
  if(s->n!=*n||s->z!=*z||s->off!=*off||s->tlo!=tlo||s->thi!=thi||
     s->klo!=*klo||s->khi!=*khi||s->slo!=slo||s->shi!=shi||
     s->zmax!=*zmax||s->lo!=lo||s->hi!=hi||*lower!=lo||*upper!=hi)
   prosp_domain_abort("domain_or_native_bounds_changed_between_iterations",cid);
  for(i=0;i<*n;i++) if(s->mjd[i]!=mjd[i])
   prosp_domain_abort("MJD_order_or_value_changed_between_iterations",cid);
 } else {
  s->configured=1;s->n=*n;s->z=*z;s->off=*off;s->tlo=tlo;s->thi=thi;
  s->klo=*klo;s->khi=*khi;s->rawlo=rawlo;s->rawhi=rawhi;
  s->lo=lo;s->hi=hi;s->slo=slo;s->shi=shi;s->zmax=*zmax;
  for(i=0;i<*n;i++) s->mjd[i]=mjd[i];
 }
 *lower=lo;*upper=hi;s->iterations++;
 printf("PROSP_PEAK_BOUNDS %s %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
  cid,*iter,*n,*z,*off,tlo,thi,*klo,*khi,rawlo,rawhi,
  lo+*off,hi+*off,lo,hi,*initial+*off,*initial-lo,hi-*initial);
 fflush(stdout);
}
/* Called only after native FCN prior-only/out-of-bound early returns, before
   its first physical mean. No likelihood penalty or coordinate clipping. */
void prosp_peak_fcn_guard__(char *cid,int *iter,double *peak,double *off,double *z) {
 if(!prosp_domain_on()) return;
 PROSP_DOMAIN_STATE *s=prosp_domain_state(cid);
 if(!s->configured||!isfinite(*peak)||!isfinite(*off)||!isfinite(*z)||
    *off!=s->off||*z!=s->z||*peak<s->lo||*peak>s->hi) {
  fprintf(stderr,"PROSP_DOMAIN_FCN_REJECT %s %d %.17g %.17g %.17g\n",cid,*iter,*peak,*off,*z);
  prosp_domain_abort("FCN_peak_or_redshift_outside_fixed_domain",cid);
 }
 double edge=fmin(*peak-s->lo,s->hi-*peak);s->fcns++;
 s->edge_peak=fmin(s->edge_peak,edge);
 s->peak_lo=fmin(s->peak_lo,*peak+*off);s->peak_hi=fmax(s->peak_hi,*peak+*off);
}
/* Unconditional in enabled mode, even without PROSP_FIT_SUPPORT. Called
   before nearest-filter, template, extinction or KCOR operations in USRFUN.
   Direct auxiliary rest-frame calls are checked in their actual phase/z. */
void prosp_physical_mean_guard__(char *cid,int *iter,int *filt,int *epoch,
 double *trest,double *tobs,double *z,double *D,double *shape,double *AV,double *RV,
 double *klo,double *khi,double *zmax,int *model,int *snoopy) {
 if(!prosp_domain_on()) return;
 double tlo,thi,slo,shi;
 if(*model!=*snoopy) prosp_domain_abort("non_SNooPy_mean",cid);
 prosp_domain_tables(cid,*klo,*khi,&tlo,&thi,&slo,&shi);
 if(!isfinite(*trest)||!isfinite(*tobs)||!isfinite(*z)||!isfinite(*D)||
    !isfinite(*shape)||!isfinite(*AV)||!isfinite(*RV)||!isfinite(*zmax)||
    *trest<tlo||*trest>thi||*shape<slo||*shape>shi||*z<0||*z>*zmax) {
  fprintf(stderr,"PROSP_DOMAIN_MEAN_REJECT %s %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
   cid,*iter,*filt,*epoch,*trest,*tobs,*z,*D,*shape,*AV,*RV,tlo,thi,slo,shi);
  prosp_domain_abort("physical_mean_outside_table_domain",cid);
 }
 PROSP_DOMAIN_STATE *s=prosp_domain_state(cid);s->means++;
 s->phase_lo=fmin(s->phase_lo,*trest);s->phase_hi=fmax(s->phase_hi,*trest);
 s->edge_phase=fmin(s->edge_phase,fmin(*trest-tlo,thi-*trest));
}
