/* Opt-in numerical start variation only. Shared INIVAL is const here;
   the only writable argument is a separate local MNPARM start value. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
static int prosp_mn_shift_initialized=0;
static double prosp_mn_shift=0.0;
static void prosp_mn_shift_abort(int cid,const char *reason) {
 fprintf(stderr,"PROSP_MNPARM_ABORT %d %s\n",cid,reason);
 fflush(stdout);fflush(stderr);exit(87);
}
void prosp_minuit_peak_start__(const int *cid,const double *iteration,
 const char *iteration_name,const double *source_initial,
 const double *step,const double *lower,const double *upper,double *local_start,int *readback_enabled) {
 if(!prosp_mn_shift_initialized) {
  const char *raw=getenv("PROSP_MINUIT_PEAK_SHIFT");
  if(raw) {
   char *end;double d=strtod(raw,&end);
   if(end==raw || *end!='\0' || !isfinite(d) || (d!=0.0 && fabs(d)!=2.0))
    prosp_mn_shift_abort(*cid,"shift_must_be_0_or_plus_minus_2_days");
   prosp_mn_shift=d;
  }
  prosp_mn_shift_initialized=1;
 }
 /* Exact default path: no arguments changed, checks or records. */
 if(prosp_mn_shift==0.0) return;
 const char *hard=getenv("PROSP_HARD_PEAK_DOMAIN");
 if(!hard || strcmp(hard,"1")!=0)
  prosp_mn_shift_abort(*cid,"requires_enabled_fixed_metadata_domain");
 if(strncmp(iteration_name,"ITER",4)!=0 ||
    (iteration_name[4]!=' ' && iteration_name[4]!='\0') ||
    !isfinite(*iteration) || *iteration<1 || *iteration>12 ||
    *iteration!=floor(*iteration) || !isfinite(*source_initial) ||
    !isfinite(*step) || !isfinite(*lower) || !isfinite(*upper) ||
    *lower>=*upper || *local_start!=*source_initial)
  prosp_mn_shift_abort(*cid,"invalid_native_peak_start_context");
 double applied=(*iteration==1.0)?prosp_mn_shift:0.0;
 if(applied!=0.0 && *step<=0.0)
  prosp_mn_shift_abort(*cid,"cannot_shift_fixed_peak");
 double proposal=*source_initial+applied;
 if(!isfinite(proposal)||proposal<*lower||proposal>*upper)
  prosp_mn_shift_abort(*cid,"MNPARM_peak_start_outside_fixed_domain");
 /* No write in later iterations; all other parameters bypass this helper. */
 if(applied!=0.0) *local_start=proposal;
 *readback_enabled=1;
 printf("PROSP_MNPARM_PEAK %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
  *cid,(int)*iteration,prosp_mn_shift,applied,*source_initial,*local_start,
  *source_initial,*lower,*upper,*local_start-*lower,*upper-*local_start);
 fflush(stdout);
}

/* MNPOUT reads the actual stored MINUIT external coordinate; output buffers
   are separate locals, so this cannot overwrite INIVAL or PARNAME. */
void prosp_minuit_peak_stored__(const int *cid,const double *iteration,
 const double *source_initial,const double *proposed,const double *stored,
 const double *lower,const double *upper,const double *stored_lower,
 const double *stored_upper,const int *internal_index,const char *name) {
 if(prosp_mn_shift==0.0) return;
 if(!isfinite(*stored)||!isfinite(*stored_lower)||!isfinite(*stored_upper)||
    *stored<*lower||*stored>*upper||*stored_lower!=*lower||
    *stored_upper!=*upper||*internal_index<1||strcmp(name,"PKMJD")!=0)
  prosp_mn_shift_abort(*cid,"invalid_MNPOUT_peak_readback");
 printf("PROSP_MNPOUT_PEAK %d %d %d %s %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
  *cid,(int)*iteration,*internal_index,name,*source_initial,*proposed,*stored,
  *source_initial,*lower,*upper,*stored_lower,*stored_upper);
 fflush(stdout);
}
