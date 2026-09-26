#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MWgaldust.h"
#define SEV_FATAL 4
char c1err[4096], c2err[4096];
void errmsg(int severity, int ignored, char *fnam, char *one, char *two) {
  fprintf(stderr, "%s: %s; %s\n", fnam, one, two); abort();
}
void concat_callfun_plus_fnam(char *one, char *two, char *out) {
  snprintf(out, 60, "%s", two);
}
double GALextinct(double RV, double AV, double WAVE, int OPT, double *PARLIST, char *callFun) {

/*** 
  
  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    RV   = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE = wavelength in angstroms

    OPT=89 => use original CCM89 :
              Cardelli, Clayton, & Mathis (1989) extinction law.

    OPT=94 => use update from O'Donell

    OPT=-99 => use Fitzpatrick 1999 (PASP 111, 63) as implemented by 
              D.Scolnic with polynomial fit to ratio of F'99/O'94 vs. lambda.
              O'94 is the opt=94 option and F'99 was computed from 
              http://idlastro.gsfc.nasa.gov/ftp/pro/astro/fm_unred.pro
              Deprecated to OPT=-99 from OPT=99 September 25 2024.
              Only reliable for RV=3.1.

    OPT=99 => use Fitzpatrick 1999 (PASP 111, 63) as implemented by S. Thorp.
                This version directly evaluates the cubic spline in inverse
                wavelength, as defined by the fm_unred.pro code. Consistent
                with extinction.py by K. Barbary, and BAYESN F99 implementation.
                Promoted to OPT=99 September 25 2024.

    OPT=203 => use Gordon et al. 2003 (ApJ 594, 279) as implemented by S. Thorp.
                This is the SMC bar dust law. No significant UV bump. Only
                defined for RV=2.74, will abort for all other values. This
                is the refined version from Gordon et al. 2016 (ApJ, 826, 104),
                based on the implementation in Gordon 2024 (JOSS 9, 7023).
                Not recommended for use as a standalone dust law.

    OPT=204 => use Fitzpatrick 2004 (ASP Conf. Ser. 309, 33) as implemented by S. Thorp.
                This uses the same curve as Fitzpatrick 99, but with updated
                behaviour in the IR. Checked against implementation by
                Gordon 2024 (JOSS 9, 7023): github.com/karllark/dust_extinction.

    OPT=208 => use Goobar 2008 (ApJ 686, L103) power law for circumstellar dust.
                This is a two parameter model, controlled by P and A.
                P is read from PARLIST[0];
                A is read from PARLIST[1].
                RV argument is ignored.
                Aborts if the required PARLIST entries are not present or
                within the valid ranges. P=-1.5, A=0.9 gives MW-like circum-
                stellar dust (G08 fit to Draine 2003). P=-2.5, A=0.8 gives
                LMC-like circumstellar dust (G08 fit to Weingartner & Draine 2001).

    OPT=214=> use Maiz Apellaniz et al. 2014 (A&A 564, A63) CCM-like curve.
                Only valid above 0.3 microns. Tested against Gordon 2024 version.

    OPT=216 => use Gordon et al. 2016 (ApJ 826, 104) as implemented by S. Thorp.
                This is a two parameter model, controlled by RVA and FA.
                RVA is read from PARLIST[0];
                FA is read from PARLIST[1].
                RV argument is ignored.
                Aborts if these are not present or within the valid ranges.
                The final curve is a mixture of Fitzpatrick 1999 and
                Gordon et al. 2003 (ApJ 594, 279), where the latter is SMC bar-like
                dust with RV=2.74 and no UV bump. FA=1 reverts to Fitzpatrick 1999 
                with RV=RVA. FA=0 gives Gordon et al. 2003 with RV=2.74. 
                Effective RV = 1/[FA/RV + (1-FA)/2.74].
                Tested against Gordon 2024 implementation.

    OPT=223 => use Gordon et al. 2023 (ApJ 950, 86) as implemented by S. Thorp.
                This is a full UV-OPT-IR extinction law parameterized by RV.
                Defined by a combination of Fitzpatrick & Massa 1990 (ApJS 72, 163)
                in the UV plus various other functions composed together at
                other wavelengths. Tested against Gordon 2024 (JOSS 9, 7023).

    OPT=225 => use Sommovigo et al. 2025 (ApJ 990, 114) as implemented by S. Thorp.
                This is a one-parameter extinction law based on a 4-parameter Pei-like
                functional form, and some scaling relations for the coefficients
                (c1, c2, c3, c4) as a function of AV. RV is ignored as the shape is
                entirely set by AV. See Eq. 2, 7, 8, 9 in the Sommovigo paper. 
                Based on fits to simulations by the Learning the Universe collaboration.

    OPT=226 => use Sommovigo et al. 2026 (arXiv:2606.10027) implemented by S. Thorp.
                This is a 4-parameter functional form derived using symbolic regression.
                B0 is read from PARLIST[0];
                B1 is read from PARLIST[1];
                B2 is read from PARLIST[2];
                B3 is read from PARLIST[3].
                Aborts if B0 or B3 is negative, or if not enough
                parameters are found.
                RV is ignored.

   PARLIST => optional set of double-precision parameters to refine calculations
              Number of PARLIST params and their meaning depend on OPT.
              OPT=208 : PARLIST[0]=P, PARLIST[1]=A;
              OPT=216 : PARLIST[0]=RVA, PARLIST[1]=FA;
              OPT=226 : PARLIST[0]=B0, PARLIST[1]=B1, PARLIST[2]=B2, PARLIST[3]=B3;

  Returns magnitudes of extinction.

 Nov 1, 2006: Add option to use new/old NIR coefficients
              (copied from Jha's MLCS code)

;     c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,    $ ;Original
;                 0.01979, -0.77530,  0.32999 ]               ;coefficients
;     c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,    $ ;from CCM89
;                -0.62251,  5.30260, -2.09002 ]

      c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,    $    ;New coefficients
                 -1.718,   -0.827,    1.647, -0.505 ]        ;from O'Donnell
      c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,    $    ;(1994)
                 11.102,    5.491,  -10.805,  3.347 ]

  Aug 4 2019 RK
   + fix subtle bug by returning XT=0 only if AV=0, and not if AV<1E-9.
     Recall that negative AV are used for warping spectra in kcor.c.
     This bug caused all AV<0 to have same mag as AV=0.

  Sep 19 2024 ST
   + add an exact Fitzpatrick 99 implementation with opt=9999.

  Sep 25 2024 S.Thorp
   + Exact F'99 spline implementation promoted to opt=99
   - Old F'99 based on F'99/O'94 ratio deprecated to opt=-99

  Oct 19 2024 S. Thorp
   + Begun adding more dust laws

  Oct 24 2024 R.Kessler
   + pass PARLIST based on sim-input key PARLIST_MWCOLORLAW: p0,p1,p2,...
     PARLIST is not used yet, but is available for future development.

  Oct 26 2024 S. Thorp
   + use the new PARLIST for Gordon '16 dust law
   + add Goobar '08 circumstellar dust law
   + add Maiz Apellaniz '14

  Feb 26 2025 S. Thorp
   + add Sommovigo '25
   + add 4-parameter Pei '92 curve (Li '08)

  Jun 26 2026 S. Thorp
   + add Sommovigo '26
 ***/

  int i, DO94  ;
  double XT, x, y, a, b, fa, fb, xpow, xx, xx2, xx3 ;
  double y2, y3, y4, y5, y6, y7, y8 ;


  char fnam[60];
  concat_callfun_plus_fnam(callFun, "GALextinct", fnam);  (void)fnam;

  // ------------------- BEGIN --------------

  XT = 0.0 ;

  if ( AV == 0.0  )  {  return XT ; }

  // -----------------------------------------
  // if selecting non-CCM89-like option,
  // bypass everything else and call S. Thorp's functions

  //  printf(" xxx %s: PARLIST = %f %f %f \n", PARLIST[0], PARLIST[1], PARLIST[2] ); fflush(stdout);
  
  if ( OPT == OPT_MWCOLORLAW_FITZ99_EXACT || OPT == OPT_MWCOLORLAW_FITZ04 || OPT == OPT_MWCOLORLAW_GORD03 ) {
    XT = GALextinct_Fitz99_exact(RV, AV, WAVE, OPT, callFun);
    return XT ;
  } else if ( OPT == OPT_MWCOLORLAW_GOOB08 ) {
      double WAVE0 = 5495.0; // reference V-band wavelength
      double P = PARLIST[0]; //extract power law index from PARLIST
      double A = PARLIST[1]; //extract prefactor from PARLIST
      // try to catch missing arguments
      if ( PARLIST[0] == -99.0 || PARLIST[1] == -99.0 ) {
          sprintf(c1err,"Found suspicious inputs: PARLIST[0]=%.1f and PARLIST[1]=%.1f",
                  PARLIST[0], PARLIST[1]);
          sprintf(c2err,"Goobar (2008) requires two values in PARLIST_MWCOLORLAW: P,A.");
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      // check parameter ranges
      if ( P > PMAX_GOOB08 || P < PMIN_GOOB08 ){
          sprintf(c1err,"Read invalid P=%.1f from PARLIST_MWCOLORLAW!", P);
          sprintf(c2err,"Goobar (2008) only recommended for %.1f<=P<=%.1f.",
                  PMIN_GOOB08, PMAX_GOOB08);
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      if ( A > 1.0 || A <= 0.0 ){
          sprintf(c1err,"Read invalid A=%.1f from PARLIST_MWCOLORLAW!", A);
          sprintf(c2err,"Goobar (2008) only valid for 0.0<A<=1.0.");
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      // check wavelength range
      if ( WAVE < WAVEMIN_GOOB08 || WAVE > WAVEMAX_GOOB08 ) {
          sprintf(c1err,"WAVE=%.1f out of range for Goobar (2008)", WAVE);
          sprintf(c2err,"Recommended limits are %.1f<=WAVE<=%.1f.", 
                  WAVEMIN_GOOB08, WAVEMAX_GOOB08);
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      // power law (eq. 3 in G08)
      XT = 1.0 - A + A*pow(WAVE/WAVE0, P);
      return AV*XT;
  } else if ( OPT == OPT_MWCOLORLAW_MAIZ14 ) {
    XT = GALextinct_Maiz14(RV, AV, WAVE, callFun);
    return XT;
  } else if ( OPT == OPT_MWCOLORLAW_GORD16 ) {
      double XTA, XTB;
      double RVB = 2.74; // R,K. -- ensure double cast for this param
      double RVA = PARLIST[0]; // extract RVA from PARLIST
      double FA  = PARLIST[1]; // extract FA from PARLIST
      // sanity check arguments from PARLIST
      // try to catch missing arguments
      if ( PARLIST[0] == -99.0 || PARLIST[1] == -99.0 ) {
          sprintf(c1err,"Found suspicious inputs: PARLIST[0]=%.1f and PARLIST[1]=%.1f",
                  PARLIST[0], PARLIST[1]);
          sprintf(c2err,"Gordon et al. (2016) requires two values in PARLIST_MWCOLORLAW: RVA,FA.");
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      // check parameter ranges
      if ( RVA > RVMAX_FITZ99 || RVA < RVMIN_FITZ99 ){
          sprintf(c1err,"Read invalid RVA=%.1f from PARLIST_MWCOLORLAW!", RVA);
          sprintf(c2err,"Gordon et al. (2016) only valid for %.1f<=RVA<=%.1f.",
                  RVMIN_FITZ99, RVMAX_FITZ99);
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      if ( FA > 1.0 || FA < 0.0 ){
          sprintf(c1err,"Read invalid FA=%.1f from PARLIST_MWCOLORLAW!", FA);
          sprintf(c2err,"Gordon et al. (2016) only valid for 0.0<=FA<=1.0.");
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      XTA = GALextinct_Fitz99_exact(RVA, AV, WAVE, OPT_MWCOLORLAW_FITZ99_EXACT, callFun);
      XTB = GALextinct_Fitz99_exact(RVB, AV, WAVE, OPT_MWCOLORLAW_GORD03,       callFun);

      return FA*XTA + (1-FA)*XTB ;
      
  } else if ( abs(OPT) == OPT_MWCOLORLAW_FITZ19_CUBIC ) {
    XT = GALextinct_Fitz19(RV, AV, WAVE, (OPT>0) ? 1 : 0,  callFun);
    return XT;
  } else if ( OPT == OPT_MWCOLORLAW_GORD23 ) {
    XT = GALextinct_Gord23(RV, AV, WAVE, callFun);
    return XT;
  } else if ( OPT == OPT_MWCOLORLAW_SOMM25 ) {
    XT = GALextinct_Somm25(AV, WAVE, callFun);
    return XT;
  } else if ( OPT == OPT_MWCOLORLAW_SOMM26 ) {
    if ( PARLIST[0] == -99.0 || PARLIST[1] == -99.0 || PARLIST[2] == -99.0 || PARLIST[3] == -99.0 ) {
      sprintf(c1err,"Found suspicious inputs: PARLIST[0]=%.1f, PARLIST[1]=%.1f, PARLIST[2]=%.1f, PARLIST[3]=%.1f.",
              PARLIST[0], PARLIST[1], PARLIST[2], PARLIST[3]);
      sprintf(c2err,"Sommovigo et al. (2026) requires four values in PARLIST_MWCOLORLAW: B0,B1,B2,B3.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    XT = GALextinct_Somm26(AV, WAVE, PARLIST[0], PARLIST[1], PARLIST[2], PARLIST[3], callFun);
    return XT;
  }
  
  // -----------------------------------------
  DO94 = (OPT == OPT_MWCOLORLAW_ODON94 ||
	  OPT == OPT_MWCOLORLAW_FITZ99_APPROX ) ;

  x = 10000./WAVE;    // inverse wavelength in microns
  y = x - 1.82;

  if (x >= 0.3 && x < 1.1) {           // IR
    xpow = pow(x,1.61) ;
    a =  0.574 * xpow ;
    b = -0.527 * xpow ;
  } 
  else if (x >= 1.1 && x < 3.3) {    // Optical/NIR

    y2 = y  * y ;
    y3 = y2 * y ;
    y4 = y2 * y2 ;
    y5 = y3 * y2;
    y6 = y3 * y3;
    y7 = y4 * y3;
    y8 = y4 * y4;

    if ( DO94 ) {
    a = 1. + 0.104*y - 0.609*y2 + 0.701*y3 + 1.137*y4
      - 1.718*y5 - 0.827*y6 + 1.647*y7 -0.505*y8 ;

    b = 1.952*y + 2.908*y2 -3.989*y3 - 7.985*y4
      + 11.102*y5 + 5.491*y6 - 10.805*y7 + 3.347*y8;
    }
    else {
    a = 1. + 0.17699*y - 0.50447*y2 - 0.02427*y3
      + 0.72085*y4 + 0.01979*y5 - 0.77530*y6 + 0.32999*y7 ;

    b = 1.41338*y + 2.28305*y2 + 1.07233*y3 - 5.38434*y4
      - 0.62251*y5 + 5.30260*y6 - 2.09002*y7 ;
    }

  } 
  else if (x >= 3.3 && x < 8.0 ) {    // UV
    if (x >= 5.9) {
      xx  = x - 5.9 ;
      xx2 = xx  * xx ;
      xx3 = xx2 * xx ;

      fa = -0.04473*xx2 - 0.009779*xx3 ;
      fb =  0.21300*xx2 + 0.120700*xx3 ;

    } else {
      fa = fb = 0.0;
    }

    xx = x - 4.67 ; xx2 = (xx*xx);
    a =  1.752 - 0.316*x - 0.104/(xx2 + 0.341) + fa;

    xx = x - 4.62 ; xx2 = (xx*xx);
    b = -3.090 + 1.825*x + 1.206/(xx2 + 0.263) + fb;
  } 
  else if (x >= 8.0 && x <= 10.0) {  // Far-UV
    xx  = x - 8.0  ;
    xx2 = xx  * xx ;
    xx3 = xx2 * xx ; 

    a = -1.073 - 0.628*xx + 0.137*xx2 - 0.070*xx3 ;
    b = 13.670 + 4.257*xx - 0.420*xx2 + 0.374*xx3 ;
  } else {
    a = b = 0.0;
  }

  XT = AV*(a + b/RV);

  // Sep 18 2013 RK/DS - Check option for Fitzptrack 99

#define NPOLY_FITZ99 11 //Dillon and Dan upped to 10, Oct 9 2021
  if ( OPT == OPT_MWCOLORLAW_FITZ99_APPROX ) {  

    double XTcor, wpow[NPOLY_FITZ99] ;    
    double F99_over_O94[NPOLY_FITZ99] = {  // Dillon and Dan, Oct 9 2021
      8.55929205e-02,  1.91547833e+00, -1.65101945e+00,  7.50611119e-01,
      -2.00041118e-01,  3.30155576e-02, -3.46344458e-03,  2.30741420e-04,
      -9.43018242e-06,  2.14917977e-07, -2.08276810e-09
    };

    if ( WAVE > WAVEMAX_FITZ99  ) {
      sprintf(c1err,"Invalid WAVE=%.1f A for Fitzpatrick 99 color law.",
	      WAVE );
      sprintf(c2err,"Avoid NIR (>%.1f), or update Fitz99 in NIR",
	      WAVEMAX_FITZ99 );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // compute powers of wavelength without using slow 'pow' function
    wpow[0]  = 1.0 ;
    wpow[1]  = WAVE/1000. ;
    wpow[2]  = wpow[1] * wpow[1] ;
    wpow[3]  = wpow[1] * wpow[2] ;
    wpow[4]  = wpow[2] * wpow[2] ;
    wpow[5]  = wpow[3] * wpow[2] ;
    wpow[6]  = wpow[3] * wpow[3] ;
    wpow[7]  = wpow[4] * wpow[3] ;
    wpow[8]  = wpow[4] * wpow[4] ;
    wpow[9]  = wpow[5] * wpow[4] ;
    wpow[10]  = wpow[5] * wpow[5] ;

    XTcor = 0.0 ;
    for(i=0; i < NPOLY_FITZ99; i++ ) 
      {  XTcor += (wpow[i] * F99_over_O94[i]) ; }
    
    XT *= XTcor ;
  }

  return XT ;

}  // end of GALextinct


// ============= EXACT F99 EXTINCTION LAW ==============
double GALextinct_Fitz99_exact(double RV, double AV, double WAVE, int OPT, char *callFun) {
/*** 
  Created by S. Thorp, Sep 19 2024

  Default Fitzpatrick (1999) implementation since Sep 25 2024

  Also used to compute Fitzpatrick (2004), Gordon et al. (2003),
  and Gordon et al. (2016) laws.

  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    RV   = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE = wavelength in angstroms
    OPT  = Option from (99, 203, 204, 216)
Returns :
    XT = magnitudes of extinction
***/

  char fnam[60];
  concat_callfun_plus_fnam(callFun, "GALextinct_Fitz99_exact", fnam);  (void)fnam;

  //Check RV=2.74 for Gordon et al. (2003)
  if ( OPT == OPT_MWCOLORLAW_GORD03 && RV != 2.74 ) {
    sprintf(c1err,"Requested OPT=%d and RV=%.2f", OPT, RV);
    sprintf(c2err,"Gordon et al. 2003 only valid for RV=2.74");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  //Check wavelengths in valid range
  if ( WAVE < WAVEMIN_FITZ99_EXACT || WAVE > WAVEMAX_FITZ99_EXACT ) {
    sprintf(c1err,"Requested WAVE=%.3f Angstroms", WAVE);
    sprintf(c2err,"F99-like curves only valid in [%.1f, %.1f]A",
	    WAVEMIN_FITZ99_EXACT, WAVEMAX_FITZ99_EXACT);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  //number of knots
  int Nk = 0 ;
  // constants
  double x02=0.0, gamma2=0.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0, c5=0.0 ;
  // target wavenumber in inverse microns
  double x = 10000.0/WAVE;
  // spline result
  double y;

  // constants
  c2 = -0.824 + 4.717/RV;
  c5 = 5.90;
  if ( OPT == OPT_MWCOLORLAW_FITZ99_EXACT ) {
    x02 = 21.123216; // 4.596*4.596
    gamma2 = 0.9801; // 0.99*0.99
    c1 = 2.03 - 3.007*c2;
    c3 = 3.23;
    c4 = 0.41;
    Nk = 9;
  } else if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
    x02 = 21.086464; // 4.592*4.592
    gamma2 = 0.850084; // 0.922*0.922
    c1 = 2.18 - 2.91*c2;
    c3 = 2.991;
    c4 = 0.319;
    Nk = 10;
  } else if ( OPT == OPT_MWCOLORLAW_GORD03 ) {
    x02 = 21.16; // 4.6*4.6
    gamma2 = 1.0;
    c1 = -4.959;
    c2 = 2.264;
    c3 = 0.389;
    c4 = 0.461;
    Nk = 11;
  } else {
    sprintf(c1err,"Requested OPT=%d", OPT);
    sprintf(c2err,"Only 99, 203, 204 are implemented!");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if (WAVE <= 2700.0) { //FM90 curve in UV
    y = GALextinct_FM90(x, c1, c2, c3, c4, c5, x02, gamma2);
    return AV * (1.0 + y/RV); 
  } else { //spline for optical/IR
    // powers of RV
    double RV2, RV3, RV4;
    
    // spline knot locations in inverse microns
    double xF[Nk];
    xF[0] = 0.0; // always put an anchor at 1/lambda = 0
    if ( OPT == OPT_MWCOLORLAW_GORD03 ) {
      xF[1] = 1.0/2.198;
      xF[2] = 1.0/1.65;
      xF[3] = 1.0/1.25;
      xF[4] = 1.0/0.81;
      xF[5] = 1.0/0.65;
      xF[6] = 1.0/0.55;
      xF[7] = 1.0/0.44;
      xF[8] = 1.0/0.37;
    } else {
      if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
	// Use FM07 knots for Fitzpatrick (2004) curve
	xF[1] = 0.5;
	xF[2] = 0.75;
	xF[3] = 1.0;
      } else {
	xF[1] = 1.0/2.65;
	xF[2] = 1.0/1.22;
      }
      xF[Nk-6] = 1.0/0.60;
      xF[Nk-5] = 1.0/0.547;
      xF[Nk-4] = 1.0/0.467; 
      xF[Nk-3] = 1.0/0.411;
    }
    // always anchor in the UV
    xF[Nk-2] = 1.0/0.270;
    xF[Nk-1] = 1.0/0.260;
    // spline knot values
    double yF[Nk];

    // RV-dependent spline knot values
    // polynomial coeffs match FM_UNRED.pro and extinction.py
    // NOTE: the optical coefficients differ from Gordon 24 implementation
    double yFNIR;
    yF[0] = -RV;
    if ( OPT == OPT_MWCOLORLAW_GORD03 ) {
      // knot values have 1 subtracted and are multiplied by RV
      yF[1] = -2.4386; //0.11*RV-RV
      yF[2] = -2.27694; //0.169*RV-RV
      yF[3] = -2.055; //0.25*RV-RV
      yF[4] = -1.18642; //0.567*RV-RV
      yF[5] = -0.54526; //0.801*RV-RV
      yF[6] = 0.0;
      yF[7] = 1.02476; //1.374*RV-RV 
      yF[8] = 1.84128; //1.672*RV-RV
    } else {
      // powers of RV
      RV2 = RV*RV;
      RV3 = RV2*RV;
      RV4 = RV2*RV2;
      if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
	yFNIR = (0.63*RV -0.84);
	yF[1] = yFNIR*pow(xF[1], 1.84) - RV;
	yF[2] = yFNIR*pow(xF[2], 1.84) - RV;
	yF[3] = yFNIR*pow(xF[3], 1.84) - RV;
      }
      else {
	yF[1] = -0.914616129*RV; // 0.26469*(RV/3.1) - RV
	yF[2] = -0.7325*RV; // 0.82925*(RV/3.1) - RV
      }
      yF[Nk-6] = -0.422809 + 0.00270*RV +  2.13572e-04*RV2;
      yF[Nk-5] = -5.13540e-02 + 0.00216*RV - 7.35778e-05*RV2;
      yF[Nk-4] =  7.00127e-01 + 0.00184*RV - 3.32598e-05*RV2;
      yF[Nk-3] =  1.19456 + 0.01707*RV - 5.46959e-03*RV2 +  
	7.97809e-04*RV3 - 4.45636e-05*RV4;
    }
    // UV knots using FM90
    yF[Nk-2] = GALextinct_FM90(xF[Nk-2], c1, c2, c3, c4, c5, x02, gamma2);
    yF[Nk-1] = GALextinct_FM90(xF[Nk-1], c1, c2, c3, c4, c5, x02, gamma2);
    
    y = GALextinct_FM_spline(x, Nk, xF, yF, 0);
    
    return AV*(1.0 + y/RV);
  }

} // end of GALextinct_Fitz99_exact

// ============= MAIZ APELLANIZ ET AL. 2014 EXTINCTION LAW ==============
double GALextinct_Maiz14(double RV, double AV, double WAVE, char *callFun) {
/*** 
  Created by S. Thorp, Oct 26 2024

  Input : 
    AV    = V band (defined to be at 5495 Angstroms) extinction
    RV    = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE  = wavelength in angstroms
  Returns :
    XT = magnitudes of extinction
***/

    char fnam[60] ;
    concat_callfun_plus_fnam(callFun, "GALextinct_Maiz14", fnam); (void)fnam;

    // Abort if out of bounds
    if ( WAVE > WAVEMAX_MAIZ14 || WAVE < WAVEMIN_MAIZ14 ) {
      sprintf(c1err,"Requested WAVE=%.3f Angstroms", WAVE);
      sprintf(c2err,"Maiz Apellaniz et al. 2014 only valid from %.0f-%.0f Angstroms",
              WAVEMIN_MAIZ14, WAVEMAX_MAIZ14);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // target wavelength in inverse microns
    double x = 10000.0/WAVE;

    // terms we'll compute
    double a, b; //a and b curves at x

    // evaluate the a and b curves
    if (x < 1.0) { // just do the IR power law above 1 micron
        a =  0.574 * pow(x, 1.61);
        b = -0.527 * pow(x, 1.61);
    } else { // do the spline
        // all knots are independent of RV
        // coefficients extracted from Gordon's SciPy spline implementation
        double xk[11] = { 1.0, 1.15, 1.81984, 2.1, 2.27015, 2.7, 3.5,
            3.9, 4.0, 4.1, 4.2 }; // knot positions
        double a3[10] = { -3.09348541,  2.28902153e-1,  5.41605406e-1,
            -6.37404842e-1,  3.52950213e-1, -5.91231605e-2,
            -5.56727269,  48.1384135, -11.6556097, -12.6892172 }; //x^3
        double a2[10] = { 5.57088021e-1, -8.34980412e-1, -3.74996957e-1,
            8.02115549e-2, -2.45151747e-1,  2.09995201e-1,  6.80996157e-2,
            -6.61262761, 7.82889643,  4.33221353 };  //x^2
        double a1[10] = { 9.24140000e-1,  8.82456141e-1,  7.19649009e-2,
            -1.06221772e-2, -3.86867508e-2, -5.37987921e-2, 1.68677061e-1,
            -2.44913414, -2.32750725, -1.11139626 }; //x^1
        double a0[10] = { 5.74000000e-1,  7.14714967e-1,  9.99971669e-1,
            1.00260970,  9.99984676e-1,  9.66090893e-1, 1.02717773,  
            7.49239041e-1,  4.86337764e-1, 3.20220393e-1 }; //x^0
        double b3[10] = { 6.11543973, -4.71924979e-1, -3.75700076,
            3.30710701, -6.80610047e-1,  4.81511488e-1, 17.8352808,
            -124.325934,  12.0120271, 48.1516935 }; //x^3
        double b2[10] = { -2.49479124e-1,  2.50246875,  1.55412607,
            -1.60355793,  8.45548471e-2, -7.93125839e-1, 3.62501733e-1,
            21.7648387, -15.5329415, -11.9293334 }; //x^2
        double b1[10] = { -8.48470000e-1, -5.10521556e-1,  2.20674792,
            2.19289909,  1.93444072,  1.62986148, 1.28536219,
            10.1362984,  10.7594881, 8.01326059 }; //x^1
        double b0[10] = { -5.27000000e-1, -6.39244171e-1, -2.26082358e-4,
            6.57384043e-1,  1.00037205,  1.79345802, 2.83628055,  
            4.54988367,  5.65683596, 6.58946738 };
       

        // find index in knot list
        int q = 0; // qmin = 0; qmax = 9
        while (q < 10) {
            if (x < xk[q+1]) { 
                break; 
            } else {
                q++;
            }
        }

        //powers of x
        double x1 = x - xk[q];
        double x2 = x1*x1;
        double x3 = x2*x1;
        
        // interpolate
        a = a3[q]*x3 + a2[q]*x2 + a1[q]*x1 + a0[q];
        b = b3[q]*x3 + b2[q]*x2 + b1[q]*x1 + b0[q];

    }
    return AV * (a + b/RV);

} // end of GALextinct_Maiz14

// ============= FITZPATRICK ET AL. 2019 EXTINCTION LAW ==============
double GALextinct_Fitz19(double RV, double AV, double WAVE, int CUBIC, char *callFun) {
/*** 
  Created by S. Thorp, Oct 20 2024

  Input : 
    AV    = V band (defined to be at 5495 Angstroms) extinction
    RV    = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE  = wavelength in angstroms
    CUBIC = if 1, uses cubic interpolation; else linear interpolation
  Returns :
    XT = magnitudes of extinction
***/

    char fnam[60] ;
    concat_callfun_plus_fnam(callFun, "GALextinct_Fitz19", fnam);  (void)fnam;

    // Abort if out of bounds
    if ( WAVE > WAVEMAX_FITZ19 || WAVE < WAVEMIN_FITZ19 ) {
      sprintf(c1err,"Requested WAVE=%.3f Angstroms", WAVE);
      sprintf(c2err,"Fitzpatrick et al. 2019 only valid from %.0f-%.0f Angstroms",
              WAVEMIN_FITZ19, WAVEMAX_FITZ19);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // target wavelength in inverse microns
    double x = 10000.0/WAVE;
    // curve at WAVE
    double y;

    // number of tabulated values
    int Nk = 102;

    //tabulated values of E(x-V)/E(B-V)
    //based on dust_extinction / Table 3 in Fitzpatrick+19
    double xk[102] = { 0.000, 0.455, 0.606, 0.800, 1.000, 1.100, 1.200, 1.250,
        1.300, 1.350, 1.400, 1.450, 1.500, 1.550, 1.600, 1.650, 1.700, 1.750,
        1.800, 1.818, 1.850, 1.900, 1.950, 2.000, 2.050, 2.100, 2.150, 2.200,
        2.250, 2.273, 2.300, 2.350, 2.400, 2.450, 2.500, 2.550, 2.600, 2.650,
        2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000, 3.100, 3.200, 3.300,
        3.400, 3.500, 3.600, 3.700, 3.800, 3.900, 4.000, 4.100, 4.200, 4.300,
        4.400, 4.500, 4.600, 4.700, 4.800, 4.900, 5.000, 5.100, 5.200, 5.300,
        5.400, 5.500, 5.600, 5.700, 5.800, 5.900, 6.000, 6.100, 6.200, 6.300,
        6.400, 6.500, 6.600, 6.700, 6.800, 6.900, 7.000, 7.100, 7.200, 7.300,
        7.400, 7.500, 7.600, 7.700, 7.800, 7.900, 8.000, 8.100, 8.200, 8.300,
        8.400, 8.500, 8.600, 8.700 }; //knot positions in inverse microns
    double k302k[102] = { -3.020, -2.747, -2.528, -2.222, -1.757, -1.567, -1.300,
        -1.216, -1.070, -0.973, -0.868, -0.750, -0.629, -0.509, -0.407, -0.320,
        -0.221, -0.133, -0.048, 0.000, 0.071, 0.188, 0.319, 0.438, 0.575, 0.665,
        0.744, 0.838, 0.951, 1.000, 1.044, 1.113, 1.181, 1.269, 1.346, 1.405,
        1.476, 1.558, 1.632, 1.723, 1.791, 1.869, 1.948, 2.009, 2.090, 2.253,
        2.408, 2.565, 2.746, 2.933, 3.124, 3.328, 3.550, 3.815, 4.139, 4.534,
        5.012, 5.560, 6.118, 6.565, 6.767, 6.681, 6.394, 6.038, 5.704, 5.432,
        5.226, 5.078, 4.978, 4.913, 4.877, 4.862, 4.864, 4.879, 4.904, 4.938,
        4.982, 5.038, 5.105, 5.181, 5.266, 5.359, 5.460, 5.569, 5.684, 5.805,
        5.933, 6.067, 6.207, 6.352, 6.502, 6.657, 6.817, 6.981, 7.150, 7.323,
        7.500, 7.681, 7.866, 8.054, 8.246, 8.441 }; //k function [R(5500)=3.02]
    double sk[102] = { -1.000, -0.842, -0.728, -0.531, -0.360, -0.284, -0.223,
        -0.198, -0.173, -0.150, -0.130, -0.110, -0.096, -0.081, -0.063, -0.048,
        -0.032, -0.017, -0.005, 0.000, 0.007, 0.013, 0.012, 0.010, 0.004, 0.003,
        0.000, 0.002, 0.001, 0.000, -0.000, 0.001, 0.001, -0.002, 0.000, -0.002,
        -0.002, -0.006, -0.009, -0.011, -0.017, -0.025, -0.029, -0.037, -0.043,
        -0.064, -0.092, -0.122, -0.161, -0.201, -0.249, -0.303, -0.366, -0.437,
        -0.517, -0.603, -0.692, -0.774, -0.843, -0.888, -0.908, -0.903, -0.880,
        -0.849, -0.816, -0.785, -0.760, -0.741, -0.729, -0.722, -0.722, -0.726,
        -0.734, -0.745, -0.760, -0.778, -0.798, -0.820, -0.845, -0.870, -0.898,
        -0.926, -0.956, -0.988, -1.020, -1.053, -1.087, -1.122, -1.158, -1.195,
        -1.232, -1.270, -1.309, -1.349, -1.389, -1.429, -1.471, -1.513, -1.555,
        -1.598, -1.641, -1.685 }; //s function
                                  
    double kRVk[102];
    for (int i=0; i<Nk; i++) { kRVk[i] = k302k[i] + sk[i]*(RV-3.10)*0.99; }

    y = GALextinct_FM_spline(x, Nk, xk, kRVk, CUBIC ? 0 : 1);

    return AV*(1.0 + y/RV);

} //end of GALextinct_Fitz19

// ============= GORDON ET AL. 2023 EXTINCTION LAW ==============
double GALextinct_Gord23(double RV, double AV, double WAVE, char *callFun) {
/*** 
  Created by S. Thorp, Oct 20 2024

  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    RV   = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE = wavelength in angstroms
  Returns :
    XT = magnitudes of extinction
***/

    char fnam[60] ;
    concat_callfun_plus_fnam(callFun, "GALextinct_Gord23", fnam);  (void)fnam;

    // target wavelength in inverse microns
    double x = 10000.0/WAVE;
    
    double x2, x3, x4; // powers of x
    x2 = x*x;
    x3 = x2*x;
    x4 = x2*x2;

    // variables for a and b part of curve
    // w = weighting function in overlap regions
    double a, b, w, f;
    a = b = w = 0.0;

    // terms for the optical part
    double x01, x02, x03, FW1, FW2; //Drude params
    double FX1, FX2, FX3, XX1, XX2, XX3, D1, D2, D3; //derived terms

    // constants for the N-MIR part
    double scale=0.38526, alpha=1.68467, alpha2=0.78791, swave=4.30578, swidth=4.78338,
        sil1_amp=0.06652, sil1_center=9.8434, sil1_fwhm=2.21205, sil1_asym=-0.24703,
        sil2_amp=0.0267, sil2_center=19.58294, sil2_fwhm=17.0, sil2_asym=-0.27;
    double mwave = WAVE / 10000.0; //wavelength in microns
    double fweight, pweight, ratio; //power law transition
    double sil1_gamma, sil2_gamma, sil1_gx2, sil2_gx2, sil1_xx, sil2_xx; //Si drude params

    // Abort if out of bounds
    if ( WAVE > WAVEMAX_GORD23 || WAVE < WAVEMIN_GORD23 ) {
      sprintf(c1err,"Requested WAVE=%.3f Angstroms; X=%.3f inv. microns", WAVE, x);
      sprintf(c2err,"Gordon et al. 2023 only valid from %.0f-%.0f Angstroms",
              WAVEMIN_GORD23, WAVEMAX_GORD23);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // UV / OPT / NIR
    if ( 1.0/0.33 <= x && x <= 1.0/0.09 ) { //UV (including UV-OPT overlap)
        // weighting function
        if ( x > 1.0/0.30 ) {
            w = 1.0;
        } else {
            // Gordon "smoothstep" function
            f = (mwave - 0.3)/0.03;
            w = 3.0 - 2.0*f;
            w = 1.0 - w*f*f;
        }
        // a component: 21.16 = 4.6*4.6; 0.9801 = 0.99*0.99
        a += w*GALextinct_FM90(x, 0.81297, 0.2775, 1.06295, 0.11303, 5.90, 21.16, 0.9801);
        // b component 
        b += w*GALextinct_FM90(x, -2.97868, 1.89808, 3.10334, 0.65484, 5.90, 21.16, 0.9801);
    }
    if ( 1.0/1.1 <= x && x < 1.0/0.3 ) { //OPT (including both overlaps)
        // weighting function
        if ( 1.0/0.9 < x && x < 1.0/0.33 ) { //internal
            w = 1.0;
        } else if ( x >= 1.0/0.33 ) { //overlap with UV
            // Gordon "smoothstep" function
            f = (mwave - 0.3)/0.03;
            w = 3.0 - 2.0*f;
            w = w*f*f;
        } else if ( x <= 1.0/0.9 ) { //overlap with IR
            // Gordon "smoothstep" function
            f = (mwave - 0.9)/0.2;
            w = 3.0 - 2.0*f;
            w = 1.0 - w*f*f;
        }

        // polynomial terms
        a += w*(-0.35848 + 0.7122*x + 0.08746*x2 - 0.05403*x3 + 0.00674*x4);
        b += w*(0.12354 - 2.68335*x + 2.01901*x2 - 0.39299*x3 + 0.03355*x4);

        //the Drude abides
        // shared terms
        x01 = 2.288;
        x02 = 2.054;
        x03 = 1.587;
        FW1 = 0.243; //FW3 = FW1
        FW2 = 0.179;
        FX1 = (FW1*FW1)/(x01*x01);
        FX2 = (FW2*FW2)/(x02*x02);
        FX3 = (FW1*FW1)/(x03*x03);
        XX1 = (x/x01 - x01/x);
        XX2 = (x/x02 - x02/x);
        XX3 = (x/x03 - x03/x);
        D1 = FX1 / (XX1*XX1 + FX1);
        D2 = FX2 / (XX2*XX2 + FX2);
        D3 = FX3 / (XX3*XX3 + FX3);
        // add contributions to a and b curves
        a += w*(0.03893*D1 + 0.02965*D2 + 0.01747*D3);
        b += w*(0.18453*D1 + 0.19728*D2 + 0.1713*D3);

    }
    if ( 1.0/35.0 <= x && x < 1.0/0.9 ) { //IR (including OPT-IR overlap)
        // weighting function
        if ( x < 1.0/1.1 ) {
            w = 1.0;
        } else {
            // Gordon "smoothstep" function
            f = (mwave - 0.9)/0.2;
            w = 3.0 - 2.0*f;
            w = w*f*f;
        }
        // a curve Gordon21 double power law
        // Gordon smoothstep
        fweight = (mwave - (swave - 0.5*swidth))/swidth;
        if (fweight < 0) {
            pweight = 0.0;
        } else if (fweight > 1) {
            pweight = 1.0;
        } else {
            pweight = (3.0 - 2.0*fweight)*fweight*fweight;
        }
        // ratio
        ratio = pow(swave, -alpha)/pow(swave, -alpha2);
        // power law 1
        a += w * scale * (1.0 - pweight) * pow(mwave, -alpha);
        // power law 2
        a += w * scale * ratio * pweight * pow(mwave, -alpha2);
        // silicate features
        sil1_gamma = 2.0 * sil1_fwhm / (1.0 + exp(sil1_asym*(mwave - sil1_center)));
        sil2_gamma = 2.0 * sil2_fwhm / (1.0 + exp(sil2_asym*(mwave - sil2_center)));
        sil1_gx2 = sil1_gamma*sil1_gamma/(sil1_center*sil1_center);
        sil2_gx2 = sil2_gamma*sil2_gamma/(sil2_center*sil2_center);
        sil1_xx = (mwave/sil1_center) - (sil1_center/mwave);
        sil2_xx = (mwave/sil2_center) - (sil2_center/mwave);
        a += w * sil1_amp * sil1_gx2 / (sil1_xx*sil1_xx + sil1_gx2);
        a += w * sil2_amp * sil2_gx2 / (sil2_xx*sil2_xx + sil2_gx2);

        // b curve power law
        b+= -1.01251 * w * pow(x, 1.06099);
    } 

    return AV * ( a + b*((1.0/RV) - (1.0/3.1)) );

} // end of GALextinct_Gord23

// ============= SOMMOVIGO ET AL. 2025 =======================
double GALextinct_Somm25(double AV, double WAVE, char *callFun) {
/*** 
  Created by S. Thorp, Feb 26 2025

  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    WAVE = wavelength in angstroms
  Returns :
    XT = magnitudes of extinction
***/

    char fnam[60] ;
    concat_callfun_plus_fnam(callFun, "GALextinct_Somm25", fnam); // return fnam
    
    // target wavelength in inverse microns
    double x = 10000.0/WAVE;
    
    // Abort if out of bounds
    if ( WAVE > WAVEMAX_SOMM25 || WAVE < WAVEMIN_SOMM25 ) {
      sprintf(c1err,"Requested WAVE=%.3f Angstroms; X=%.3f inv. microns", WAVE, x);
      sprintf(c2err,"Sommovigo et al. 2025 only valid from %.0f-%.0f Angstroms",
              WAVEMIN_SOMM25, WAVEMAX_SOMM25);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // coefficients as a function of AV
    double c1, c2, c3, c4;
    double logc1, logc4, logAV;
    logAV = log10(AV);
    logc1 = -0.37*logAV + 0.75; // Eq. 7
       c1 = pow(10.0, logc1);
       c2 =  1.88;              // median for TNG galaxies
       c3 =  1.21*logc1 - 1.33; // Eq. 8
    logc4 = -0.59*logAV - 1.42; // Eq. 9
       c4 = pow(10.0, logc4);
    
       return AV*GALextinct_Pei4(x, c1, c2, c3, c4);

} // end of GALextinct_Somm25

// ============= SOMMOVIGO ET AL. 2026 =======================
double GALextinct_Somm26(double AV, double WAVE, 
        double B0, double B1, double B2, double B3, char *callFun) {
/*** 
  Created by S. Thorp, Jun 26 2026

  Four-parameter attenuation curve based on symbolic regression
  analysis by Learning the Universe team.

  Input : 
    AV   = V band (defined to be at 5542 Angstroms) extinction
    WAVE = wavelength in angstroms
    B0   ~ bump strength
    B1   ~ FUV slope
    B2   ~ UV-optical curvature
    B3   ~ UV-optical slope
  Returns :
    XT = magnitudes of extinction
***/
    char fnam[60] ;
    concat_callfun_plus_fnam(callFun, "GALextinct_Somm26", fnam) ;

    // Abort if out of bounds
    if ( WAVE > WAVEMAX_SOMM26 || WAVE < WAVEMIN_SOMM26 ) {
      sprintf(c1err,"Requested WAVE=%.3f", WAVE);
      sprintf(c2err,"Sommovigo et al. 2026 only valid from %.0f-%.0f Angstroms",
              WAVEMIN_SOMM26, WAVEMAX_SOMM26);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( B0 < BMIN_SOMM26 || B3 < BMIN_SOMM26 ){
          sprintf(c1err,"Read invalid B0=%.1f or B3=%.1f from PARLIST_MWCOLORLAW!", B0, B3);
          sprintf(c2err,"Sommovigo et al. (2026) only recommended for positive B0 and B3.");
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    double lV = 5542.0;
    double c0 = 0.4002;
    double c1 = 285.6 ;
    double c2 = 0.2092;
    double c3 = 9.223 ;
    double c4 = 1.016 ;
    double k          ;
    double x  = WAVE/lV;

    k =  B0 * (exp(-c1*(x - c0)*(x - c0)) - exp(-c1*(1.0 - c0)*(1.0 - c0))); // bump term
    k += (B1 + B2*(x - c2)) * (x - c2) * (exp(-c3*x) - exp(-c3));
    k += exp(B3 * (tanh(c4) - tanh(c4*x)));

    return AV*k;
} // end of GALextinct_Somm26
 
// ============= FITZPATRICK & MASSA 1990 ====================
double GALextinct_FM90(double x, double c1, double c2, double c3, double c4, 
        double c5, double x02, double g2) {
  /*
  Created by S. Thorp, Oct 20 2024

  Input : 
    x   = wavenumber (inverse microns)
    c1  = y-intercept of linear component
    c2  = slope of linear component
    c3  = bump amplitude
    c4  = FUV rise amplitude
    c5  = FUV transition point
    x02 = x0*x0, centroid of bump squared
    g2  = gamma*gamma, width of bump squared
Returns :
    E(x-V)/E(B-V)
  */

  double x2, y, y2, b, k;
  char fnam[] = "GALextinct_FM90" ;  (void)fnam;

  x2 = x*x;
  b = x2 / ((x2-x02)*(x2-x02) + x2*g2);
  k = c1 + c2*x + c3*b;
  if (x >= c5) {
    y = x - c5;
    y2 = y * y;
    k += c4 * (0.5392*y2 + 0.05644*y2*y);
  }
  return k;

} // end of GALextinct_FM90

// ============= PEI 1992 / LI ET AL. 2008 ===================
double GALextinct_Pei4(double x, double c1, double c2, double c3, double c4) {
  /*
  Created by S. Thorp, Feb 26 2025

  Four-parameter version of the Pei 1992 (ApJ 395, 130) extinction curve.
  From Li et al. 2008 (ApJ 685, 1046) and Sommovigo et al. 2025 (arXiv:2502.13240).

  Input : 
    x   = wavenumber (inverse microns)
    c1  = UV rise
    c2  = slope
    c3  = FUV shape
    c4  = bump strength
  Returns :
    Ax/AV
  */

    char fnam[] = "GALextinct_Pei4" ;  (void)fnam;

    double x08, x046, x2175;
    double y08, y046, y2175;
    double k, b; 
    x08 = x*0.08;
    x046 = x*0.046;
    x2175 = x*0.2175;
    y08 = pow(x08, c2);
    y046 = x046*x046;
    y2175 = x2175*x2175;
    b = pow(0.145, c2);

    k = c1 / (y08 + 1.0/y08 + c3);
    k += 233.0*(1.0 - c4/4.60 - c1/(b + 1.0/b + c3)) / (y046 + 1.0/y046 + 90.0);
    k += c4 / (y2175 + 1.0/y2175 - 1.95);
    return k;

} // end of GALextinct_Pei4

// ============= FM_UNRED SPLINE ====================
double GALextinct_FM_spline(double x, int Nk, double *xk, double *yk, int lin) {
  /*
  Created by S. Thorp, Oct 22 2024

  Natural cubic spline (a la FM_UNRED).

  Option to return after linear term computed.

  Input :
    x   =  Target to evaluate curve at (inverse microns)
    Nk  =  Number of spline knots (including edges)
    xk  =  Locations of spline knots (inverse microns)
    yk  =  Value of curve at knot locations
    lin =  If 1, quit early and return a linear interpolation.
Returns :
    y   =  Value of curve at x.
  */

    char fnam[] = "GALextinct_FM_spline" ; (void)fnam;

    // abort on x out of knot range
    if (x < xk[0] || x > xk[Nk-1]) {
      sprintf(c1err,"Spline interpolation out of bounds!");
      sprintf(c2err,"Requested %.3f. Limits are [%.3f, %.3f].", 
              x, xk[0], xk[Nk-1]);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    int j, q; //indexes
    double A, B, C, D; //coefficients
    double deltax, deltax2; // gap between bounding knots
    double y; //result

    // find index in knot list
    q = 0; // qmin = 0; qmax = Nk-2
    while (q < Nk-1) {
        if (x < xk[q+1]) { 
            break; 
        } else {
            q++;
        }
    }

    // evaluate the spline (linear part & coeffs)
    deltax = xk[q+1] - xk[q];
    deltax2 = deltax * deltax;
    A = (xk[q+1] - x) / deltax;
    B = 1.0 - A;
    y = A*yk[q] + B*yk[q+1];
    if ( lin == 1 ) { return y; } //stop at linear part
    // cubic part below
    C = (A*A*A - A) * deltax2 / 6.0;
    D = (B*B*B - B) * deltax2 / 6.0;

    // compute 2nd derivatives
    // tridiagonal solve using Thomas algorithm
    double d2yq = 0.0;
    double d2yq1 = 0.0;
    double Kb[Nk-2]; //main diagonal in tridiagonal system
    double Kc[Nk-3]; //off-diagonal in tridiagonal system
    double Vd[Nk-2]; //right hand side
    double wj; //scratch variable
    // fill vectors
    for (j=0; j<Nk-2; j++) {
        Kb[j] = (xk[j+2] - xk[j])/3.0;
        if (j<Nk-3) { Kc[j] = (xk[j+2] - xk[j+1])/6.0; }
        Vd[j] = (yk[j+2] - yk[j+1])/(xk[j+2] - xk[j+1]) - (yk[j+1] - yk[j])/(xk[j+1] - xk[j]);
    }
    // forward substitution
    for (j=1; j<Nk-2; j++) {
        wj = Kc[j-1]/Kb[j-1]; //w factor
        Kb[j] -= wj*Kc[j-1]; //update diagonal
        Vd[j] -= wj*Vd[j-1]; //update rhs
    } // forward substitution complete
    // back substitution (stop at q)
    d2yq = Vd[Nk-3]/Kb[Nk-3]; //final element of solution
    // if q=Nk-2, terminate straight away
    // otherwise enter, the loop and work back to q
    for (j=Nk-4; j>q-2; j--) { //loop backwards
        d2yq1 = d2yq; //shift previous element
        if (j<0) { //we've gone far enough
            d2yq = 0;
            break;
        } else { //find next element
            d2yq = (Vd[j] - Kc[j]*d2yq1)/Kb[j];
        }
    }
    y += C*d2yq + D*d2yq1;
    return y;
} //end GALextinct_FM_spline


double GALextinct_historical(double RV, double AV, double WAVE, int OPT) {

/*** 
  
  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    RV   = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE = wavelength in angstroms

    OPT=89 => use original CCM89 :
              Cardelli, Clayton, & Mathis (1989) extinction law.

    OPT=94 => use update from O'Donell

    OPT=99 => use Fitzpatrick 99 (PASP 111, 63) as implemented by 
              D.Scolnic with polynomial fit to ratio of F'99/O'94 vs. lambda.
              O'94 is the opt=94 option and F'99 was computed from 
              http://idlastro.gsfc.nasa.gov/ftp/pro/astro/fm_unred.pro

  Returns magnitudes of extinction.

 Nov 1, 2006: Add option to use new/old NIR coefficients
              (copied from Jha's MLCS code)

;     c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,    $ ;Original
;                 0.01979, -0.77530,  0.32999 ]               ;coefficients
;     c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,    $ ;from CCM89
;                -0.62251,  5.30260, -2.09002 ]

      c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,    $    ;New coefficients
                 -1.718,   -0.827,    1.647, -0.505 ]        ;from O'Donnell
      c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,    $    ;(1994)
                 11.102,    5.491,  -10.805,  3.347 ]

  Sep 18 2013 RK 
    - add opt=99 option to use Fitzpatrick 99 update.
    - reduce/remove pow calls to save CPU
    - rename CCMextinct -> GALextinct

  Aug 4 2019 RK
   + fix subtle bug by returning XT=0 only if AV=0, and not if AV<1E-9.
     Recall that negative AV are used for warping spectra in kcor.c.
     This bug caused all AV<0 to have same mag as AV=0.

 ***/

  int i, DO94  ;
  double XT, x, y, a, b, fa, fb, xpow, xx, xx2, xx3 ;
  double y2, y3, y4, y5, y6, y7, y8 ;
  char fnam[] = "GALextinct" ;

  // ------------------- BEGIN --------------

  XT = 0.0 ;

  if ( AV == 0.0  )  {  return XT ; }

  // -----------------------------------------
  DO94 = (OPT == 94 || OPT == 99 ) ;

  x = 10000./WAVE;    // inverse wavelength in microns
  y = x - 1.82;

  if (x >= 0.3 && x < 1.1) {           // IR
    xpow = pow(x,1.61) ;
    a =  0.574 * xpow ;
    b = -0.527 * xpow ;
  } 
  else if (x >= 1.1 && x < 3.3) {    // Optical/NIR

    y2 = y  * y ;
    y3 = y2 * y ;
    y4 = y2 * y2 ;
    y5 = y3 * y2;
    y6 = y3 * y3;
    y7 = y4 * y3;
    y8 = y4 * y4;

    if ( DO94 ) {
    a = 1. + 0.104*y - 0.609*y2 + 0.701*y3 + 1.137*y4
      - 1.718*y5 - 0.827*y6 + 1.647*y7 -0.505*y8 ;

    b = 1.952*y + 2.908*y2 -3.989*y3 - 7.985*y4
      + 11.102*y5 + 5.491*y6 - 10.805*y7 + 3.347*y8;
    }
    else {
    a = 1. + 0.17699*y - 0.50447*y2 - 0.02427*y3
      + 0.72085*y4 + 0.01979*y5 - 0.77530*y6 + 0.32999*y7 ;

    b = 1.41338*y + 2.28305*y2 + 1.07233*y3 - 5.38434*y4
      - 0.62251*y5 + 5.30260*y6 - 2.09002*y7 ;
    }

  } 
  else if (x >= 3.3 && x < 8.0 ) {    // UV
    if (x >= 5.9) {
      xx  = x - 5.9 ;
      xx2 = xx  * xx ;
      xx3 = xx2 * xx ;

      fa = -0.04473*xx2 - 0.009779*xx3 ;
      fb =  0.21300*xx2 + 0.120700*xx3 ;

    } else {
      fa = fb = 0.0;
    }

    xx = x - 4.67 ; xx2 = (xx*xx);
    a =  1.752 - 0.316*x - 0.104/(xx2 + 0.341) + fa;

    xx = x - 4.62 ; xx2 = (xx*xx);
    b = -3.090 + 1.825*x + 1.206/(xx2 + 0.263) + fb;
  } 
  else if (x >= 8.0 && x <= 10.0) {  // Far-UV
    xx  = x - 8.0  ;
    xx2 = xx  * xx ;
    xx3 = xx2 * xx ; 

    a = -1.073 - 0.628*xx + 0.137*xx2 - 0.070*xx3 ;
    b = 13.670 + 4.257*xx - 0.420*xx2 + 0.374*xx3 ;
  } else {
    a = b = 0.0;
  }

  XT = AV*(a + b/RV);

  // Sep 18 2013 RK/DS - Check option for Fitzptrack 99

#define NPOLY_FITZ99 11 //Dillon and Dan upped to 10, Oct 9 2021
  if ( OPT == 99 ) {  

    double XTcor, wpow[NPOLY_FITZ99] ;

    // xxx mark delete double F99_over_O94[NPOLY_FITZ99] = {  // From D.Scolnic, Sep 18 2013
    //  0.485382, 0.791117, -0.534349, 0.191105,
    //  -0.0380031, 0.00416853,  -0.000235077, 5.31309e-06 
    //} ;
    
    double F99_over_O94[NPOLY_FITZ99] = {  // Dillon and Dan, Oct 9 2021
      8.55929205e-02,  1.91547833e+00, -1.65101945e+00,  7.50611119e-01,
      -2.00041118e-01,  3.30155576e-02, -3.46344458e-03,  2.30741420e-04,
      -9.43018242e-06,  2.14917977e-07, -2.08276810e-09
    };

    if ( WAVE > WAVEMAX_FITZ99  ) {
      sprintf(c1err,"Invalid WAVE=%.1f A for Fitzpatrick 99 color law.",
	      WAVE );
      sprintf(c2err,"Avoid NIR (>%.1f), or update Fitz99 in NIR",
	      WAVEMAX_FITZ99 );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // compute powers of wavelength without using slow 'pow' function
    wpow[0]  = 1.0 ;
    wpow[1]  = WAVE/1000. ;
    wpow[2]  = wpow[1] * wpow[1] ;
    wpow[3]  = wpow[1] * wpow[2] ;
    wpow[4]  = wpow[2] * wpow[2] ;
    wpow[5]  = wpow[3] * wpow[2] ;
    wpow[6]  = wpow[3] * wpow[3] ;
    wpow[7]  = wpow[4] * wpow[3] ;
    wpow[8]  = wpow[4] * wpow[4] ;
    wpow[9]  = wpow[5] * wpow[4] ;
    wpow[10]  = wpow[5] * wpow[5] ;

    XTcor = 0.0 ;
    for(i=0; i < NPOLY_FITZ99; i++ ) 
      {  XTcor += (wpow[i] * F99_over_O94[i]) ; }
    
    XT *= XTcor ;
  }

  return XT ;

}  // end of GALextinct



