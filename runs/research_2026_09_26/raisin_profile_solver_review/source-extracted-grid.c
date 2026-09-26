#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MXPATHLEN 1024
#define MXFILTINDX 100
#define SEV_FATAL 1
typedef void fitsfile;
char c1err[2048],c2err[2048];
void errmsg(int a,int b,char*c,char*d,char*e){fprintf(stderr,"%s: %s %s\n",c,d,e);abort();}
/***********
 Created Nov 2010 by R.Kessler

 
 parameters and global variables used to write/read 
 the GRID in fits format.

 Mar 16, 2016: add IFILTOBS[ifilt]
 Aug 27, 2017:  NON1A_NAME -> 200 chars instead of 40
 Sep 15, 2027:  MXGRIDGEN-> 500 (was 200)

**************/

// define internal version for backward-compatibility
//#define  IVERSION_GRID_WRITE 2    // Apr 3 2013
//#define  IVERSION_GRID_WRITE 3    // Aug 30 2013: add WGT and MAGOFF
#define  IVERSION_GRID_WRITE 4    // Jan 08 2017: add ITYPE_USER

#define MXGRIDGEN       500       // max number of grid point in any dimension
#define MAXMAG_GRIDGEN  32.0         // max mag to store with 16 bits
#define MAGPACK_GRIDGEN  1000.0       // I2MAG = MAG * MAGPACKSCALE
#define FLUXPACK_GRIDGEN 20000.      // for snoopy fluxes
#define NPADWD_LCBEGIN  2            // # pad words for each light curve
#define NPADWD_LCEND    2            // # pad words for each light curve
#define MARK_GRIDGEN_LCBEGIN -1111    // begin LC marker
#define MARK_GRIDGEN_LCEND   -9999    // end of LC marker

#define IPAR_GRIDGEN_LOGZ      1   // log10(redshift)
#define IPAR_GRIDGEN_COLORPAR  2   // AV or SALT2c
#define IPAR_GRIDGEN_COLORLAW  3   // RV or beta
#define IPAR_GRIDGEN_SHAPEPAR  4   // x1, Delta, dm15, nonIa-index ..
#define IPAR_GRIDGEN_FILTER    5   // filter index
#define IPAR_GRIDGEN_TREST     6   // rest-frame epoch
#define NPAR_GRIDGEN           6

#define EXTNAME_SNPAR_INFO   "SNPAR-INFO"
#define EXTNAME_LOGZ         "LOGZ-GRID" 
#define EXTNAME_COLORPAR     "COLOR-GRID" 
#define EXTNAME_COLORLAW     "RV/BETA-GRID" 
#define EXTNAME_SHAPEPAR     "SHAPE-GRID" 
#define EXTNAME_FILTER       "FILTER-GRID" 
#define EXTNAME_TREST        "TREST-GRID" 
#define EXTNAME_NONIa_INFO   "NONIa-INFO" 
#define EXTNAME_PTRI2LCMAG   "PTR_I2LCMAG" 
#define EXTNAME_I2LCMAG      "I2LCMAG" 
char    EXTNAME_GRIDGEN[NPAR_GRIDGEN+1][40] ;

#define SNTYPE_GRIDGEN_Ia      1
#define SNTYPE_GRIDGEN_NONIa   2
#define SNTYPE_GRIDGEN_SIMSED  3

#define OPT_GRIDGEN_FORMAT_TEXT 1
#define OPT_GRIDGEN_FORMAT_FITS 2
int     OPT_GRIDGEN_FORMAT ;

FILE     *fp_GRIDGEN_TEXT ;
fitsfile *fp_GRIDGEN_FITS ;

#define NCOMMENT_GRIDGEN  8
char COMMENT_GRIDGEN[NCOMMENT_GRIDGEN][80];


float GRIDGEN_I2LCPACK;   // stored I*2 value = mag(or flux) x this value


// define user inputs
struct GRIDGEN_INPUTS {
  int  NBIN[NPAR_GRIDGEN+1] ;
  char FORMAT[20];  // from sim-input file: "TEXT" or "FITS"
} GRIDGEN_INPUTS ;



typedef struct {

  int   IVERSION ;          // internal version
  char  UNIQUE_KEY[500];

  char  SURVEY[60];
  char  MODEL[MXPATHLEN];
  char  FILTERS[MXFILTINDX] ;
  int   IFILTOBS[MXFILTINDX] ;  // Mar 2016, computed on readback only

  // start with info vs. grid parameter
  int   NBIN[NPAR_GRIDGEN+1];                // Number of bins
  float BINSIZE[NPAR_GRIDGEN+1] ;            // bin size
  float VALUE[NPAR_GRIDGEN+1][MXGRIDGEN] ;   // value at each grid point
  char  NAME[NPAR_GRIDGEN+1][20];        // name of param (i.e, DELTA, AV)
  float VALMIN[NPAR_GRIDGEN+1] ;   // same as VALUE[1]
  float VALMAX[NPAR_GRIDGEN+1] ;   // same as VALUE[NBIN]

  float FILTER_LAMAVG[MXGRIDGEN];  // mean wavelength for each filter

  // ilc = 1 + sum ILCOFF * (indx[IPAR] -1)
  int ILCOFF[NPAR_GRIDGEN+1] ; // to get abs LC index from iz,ic,ilum, etc ...

  int   *PTR_VALUE[NPAR_GRIDGEN+1];    // VALUE[ilc] = value for this 'ilc'
  float CURRENT_VALUE[NPAR_GRIDGEN+1];   // current par values for each LC

  // now the LC info
  int NGRIDGEN_LC  ;    // total number of lightcurves to generate
  int NGRIDGEN_PER_LC ; // number of measures per LC

  long int SIZEOF_GRIDGEN ; // total size of GRIDGEN (excluding header)
  int      NWD_I2GRIDGEN;

  short *I2GRIDGEN_LCMAG ;
  short *I2GRIDGEN_LCERR ;
  int   *PTR_GRIDGEN_LC ;

  // extra NONIA-INFO from simgen-input file
  int   NON1A_INDEX[MXGRIDGEN];     // SNANA index vs. sparse nonIa index
  char  NON1A_NAME[MXGRIDGEN][200];  // full SN name vs. idem
  char  NON1A_CTYPE[MXGRIDGEN][20]; // string Type(Ib,II..) vs. idem
  int   NON1A_ITYPE_AUTO[MXGRIDGEN];  // 1=Ib,Ic,Ibc, etc ... 2=II, IIP, IIL..
  int   NON1A_ITYPE_USER[MXGRIDGEN];  // SNTAG = user type (Jan 2017)
  float NON1A_WGT[MXGRIDGEN];       // relative rate (Aug 30 2013)
  float NON1A_WGTSUM[MXGRIDGEN];    // sumulative sum
  float NON1A_MAGOFF[MXGRIDGEN];    // mag offset
  float NON1A_MAGSMEAR[MXGRIDGEN] ; // mag smear in sim (not used to make grid)

  // Aug 2016: define info for PEC1A to allow different z-dependent rate
  int   ISPEC1A[MXGRIDGEN];         // flag peculiar Ia (Aug 2016)
  float FRAC_PEC1A ;

} SNGRID_DEF ;

SNGRID_DEF  SNGRID_WRITE ; // used by sim to write GRID

int OPT_SNOOPY_FLUXPACK ;

int NROW_WRITE_TOT ;

// ==========================


// ===================================
//     GRID  function declarations
// ===================================

// grid-write utils
void   init_GRIDsource(int opt);     // init stuff for GRID
void   init0_GRIDsource(void);   
void   init1_GRIDsource(void);    // init after genmodel is init'ed

void   gen_GRIDevent(int ilc);
void   update_GRIDarrays(void);

void   wr_GRIDfile(int OPT, char *GRIDfile);
void   wrhead_GRIDfile_text(void);
void   wrhead_GRIDfile_fits(void);
void   append_GRIDfile_text(void);
void   append_GRIDfile_fits(void);

void   end_GRIDfile(void);
void   get_GRIDKEY(void) ;

int    SNTYPE_GRIDGEN(void); 
void   load_EXTNAME_GRIDGEN(void);  // load EXTNAME_GRIDGEN array

// read back utils
void   load_EXTNAME_GRIDREAD(int IVERSION); // idem for different ifdef 
void   fits_read_SNGRID(int OPTMASK, char *sngridFile, 
			SNGRID_DEF *SNGRID); 

void   check_fitserror(char *comment, int status);

int  INDEX_GRIDGEN(int ipar, double parval, SNGRID_DEF *SNGRID);
int  get_NON1A_ITYPE_SNGRID(char *NON1A_CTYPE );

void renorm_wgts_SNGRID(SNGRID_DEF *SNGRID ) ;
void sort_by_PEC1A_SNGRID(SNGRID_DEF *SNGRID ) ;
void copy_NON1A_SNGRID(int IROW1, int IROW2, 
		       SNGRID_DEF *SNGRID1, SNGRID_DEF *SNGRID2 );

void dump_SNGRID(SNGRID_DEF *SNGRID ) ;

// ============= END ==========

SNGRID_DEF SNGRID_SNOOPY;
void load_arrays(int ns,float*s,int nt,float*t,int*ptr,short*mag,short*err,float pack){
 memset(&SNGRID_SNOOPY,0,sizeof(SNGRID_SNOOPY));
 int pars[2]={IPAR_GRIDGEN_SHAPEPAR,IPAR_GRIDGEN_TREST};
 int sizes[2]={ns,nt};float*vals[2]={s,t};
 for(int q=0;q<2;q++){int p=pars[q],n=sizes[q];
  SNGRID_SNOOPY.NBIN[p]=n;
  for(int j=0;j<n;j++)SNGRID_SNOOPY.VALUE[p][j+1]=vals[q][j];
  SNGRID_SNOOPY.VALMIN[p]=vals[q][0];SNGRID_SNOOPY.VALMAX[p]=vals[q][n-1];
  float dif=vals[q][n-1]-vals[q][0];SNGRID_SNOOPY.BINSIZE[p]=dif/(float)(n-1);
 }
 SNGRID_SNOOPY.PTR_GRIDGEN_LC=ptr;
 SNGRID_SNOOPY.I2GRIDGEN_LCMAG=mag;SNGRID_SNOOPY.I2GRIDGEN_LCERR=err;
 GRIDGEN_I2LCPACK=pack;
}
int INDEX_GRIDGEN(int ipar, double parval, SNGRID_DEF *SNGRID) {

  // return index of ipar-parameter with value parval.
  double valmin, valmax, valbin, dif, ratio ;
  int indx ;
  // ---------- BEGIN
  indx = -9;
  valmin = (double)SNGRID->VALMIN[ipar] ;
  valmax = (double)SNGRID->VALMAX[ipar] ;
  valbin = (double)SNGRID->BINSIZE[ipar] ;

  if ( parval <= valmin ) 
    { indx = 1 ; }
  else if ( parval >= valmax ) 
    { indx = SNGRID->NBIN[ipar] ; }
  else  { 
    dif   = parval - valmin ;
    ratio = dif / valbin ;
    indx  = 1 + (int)ratio  ;  
  }

  return(indx) ;

} void gridinterp_snoopy(int ifilt, double shape, 
		       int nobs, double *Trest_list,
		       double *mag_list, double *magerr_list ) {


  // Interpolate mag and magerr on snoopy grid.
  // Note that fits_read_SNGRID() must be called before
  // calling this function.

  int  i, ioff, im, index_shape, index_Trest ;
  int  ILC, IPTRLC_OFF[2], IPTRLC, NBIN_TREST, NBIN_SHAPE ;    

  double  
    Trest, ratio_Trest, ratio_shape
    ,gridFlux2[2][2]     // [Trest][shape]
    ,gridFlux1[2]
    ,gridErr2[2][2]     // [Trest][shape]
    ,gridErr1[2]
    ,dif
    ,TREST_MIN,  SHAPE_MIN 
    ,TREST_MAX,  SHAPE_MAX
    ,TREST_BIN,  SHAPE_BIN
    ;

  short I2TMP, I2ERR ;

  char fnam[] = "genmag_snoopy";

  // ----------- BEGIN -------------

  NBIN_TREST    = SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_TREST] ;
  NBIN_SHAPE    = SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_SHAPEPAR] ;

  TREST_BIN  = SNGRID_SNOOPY.BINSIZE[IPAR_GRIDGEN_TREST] ;
  SHAPE_BIN  = SNGRID_SNOOPY.BINSIZE[IPAR_GRIDGEN_SHAPEPAR] ;

  // get shape grid-index

  index_shape = INDEX_GRIDGEN(IPAR_GRIDGEN_SHAPEPAR, shape, &SNGRID_SNOOPY );
  if ( index_shape >= NBIN_SHAPE ) 
    { index_shape-- ; }

  ILC           = index_shape; 
  IPTRLC_OFF[0] = SNGRID_SNOOPY.PTR_GRIDGEN_LC[ILC+0] ; 
  IPTRLC_OFF[1] = SNGRID_SNOOPY.PTR_GRIDGEN_LC[ILC+1] ;

  SHAPE_MIN   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][index_shape] ;
  SHAPE_MAX   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][index_shape+1] ; (void)SHAPE_MAX;
  ratio_shape = (shape -  SHAPE_MIN)/SHAPE_BIN ;

  // make sure that 1st word is BEGIN-LC marker
  I2TMP = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC_OFF[0]];
  if ( I2TMP != MARK_GRIDGEN_LCBEGIN ){
    sprintf(c1err,"First I*2 word of ILC=%d is %d .", ILC, I2TMP );
    sprintf(c2err,"But expected %d", MARK_GRIDGEN_LCBEGIN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // make sure that 2nd word is first 8 bits of ILC
  I2TMP = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC_OFF[0]+1];
  if ( I2TMP != ( ILC & 127 ) ){
    sprintf(c1err,"2nd I*2=%d  for ILC=%d, PTRLC_OFF=%d .", 
	    I2TMP, ILC, IPTRLC_OFF[0]+1  );
    sprintf(c2err,"But expected ILC&127 = %d", (ILC & 127) );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  ioff       = ifilt*NBIN_TREST + NPADWD_LCBEGIN-1 ;

  for ( i=0; i < nobs; i++ ) {

    Trest       = Trest_list[i];
    index_Trest = INDEX_GRIDGEN(IPAR_GRIDGEN_TREST,Trest, &SNGRID_SNOOPY );
    if ( index_Trest >= NBIN_TREST ) 
      { index_Trest-- ; }

    TREST_MIN   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][index_Trest] ;
    TREST_MAX   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][index_Trest+1] ; (void)TREST_MAX;
    ratio_Trest = (Trest - TREST_MIN)/TREST_BIN ;

    /*
    if ( ratio_Trest < 0.0 || ratio_Trest > 1.000000001 ) {
      sprintf(c1err, "Bad ratio_Trest = %f at Trest(%s)=%6.2f  dm15=%6.2f", 
	      ratio_Trest, SNOOPY_TEMPLATE[ifilt].FILTERNAME, Trest, dm15 );
      sprintf(c2err,"index_Trest=%d  TREST(MIN,MAX,BIN)=%6.2f,%6.2f,%6.2f", 
	      index_Trest, TREST_MIN, TREST_MAX, TREST_BIN );
      MADABORT(c1err,c2err);
    }
    */

    // get gridFlux2 at the 4 corners bounding Trest,dm15
    for ( im=0; im <=1; im++ ) {    // dm15 edges
      IPTRLC      = IPTRLC_OFF[im] + ioff + index_Trest;
      I2TMP       = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC] ;
      I2ERR       = SNGRID_SNOOPY.I2GRIDGEN_LCERR[IPTRLC] ;
      // xxx      if ( I2ERR > 999 ) { I2TMP = 30000; } // temp hack, May 2020
      gridFlux2[0][im]  = (double)I2TMP / GRIDGEN_I2LCPACK ;
      gridErr2[0][im]   = (double)I2ERR / GRIDGEN_I2LCPACK ;

      IPTRLC      = IPTRLC_OFF[im] + ioff + (index_Trest+1) ;
      I2TMP       = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC] ;
      I2ERR       = SNGRID_SNOOPY.I2GRIDGEN_LCERR[IPTRLC] ;
      // xxxx      if ( I2ERR > 999 ) { I2TMP = 30000; } // temp hack, May 2020
      gridFlux2[1][im]  = (double)I2TMP / GRIDGEN_I2LCPACK ;
      gridErr2[1][im]   = (double)I2ERR / GRIDGEN_I2LCPACK ;

      dif = gridFlux2[1][im] - gridFlux2[0][im] ;
      gridFlux1[im] = gridFlux2[0][im] + (dif * ratio_Trest);

      dif = gridErr2[1][im] - gridErr2[0][im] ;
      gridErr1[im] = gridErr2[0][im] + (dif * ratio_Trest);
    }

    dif = gridFlux1[1] - gridFlux1[0] ;
    mag_list[i]  = gridFlux1[0] + (dif * ratio_shape);

    dif = gridErr1[1] - gridErr1[0] ;
    magerr_list[i]  = gridErr1[0] + (dif * ratio_shape);

  } // i-loop over observations
  

} // end of gridinterp_snoopy

