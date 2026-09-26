#!/usr/bin/env python
import numpy as np
import pylab as plt
import f90nml
import argparse
import os
import glob
import snana
from raisin_cosmo import get_vpec
from scipy.optimize import minimize
from txtobj import txtobj
import copy
import cosmo
import scipy.stats
from scipy.stats import binned_statistic
from raisin_cosmo import ovdatamc
import yaml

_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)))

_cosmomc_batch = """#!/bin/bash
#SBATCH --job-name=RAISIN_stat
#SBATCH --time=34:00:00
###SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=1-2
#SBATCH --cpus-per-task=1
#SBATCH --partition=broadwl-lc
#SBATCH --output=/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_stat.log
#SBATCH --account=pi-rkessler
#SBATCH --mem=20GB

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module unload openmpi
module load intelmpi/5.1+intel-16.0
module load cfitsio/3
module load mkl
module load gcc/6.1
source /project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/code/plc_3.0/plc-3.01/bin/clik_profile.sh

PARAMS=`expr ${{SLURM_ARRAY_TASK_ID}} - 1`

INI_FILES=({} {})
DONE_FILES=(done_0.txt done_1.txt)

cd /scratch/midway2/rkessler/djones/cosmomc/chains_2015/
mpirun /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/cosmomc ${{INI_FILES[$PARAMS]}}
#mpirun $RAISIN_ROOT/CosmoMC/cosmomc ${{INI_FILES[$PARAMS]}}

if [ $? -eq 0 ]; then
    echo "SUCCESS" > ${{DONE_FILES[$PARAMS]}}
else
    echo "FAILURE" > ${{DONE_FILES[$PARAMS]}}
fi

"""

_cosmosis_batch = """#!/bin/bash -l
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=34:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=5000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --job-name=RAISIN_stat

mpirun -n 2 cosmosis --mpi {}
"""


_cosmomc_ini = """DEFAULT(batch2/JLA.ini)
DEFAULT(batch2/GAUSS.ini)
DEFAULT(batch2/common.ini)
accuracy_level = 1.35
high_accuracy_default = F
CMB_lensing = F
lmin_store_all_cmb = 0
stop_on_error = F

MPI_Converge_Stop = 0.002
MPI_Limit_Converge = 0.002
MPI_Limit_Converge_Err = 0.13
MPI_Max_R_ProposeUpdate = 0.03
indep_sample = 5

propose_matrix= planck_covmats/base_TT_lowTEB_plik.covmat
root_dir = chains/

action = 0
num_threads = 1

start_at_bestfit = F
feedback = 1
use_fast_slow = F
checkpoint = T

#VIN LETS TRY SIMPLE METROPOLIS WITH THAT LIKELIHOOD
#sampling_method = 1
#sampling_method=7 is a new fast-slow scheme good for Planck
sampling_method = 7
dragging_steps  = 3
propose_scale = 2

#these are just small speedups for testing
get_sigma8=F

param[wa]=0
param[w]=-0.995 -1.5 -0.5 0.001 0.001
param[omegak]=0
param[omegam]=0.3 0.2 0.4 0.002 0.002
compute_tensors=F
param[r]=0.0
param[calPlanck]=1
# VERY IMPORTANT TO FIX PARAMETERS WE DONT CONSTRAIN WITH OUR APPROXIMATION
param[tau]=0.078
param[logA]=3.090388
param[alpha_JLA]=0.14
param[beta_JLA]=3.1
file_root={}
jla_dataset={}"""

_cosmosis_ini="""
[runtime]
sampler = importance

[importance]
; Chain of input samples (Planck likelihoods in this case).
input = planck_samples/chain_p-TTTEEE-lowE_SNmag_wcdm.txt
; Number of samples to do between saving output.
nstep = 1
; Include the old likelihood in the old likelihood; i.e. P'=P*P_new.
add_to_likelihood = T

[output]
filename = raisin_samples/{}.txt
format = text
verbosity= debug

[pipeline]
; We use two likelihoods, the JLA (for high redshift) and
; Riess 2011 to anchor H0, which is otherwise degenerate
; with the nuisance parameter M
modules = consistency camb pantheon
values = raisin_values.ini
extra_output =
likelihoods = pantheon
; jla
quiet=T
debug=F
timing=F

[camb]
; For background-only data we do not need a full
; Boltzmann evaluation, just D(z), etc.
; Setting mode=background means we get this.
file = cosmosis-standard-library/boltzmann/camb/camb.so
mode=background
feedback=0

[pantheon]
file = /n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/cosmosis/cosmosis-standard-library/likelihood/pantheon/pantheon.py
data_file = /n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/cosmo/{}
covmat_file = /n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/cosmo/{}

[consistency]
file = cosmosis-standard-library/utility/consistency/consistency_interface.py
"""

_datasettext = """#Settings for the joint SNLS/SDSS data analysis
name = JLA
data_file = {}
pecz = 0
intrinsicdisp = 0
intrinsicdisp0 = 0
intrinsicdisp1 = 0
intrinsicdisp2 = 0
intrinsicdisp3 = 0
twoscriptmfit = F
#scriptmcut = 10.0
has_mag_covmat = {}
mag_covmat_file = {}
has_stretch_covmat = F
has_colour_covmat = F
has_mag_stretch_covmat = F
has_mag_colour_covmat = F
has_stretch_colour_covmat = F"""

_nir_nml = ['$RAISIN_ROOT/cosmo/fit/CSP_RAISIN.nml',
            '$RAISIN_ROOT/cosmo/fit/PS1_RAISIN.nml',
            '$RAISIN_ROOT/cosmo/fit/DES_RAISIN.nml']
_outdirs = ['$RAISIN_ROOT/cosmo/output/fit_nir_sys/CSP_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/PS1_RAISIN',
            '$RAISIN_ROOT/cosmo/output/fit_nir_sys/DES_RAISIN']
_data_dirs = ['$RAISIN_ROOT/cosmo/data/Photometry/CSPDR3_RAISIN',
              '$RAISIN_ROOT/cosmo/data/Photometry/PS1_RAISIN',
              '$RAISIN_ROOT/cosmo/data/Photometry/DES_RAISIN']

_sysgroupdict = {'photcal':('CSP_Y_SURVCAL','CSP_J_SURVCAL','CSP_H_SURVCAL','HST_CAL'),
                 'hstcal':('HST_CAL',),
                 'lowzcal':('CSP_Y_SURVCAL','CSP_J_SURVCAL','CSP_H_SURVCAL'),
                 'massdivide':('MASS_DIVIDE','MASS_STEP'),
                 'biascor':('BIASCOR_SHAPE_LOWZ','BIASCOR_AV_1','BIASCOR_SHAPE_HIGHZ','BIASCOR_AV_2'),
                 'biascorav':('BIASCOR_AV_1','BIASCOR_AV_2'),
                 'biascorst':('BIASCOR_SHAPE_LOWZ','BIASCOR_SHAPE_HIGHZ'),
                 'pecvel':('VPEC',),
                 'mwebv':('MWEBV',),
                 'kcor':('KCOR',),
                 'tmpl':('TMPL',),
                 'lcfitter':('LCFITTER',)}


_fitopt_dict = {'MWEBV':('MWEBV_SCALE 0.95','MWEBV_SCALE 0.95','MWEBV_SCALE 0.95'),
                'HST_CAL':('MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0',
                           'MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0',
                           'MAGOBS_SHIFT_ZP_PARAMS 0 0.00714 0'),
                'VPEC':('HEADER_OVERRIDE_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\'',
                        'HEADER_OVERRIDE_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\'',
                        'HEADER_OVERRIDE_FILE \'$RAISIN_ROOT/cosmo/vpec_sys_raisin.list\''),
                'MASS_DIVIDE':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'MASS_STEP':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'LCFITTER':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'CSP_Y_SURVCAL':('MAGOBS_SHIFT_ZP \'Y 0.03 y 0.03\'','->FITOPT000','->FITOPT000'),
                'CSP_J_SURVCAL':('MAGOBS_SHIFT_ZP \'J 0.02 j 0.02\'','->FITOPT000','->FITOPT000'),
                'CSP_H_SURVCAL':('MAGOBS_SHIFT_ZP \'H 0.02\'','->FITOPT000','->FITOPT000'),
                # Bomngri
                'CSP_B_SURVCAL':('MAGOBS_SHIFT_ZP \'B 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_V_SURVCAL':('MAGOBS_SHIFT_ZP \'o 0.01 m 0.01 n 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_g_SURVCAL':('MAGOBS_SHIFT_ZP \'g 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_r_SURVCAL':('MAGOBS_SHIFT_ZP \'r 0.01\'','->FITOPT000','->FITOPT000'),
                'CSP_i_SURVCAL':('MAGOBS_SHIFT_ZP \'i 0.01\'','->FITOPT000','->FITOPT000'),
                'PS1_g_SURVCAL':('->FITOPT000','MAGOBS_SHIFT_ZP \'g 0.003\'','->FITOPT000'),
                'PS1_r_SURVCAL':('->FITOPT000','MAGOBS_SHIFT_ZP \'r 0.003\'','->FITOPT000'),
                'PS1_i_SURVCAL':('->FITOPT000','MAGOBS_SHIFT_ZP \'i 0.003\'','->FITOPT000'),
                'PS1_z_SURVCAL':('->FITOPT000','MAGOBS_SHIFT_ZP \'z 0.003\'','->FITOPT000'),
                'DES_g_SURVCAL':('->FITOPT000','->FITOPT000','MAGOBS_SHIFT_ZP \'g 0.006\''),
                'DES_r_SURVCAL':('->FITOPT000','->FITOPT000','MAGOBS_SHIFT_ZP \'r 0.006\''),
                'DES_i_SURVCAL':('->FITOPT000','->FITOPT000','MAGOBS_SHIFT_ZP \'i 0.006\''),
                'DES_z_SURVCAL':('->FITOPT000','->FITOPT000','MAGOBS_SHIFT_ZP \'z 0.006\''),
                'BIASCOR_SHAPE_LOWZ':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_AV_1':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_SHAPE_HIGHZ':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'BIASCOR_AV_2':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'KCOR':('KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_CSPDR3_BD17_sys.fits\'',
                        'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_PS1MD_NIR_sys.fits\'',
                        'KCOR_FILE \'$RAISIN_ROOT/cosmo/kcor/kcor_DES_NIR_sys.fits\''),
                'TMPL':('->FITOPT000','->FITOPT000','->FITOPT000'),
                'RV':('INIVAL_RV 3.1','INIVAL_RV 3.1','INIVAL_RV 3.1')
}

_filters_fit_dict = {'CSP':{'optical':'Bomngri',
                            'opticalnir':'BomngriYyJjH',
                            'nir':'YyJjH'},
                     'PS1':{'optical':'griz',
                            'opticalnir':'grizJH',
                            'nir':'JH'},
                     'DES':{'optical':'griz',
                            'opticalnir':'grizJH',
                            'nir':'JH'}}

def writecov(covmat,covmatfile):

    fout = open(covmatfile,'w')
    print('%i'%np.shape(covmat)[0],file=fout)
    for i in range(np.shape(covmat)[0]):
        for j in range(np.shape(covmat)[1]):
            if i != j:
                print('%8.5e'%covmat[j,i],file=fout)
            else:
                print('%8.5e'%0,file=fout)
    fout.close()

class biascor:
    def __init__(self,options=None):

        self.nirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SIM',
                                 '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SIM',
                                 '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SIM']
        self.nirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.nirsimfitreslist]
        
        self.opticalnirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPTNIR_SIM',
                                        '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPTNIR_SIM',
                                        '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPTNIR_SIM']
        self.opticalnirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.opticalnirsimfitreslist]

        self.opticalsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPT_SIM',
                                     '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPT_SIM',
                                     '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPT_SIM']
        self.opticalsimfitreslist = [os.path.expandvars(filepath) for filepath in self.opticalsimfitreslist]

        

        self.nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
        self.nirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.nirdatafitreslist]

        self.opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
                                         '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
                                         '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']
        self.opticalnirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticalnirdatafitreslist]

        if options is not None:
            self.options=options
        
    def apply_all_cuts(self,fr,fropt,restrict_to_good_list=False):

        if not restrict_to_good_list:
            # AV
            iGoodAV = np.zeros(len(fr.CID),dtype=bool)
            for j,i in enumerate(fr.CID):
                if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
                    iGoodAV[j] = True

            # reasonable stretch
            iGoodSt = np.zeros(len(fr.CID),dtype=bool)
            for j,i in enumerate(fr.CID):
                if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.75 and fropt.STRETCH[fropt.CID == i] < 1.185: #175:
                    iGoodSt[j] = True

            for k in fr.__dict__.keys():
                fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]
                
        if restrict_to_good_list:
            # then case-by-case removal of bad stuff
            # because some RAISIN things have bad photometry
            iGood = np.array([],dtype=int)
            for j,i in enumerate(fr.CID):
                if i in _goodcids: iGood = np.append(iGood,j)
            for k in fr.__dict__.keys():
                fr.__dict__[k] = fr.__dict__[k][iGood]
            
        return fr

    def apply_biascor(self,fitopt,sys=None):
        if fitopt == 0:
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            fig2 = plt.figure()
        print(fitopt,sys)

        if self.options.version == 'nir':
            simfitreslist = self.nirsimfitreslist
        elif self.options.version == 'opticalnir':
            simfitreslist = self.opticalnirsimfitreslist
        elif self.options.version == 'optical':
            simfitreslist = self.opticalsimfitreslist
        datafitresdirs = [o.replace('_nir_','_'+self.options.version+'_') for o in _outdirs]

        # can interpolate for each SN individually because samples are so small
        for nirdatadir,opticalnirdatafitres,nirsimfitres,opticalnirsimfitres,name,idx in zip(
                datafitresdirs,self.opticalnirdatafitreslist,simfitreslist,self.opticalnirsimfitreslist,['CSP','PS1','DES'],range(4)):

            if name == 'CSP': simname = 'CSP_RAISIN_SIM'
            elif name == 'PS1': simname = 'PS1_RAISIN_SIM'
            elif name == 'DES': simname = 'DES_RAISIN_SIM'
            print(nirdatadir)
            frdata = txtobj(glob.glob(os.path.expandvars("%s/*/FITOPT%03i.FITRES.gz"%(nirdatadir,fitopt)))[0],fitresheader=True)
            fropt = txtobj(opticalnirdatafitres,fitresheader=True)
            with open(os.path.expandvars(f"{nirdatadir}/SUBMIT.INFO")) as fin:
                yml = yaml.load(fin, Loader=yaml.FullLoader)
                fitoptstr = [' '.join(f) for f in yml['FITOPT_LIST']]

            sys = fitoptstr[fitopt]
            if 'BIASCOR_AV_1' in sys: syskey = '_AVSYS'; print(syskey,fitopt)
            elif 'BIASCOR_SHAPE_LOWZ' in sys and 'CSP_RAISIN' in nirdatadir: syskey = '_STSYS'; print(syskey,fitopt)
            elif 'BIASCOR_AV_2' in sys: syskey = '_AVSYS_2'; print(syskey,fitopt)
            elif 'BIASCOR_SHAPE_HIGHZ' in sys and 'CSP_RAISIN' not in nirdatadir: syskey = '_STSYS'; print(syskey,fitopt)
            else: syskey = ''


            frsim = txtobj(glob.glob(os.path.expandvars(f'{nirsimfitres}/{simname}{syskey}/FITOPT000.FITRES.gz'))[0],
                           fitresheader=True)
            froptsim = txtobj(glob.glob(os.path.expandvars(f'{opticalnirsimfitres}/{simname}{syskey}/FITOPT000.FITRES.gz'))[0],
                              fitresheader=True)
            
            frdata = self.apply_all_cuts(frdata,fropt,restrict_to_good_list=True)
            fropt = self.apply_all_cuts(fropt,fropt,restrict_to_good_list=True)
            frsim = self.apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)
            froptsim = self.apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)

            def get_sigint(resid,residerr):
                chi2_redlist = []
                sigint_testlist = np.arange(0,0.3,0.005)
                for sigint_test in sigint_testlist:
                    chi2_redlist += [np.sum(resid**2./(residerr**2.+sigint_test**2.))/float(len(resid)-1)]
                sigint = sigint_testlist[np.abs(np.array(chi2_redlist)-1) == np.min(np.abs(np.array(chi2_redlist)-1))][0]
                return sigint
                
            sigint = get_sigint(frsim.DLMAG-cosmo.mu(frsim.zHD),frsim.DLMAGERR)
            optsigint = get_sigint(froptsim.DLMAG-cosmo.mu(froptsim.zHD),froptsim.DLMAGERR)
            frsim.DLMAGERR = np.sqrt(frsim.DLMAGERR**2. + sigint**2.)
            froptsim.DLMAGERR = np.sqrt(froptsim.DLMAGERR**2. + optsigint**2.)
            
            frdata.DLMAG_biascor = np.array([-99.]*len(frdata.CID))
            for j,i in enumerate(frdata.CID):
                if name not in ['PS1','DES']:
                    iBias = np.where(np.abs(frsim.zHD - frdata.zHD[j]) < 0.005)[0]
                    #print('hack for cepheids!!')
                    #iBias = np.where(np.abs(frsim.zHD - frdata.zHD[j]) < 0.01)[0]
                else: iBias = np.where(np.abs(frsim.zHD - frdata.zHD[j]) < 0.03)[0]
                if len(iBias) < 60: #95:
                    import pdb; pdb.set_trace()
                    raise RuntimeError('not enough biascor events!')
                bias_j = np.average(frsim.DLMAG[iBias]-frsim.SIM_DLMAG[iBias],weights=1/(frsim.DLMAGERR[iBias]**2.))
                frdata.DLMAG_biascor[j] = bias_j
                frdata.DLMAG[j] -= bias_j

            if fitopt == 0:
                zbins = np.linspace(np.min(frsim.zHD),np.max(frsim.zHD),10)
                biasbinned = binned_statistic(frsim.zHD,frsim.DLMAG-frsim.SIM_DLMAG,bins=zbins,statistic='median').statistic
                ax.plot((zbins[1:]+zbins[:-1])/2.,biasbinned,'o-',label=name)
                fig2.clf()
                hist = ovdatamc.ovhist()
                parser = hist.add_options(usage='')
                options,  args = parser.parse_args()
                hist.options = options
                hist.options.journal = True
                hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.75,1.3),('AV',-0.5,0.93)]
                hist.options.nbins = 10
                hist.options.clobber = True
                hist.options.outfile = 'figs/sim_%s.png'%name
                hist.options.histvar = ['MURES','AV','STRETCH']
                hist.main(data=fropt,sim=froptsim)
                
            frdata.writefitres(f"output/fitres_cosmo_{self.options.version}/{name}.FITRES",clobber=True)
            if idx == 0:
                frdata_combined = copy.deepcopy(frdata)
            else:
                if 'PKMJDINI' in frdata_combined.__dict__.keys(): del frdata_combined.__dict__['PKMJDINI']
                if 'ERRFLAG_FIT' in frdata_combined.__dict__.keys(): del frdata_combined.__dict__['ERRFLAG_FIT']
                for k in frdata_combined.__dict__.keys():
                    frdata_combined.__dict__[k] = np.append(frdata_combined.__dict__[k],frdata.__dict__[k])


        frdata_combined.writefitres(f'output/cosmo_fitres_{self.options.version}/RAISIN_combined_FITOPT%03i.FITRES'%fitopt,clobber=True)
        frdata_combined.DLMAG += frdata_combined.DLMAG_biascor
        frdata_combined.writefitres(f'output/cosmo_fitres_{self.options.version}/RAISIN_combined_FITOPT%03i_nobiascor.FITRES'%fitopt,clobber=True)
        if fitopt == 0:
            ax.legend()
            fig.savefig('figs/biascor_applied.png')
        #import pdb; pdb.set_trace()
        
class cosmo_sys:
    def __init__(self):
        pass

    def add_arguments(self,parser=None, usage=None, config=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

        # The basics
        parser.add_argument('-v', '--verbose', action="count", dest="verbose",
                            default=0,help='verbosity level')
        parser.add_argument('-n','--make_nml', default=False,action="store_true",
                            help='NML file with FITOPT')
        parser.add_argument('-c','--make_covmat', default=False, action="store_true",
                            help='covmat & cosmoMC inputs')
        parser.add_argument('--clobber', default=False,action="store_true",
                            help='clobber flag')
        parser.add_argument('--version', default='nir',type=str,
                            help='optical, opticalnir, or nir')

        return parser

    def get_vpec(self):
        with open(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_sys_raisin.list'),'w') as foutsys,\
             open(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_baseline_raisin.list'),'w') as foutbase:
            print("VARNAMES: CID VPEC VPEC_ERR",file=foutsys)
            print("VARNAMES: CID VPEC VPEC_ERR",file=foutbase)
            for d in _data_dirs:
                listfile = glob.glob(os.path.expandvars(f"{d}/*LIST"))[0]
                files = np.loadtxt(listfile,unpack=True,dtype='str')
                for f in files:
                    sn = snana.SuperNova(os.path.expandvars(f"{d}/{f}"))
                    if 'DECL' not in sn.__dict__.keys():
                        sn.DECL = sn.DEC
                    vpec,vpec_sys = get_vpec.main(float(sn.RA.split()[0]),float(sn.DECL.split()[0]),float(sn.REDSHIFT_FINAL.split('+-')[0]))
                    print(f"SN: {sn.SNID} {vpec+vpec_sys} 250",file=foutsys)
                    print(f"SN: {sn.SNID} {vpec} 250",file=foutbase)
                    
    def mk_nml(self):

        # peculiar velocity list
        if self.options.clobber or not os.path.exists(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_sys_raisin.list')) or not\
           os.path.exists(os.path.expandvars('$RAISIN_ROOT/cosmo/vpec_baseline_raisin.list')):
            self.get_vpec()
        
        for i,nml in enumerate(_nir_nml):
            nml = os.path.expandvars(nml)
            survey = nml.split('/')[-1].split('_')[0] ###'$RAISIN_ROOT/cosmo/fit/CSP_RAISIN.nml'
            with open(nml.replace('.nml',f'_{self.options.version}_sys.nml'),'w') as fout:
                print(f"OUTDIR: {_outdirs[i].replace('_nir_','_'+self.options.version+'_')}",file=fout)
                with open(nml) as fin:
                    for line in fin:
                        if not line.startswith('OUTDIR') and not line.startswith('APPEND_FITRES'):
                            if 'optical' in self.options.version and 'INISTP' in line:
                                print('!'+line.replace('\n',''),file=fout)
                            elif 'FILTLIST_FIT' in line:
                                filters = _filters_fit_dict[survey][self.options.version]
                                print(f"         FILTLIST_FIT = '{filters}'",file=fout)
                            else:
                                print(line.replace('\n',''),file=fout)
                print('',file=fout)
                for k in _fitopt_dict.keys():
                    # YSE kcor variants
                    # if 'KCOR' in k: continue
                    
                    print(f'FITOPT: [{k}] {_fitopt_dict[k][i]}',file=fout)

    def bias_correct(self):

        bc = biascor(options=self.options)
        with open(os.path.expandvars(f"{_outdirs[0]}/SUBMIT.INFO")) as fin:
            yml = yaml.load(fin, Loader=yaml.FullLoader)
            fitoptcount = len(yml['FITOPT_LIST'])
            #fitoptcount = 0
            #for line in fin:
            #    if line.startswith('FITOPT'): fitoptcount += 1
        for i in range(fitoptcount):
            bc.apply_biascor(i)

    def nuisance_params(self,fitopt):

        # simultaneously fit for mass step, 3 sigint bins and 3 distance bins
        fr = txtobj(f'output/cosmo_fitres_{self.options.version}/RAISIN_combined_FITOPT%03i.FITRES'%fitopt,fitresheader=True)
        fr.mures = fr.DLMAG - cosmo.mu(fr.zHD)
        fr.mures -= np.median(fr.mures)

        # need to apply PV terms
        #print('hack!')
        peczerr=0.00083
        zerr = peczerr*5.0/np.log(10)*(1.0+fr.zHD)/(fr.zHD*(1.0+fr.zHD/2.0))
        fr.DLMAGERR = np.sqrt(fr.DLMAGERR**2. + zerr**2. + 0.055**2.*fr.zHD**2.)

        
        with open(os.path.expandvars(f"{_outdirs[0]}/SUBMIT.INFO")) as fin:
            yml = yaml.load(fin, Loader=yaml.FullLoader)
            fitoptstr = [' '.join(f) for f in yml['FITOPT_LIST']]

        #with open(os.path.expandvars(f"{_outdirs[0]}/FITOPT.README")) as fin:
        #    fitoptstr = fin.readlines()[1:]

        sys = fitoptstr[fitopt]
        if 'MASS_DIVIDE' in sys:
            msteploc = 10.44
        else:
            msteploc = 10
        p_lm = np.zeros(len(fr.CID))-99.
        fr.HOST_LOGMASS_ERR = np.sqrt(fr.HOST_LOGMASS_ERR**2.+0.02**2.)
        for i in range(len(fr.CID)):
            if fr.HOST_LOGMASS_ERR[i] == 10: fr.HOST_LOGMASS_ERR[i] = 0.2

            if fr.HOST_LOGMASS[i] > msteploc:
                p_lm[i] = scipy.stats.norm.cdf(
                    msteploc,fr.HOST_LOGMASS[i],
                    fr.HOST_LOGMASS_ERR[i])*100.
            else:
                p_lm[i] = scipy.stats.norm.cdf(
                    msteploc,fr.HOST_LOGMASS[i],
                    fr.HOST_LOGMASS_ERR[i])*100.

        iCSP = fr.IDSURVEY == 5
        iPS1_midz = (fr.IDSURVEY == 15) & (fr.zHD < 0.4306)
        iPS1_highz = (fr.IDSURVEY == 15) & (fr.zHD >= 0.4306)
        iDES_midz = (fr.IDSURVEY == 10) & (fr.zHD < 0.4306)
        iDES_highz = (fr.IDSURVEY == 10) & (fr.zHD >= 0.4306)
        
        def neglnlikefunc(x,p_lm=None,mu_i=None,sigma_i=None):

            mu_lowz,mu_midz,mu_highz = x[0],x[1],x[2]
            sigint_csp,sigint_ps1,sigint_des = x[3],x[4],x[5]
            mass_step = x[6]
            
            # sigint split by sample, but distance split by redshift
            # each one with a low-mass and high-mass component
            loglike_csp = -np.sum(np.logaddexp(-(mu_i[iCSP]+mass_step-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                               np.log((1-0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.))),
                                               -(mu_i[iCSP]-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                               np.log((0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.)))))
            
            loglike_ps1_midz = -np.sum(np.logaddexp(-(mu_i[iPS1_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                                    np.log((1-0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.))),
                                               -(mu_i[iPS1_midz]-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                                    np.log((0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.)))))
            
            loglike_ps1_highz = -np.sum(np.logaddexp(-(mu_i[iPS1_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                                     np.log((1-0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.))),
                                               -(mu_i[iPS1_highz]-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                                     np.log((0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.)))))
            
            loglike_des_midz = -np.sum(np.logaddexp(-(mu_i[iDES_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                                    np.log((1-0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.))),
                                                    -(mu_i[iDES_midz]-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                                    np.log((0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.)))))
            
            loglike_des_highz = -np.sum(np.logaddexp(-(mu_i[iDES_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                                     np.log((1-0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.))),
                                                     -(mu_i[iDES_highz]-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                                     np.log((0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.)))))

            return loglike_csp + loglike_ps1_midz + loglike_ps1_highz + loglike_des_midz + loglike_des_highz
        
        md = minimize(neglnlikefunc,[0,0.01,0.02,0.09,0.1,0.11,0.1],args=(p_lm,fr.mures,fr.DLMAGERR))


        # apply sigint and host mass, along with lensing uncertainty
        # write the new fitres files
        fr.DLMAGERR[fr.IDSURVEY == 15] = np.sqrt(fr.DLMAGERR[fr.IDSURVEY == 15]**2. + md.x[4]**2.)
        fr.DLMAGERR[fr.IDSURVEY == 5] = np.sqrt(fr.DLMAGERR[fr.IDSURVEY == 5]**2. + md.x[3]**2.)
        fr.DLMAGERR[fr.IDSURVEY == 10] = np.sqrt(fr.DLMAGERR[fr.IDSURVEY == 10]**2. + md.x[5]**2.)

        
        fr.MASS_CORR = np.zeros(len(fr.DLMAG))
        if 'MASS_STEP' not in sys:
            fr.DLMAG[fr.HOST_LOGMASS > 10] += md.x[6]/2.
            fr.DLMAG[fr.HOST_LOGMASS <= 10] -= md.x[6]/2.
            fr.MASS_CORR[fr.HOST_LOGMASS > 10] = md.x[6]/2.
            fr.MASS_CORR[fr.HOST_LOGMASS <= 10] = -md.x[6]/2.
        else:
            fr.DLMAG[fr.HOST_LOGMASS > 10] += (md.x[6]+np.sqrt(md.hess_inv[6,6]))/2.
            fr.DLMAG[fr.HOST_LOGMASS <= 10] -= (md.x[6]+np.sqrt(md.hess_inv[6,6]))/2.
            fr.MASS_CORR[fr.HOST_LOGMASS > 10] = (md.x[6]+np.sqrt(md.hess_inv[6,6]))/2.
            fr.MASS_CORR[fr.HOST_LOGMASS <= 10] = -(md.x[6]+np.sqrt(md.hess_inv[6,6]))/2.
        print(f'mass step = {md.x[6]:.3f} +/- {np.sqrt(md.hess_inv[6,6]):.3f}')
        fr.writefitres(f'output/cosmo_fitres_{self.options.version}/RAISIN_combined_FITOPT%03i_new.FITRES'%fitopt,clobber=True) 
        print(msteploc)
        
    def sys_covmat(self):
        syslist = ['stat','all','photcal','hstcal','lowzcal',
                   'kcor',
                   'massdivide',
                   'biascor','biascorav','biascorst','pecvel','mwebv','tmpl','lcfitter']
        for sys in syslist:
            count = 0
            with open(os.path.expandvars(f"{_outdirs[0]}/SUBMIT.INFO")) as fin:
                yml = yaml.load(fin, Loader=yaml.FullLoader)
                syslines = [' '.join(f) for f in yml['FITOPT_LIST']]

            for sysline,i in zip(syslines,range(len(syslines))):
                
                sysname = sysline.split()[1] #split('[')[-1].split(']')[0]
                if sys == 'all' or sys == 'stat' or sysname in _sysgroupdict[sys]:
                    if count == 0:
                        baselinefitres = f'output/cosmo_fitres_{self.options.version}/RAISIN_combined_FITOPT%03i_new.FITRES'%0
                        frbase = txtobj(baselinefitres,fitresheader=True)
                        covshape = len(frbase.CID)
                        
                        basecov = np.zeros([covshape,covshape])
                        for j in range(covshape):
                            basecov[j,j] = frbase.DLMAGERR[j]**2.

                        outcov = basecov[:]
                        if sys == 'stat': break
                    count += 1
                    
                    fr = txtobj(f'output/cosmo_fitres_{self.options.version}/RAISIN_combined_FITOPT%03i_new.FITRES'%i,fitresheader=True)

                    try:
                        global_off = np.average(fr.DLMAG-frbase.DLMAG,weights=1/(fr.DLMAGERR**2.))
                    except:
                        import pdb; pdb.set_trace()
                    dm2=(fr.DLMAG-cosmo.mu(fr.zHD))-global_off-(frbase.DLMAG-cosmo.mu(frbase.zHD))
                    dm2t=np.matrix(dm2)
                    dm2t=dm2t.T
                    dmm=dm2t*np.matrix(dm2)
                    outcov = np.add(outcov,dmm)

            writecov(outcov,f'output/cosmo_fitres_{self.options.version}/RAISIN_%s.covmat'%sys)

            with open(f'output/cosmo_fitres_{self.options.version}/RAISIN_%s_lcparams_cosmosis.txt'%sys,'w') as fout:
                print('# name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor snana',file=fout)
                for i in range(covshape):
                    print('0.0 %.6f %.6f 0.000000 %.6f %.6f 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.0000e+00 0.0000e+00 0.0000e+00 0.0000 0.0000000 0.0000000 0.0000'%(
                        frbase.zHD[i],frbase.zHD[i],frbase.DLMAG[i]-19.36,np.sqrt(outcov[i,i])),file=fout)
                fout.close()

    def cosmomc_inputs(self,dosubmit=False):
        syslist = ['stat','all','photcal','hstcal','lowzcal',
                   'massstep','massdivide',
                   'biascor','biascorav','biascorst','pecvel','mwebv','kcor','tmpl','lcfitter']

        for sys in syslist:
            batchfile = f'cosmomc/RAISIN_{sys}.sbatch'
            inifile = f'cosmomc/RAISIN_{sys}.ini'
            datasetfile = f'cosmomc/RAISIN_{sys}.dataset'
            lcparfile = f'output/cosmo_fitres_{self.options.version}/RAISIN_{sys}_lcparams.txt'
            covfile = f'output/cosmo_fitres_{self.options.version}/RAISIN_{sys}.covmat'
            root = f'raisin_{sys}'
            
            # batch file
            with open(batchfile,'w') as fout:
                print(_cosmomc_batch.format(inifile,inifile),file=fout)
            
            # ini file
            with open(inifile,'w') as fout:
                print(_cosmomc_ini.format(root,datasetfile),file=fout)
            
            # dataset file
            has_covmat = 'T'
            with open(datasetfile,'w') as fout:
                print(_datasettext.format(lcparfile,has_covmat,covfile,sys),file=fout)
            
            # lcparams and covmats should already be written

            if dosubmit:
                os.system(f'qsub {batchfile}')

    def cosmosis_inputs(self,dosubmit=False):
        syslist = ['stat','all','photcal','hstcal','lowzcal',
                   'massstep','massdivide',
                   'biascor','biascorav','biascorst','pecvel','mwebv','kcor','tmpl','lcfitter']

        for sys in syslist:
            batchfile = f'cosmosis_{self.options.version}/RAISIN_{sys}.sbatch'
            inifile = f'cosmosis_{self.options.version}/RAISIN_{sys}.ini'
            lcparfile = f'output/cosmo_fitres_{self.options.version}/RAISIN_{sys}_lcparams_cosmosis.txt'
            covfile = f'output/cosmo_fitres_{self.options.version}/RAISIN_{sys}.covmat'
            root = f'raisin_{sys}'
            
            # batch file
            with open(batchfile,'w') as fout:
                print(_cosmosis_batch.format(inifile.replace(f'cosmosis_{self.options.version}/','')),file=fout)
            
            # ini file
            with open(inifile,'w') as fout:
                print(_cosmosis_ini.format(root,lcparfile,covfile),file=fout)
            
            if dosubmit:
                os.system(f'sbatch {batchfile}')

                
    def mk_cosmo_inputs(self):

        # how many fit options?
        with open(os.path.expandvars(f"{_outdirs[0]}/SUBMIT.INFO")) as fin:
            yml = yaml.load(fin, Loader=yaml.FullLoader)
            fitoptcount = len(yml['FITOPT_LIST'])

        # copy over the template systematics fitting files

        datafitresdirs = [o.replace('_nir_','_'+self.options.version+'_') for o in _outdirs]
        for o in datafitresdirs:
            with open(os.path.expandvars(f"{o}/SUBMIT.INFO")) as fin:
                yml = yaml.load(fin, Loader=yaml.FullLoader)
                fitoptstr = [' '.join(f) for f in yml['FITOPT_LIST']]

            for line in fitoptstr:
                print(line)
                if 'TMPL' in line:
                    fitopt = line.split()[0].replace('FITOPT','')
                    # no point worrying about templates for opt+nir or opt only
                    if self.options.version != 'nir': continue

                    if f'fit_{self.options.version}_sys/PS1_RAISIN' in o: surveydir = 'PS1_RAISIN'
                    elif f'fit_{self.options.version}_sys/DES_RAISIN' in o: surveydir = 'DES_RAISIN'
                    else: 
                        print(o)
                        continue

                    fr = txtobj(f"output/fit_nir/{surveydir}_TMPLSYS.FITRES.TEXT",fitresheader=True)
                    fr0 = txtobj(f"{os.path.expandvars(o)}/{surveydir}/FITOPT000.FITRES.gz",fitresheader=True)
                    idx = np.array([],dtype=int)
                    for j,i in enumerate(fr0.CID):
                        idx = np.append(idx,np.where(fr.CID == fr0.CID[j])[0][0])
                    for k in fr.__dict__.keys():
                        fr.__dict__[k] = fr.__dict__[k][idx]
                    fr.writefitres(f"{os.path.expandvars(o)}/{surveydir}/FITOPT{fitopt}.FITRES")

                    os.chdir(f"{os.path.expandvars(o)}/{surveydir}/")
                    os.system(f"rm FITOPT{fitopt}.FITRES.gz")
                    os.system(f"gzip FITOPT{fitopt}.FITRES")
                    os.chdir("../../../../")
                if 'LCFITTER' in line:
                    fitopt = line.split()[0].replace('FITOPT','')
                    if self.options.version != 'nir': continue
                    if 'fit_nir_sys/PS1_RAISIN' in o:
                        lcfittingfile = 'output/BayeSN/RAISIN_BayeSN_theta_-1_NIR_dists_jones_tmax.txt'
                        surveydir = 'PS1_RAISIN'
                    elif 'fit_nir_sys/DES_RAISIN' in o:
                        lcfittingfile = 'output/BayeSN/RAISIN_BayeSN_theta_-1_NIR_dists_jones_tmax.txt'
                        surveydir = 'DES_RAISIN'
                    elif 'fit_nir_sys/CSP_RAISIN' in o:
                        lcfittingfile = 'output/BayeSN/CSP_RAISIN_BayeSN_NIR_theta_-1_dists_jones_snoopy_tmax.txt'
                        surveydir = 'CSPDR3_RAISIN'
                    else: 
                        raise RuntimeError(f"confused by directory {o}")
                        continue

                    fr = txtobj(lcfittingfile)
                    fr0 = txtobj(f"{os.path.expandvars(o)}/{surveydir}/FITOPT000.FITRES.gz",fitresheader=True)
                    fro = txtobj(f"{os.path.expandvars(o.replace('fit_nir_sys','fit_nir_sys_oldpkmjds'))}/{surveydir}/FITOPT000.FITRES.gz",fitresheader=True)
                    idx,idx2 = np.array([],dtype=int),np.array([],dtype=int)
                    for j,i in enumerate(fr0.CID):
                        if i in fr.sn and i in fro.CID:
                            fr0.DLMAG[j] = fr0.DLMAG[j] + (fr.__dict__['mu_mu+eta'][fr.sn == i][0]-fro.DLMAG[fro.CID == i][0])
                            fr0.DLMAGERR[j] = fr.__dict__['muerr_mu+eta'][fr.sn == i][0]
                        elif i in fr.sn and i not in fro.CID:
                            fr0.DLMAG[j] = fr.__dict__['mu_mu+eta'][fr.sn == i][0]
                            fr0.DLMAGERR[j] = fr.__dict__['muerr_mu+eta'][fr.sn == i][0]
                        elif i not in fr.sn:
                            fr0.DLMAG[j] = fr0.DLMAG[j]
                            fr0.DLMAGERR[j] = fr0.DLMAGERR[j]

                    fr0.writefitres(f"{os.path.expandvars(o)}/{surveydir}/FITOPT{fitopt}.FITRES")
                    os.chdir(f"{os.path.expandvars(o)}/{surveydir}/")
                    os.system(f"rm FITOPT{fitopt}.FITRES.gz")
                    os.system(f"gzip FITOPT{fitopt}.FITRES")
                    os.chdir("../../../../")


        # make files for every FITOPT, bias-correct, compute sigint, apply mass step
        self.bias_correct()
        for i in range(fitoptcount):
            self.nuisance_params(i)
        
        # sys covmat from FITRES files
        # covmat for each individual systematic shift also
        self.sys_covmat()

        # cosmomc inputs
        #self.cosmomc_inputs()
        self.cosmosis_inputs()
        
    def main(self):
        
        if self.options.make_nml:
            self.mk_nml()
        elif self.options.make_covmat:
            self.mk_cosmo_inputs()
        else:
            raise RuntimeWarning('doing nothing!')
        
if __name__ == "__main__":

    cs = cosmo_sys()
    parser = cs.add_arguments()
    args = parser.parse_args()
    cs.options = args
    cs.main()
