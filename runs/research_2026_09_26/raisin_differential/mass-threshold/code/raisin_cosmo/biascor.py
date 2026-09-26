#!/usr/bin/env python
# D. Jones - 2/27/20

import numpy as np
from scipy.stats import binned_statistic
from txtobj import txtobj
import pylab as plt
plt.ion()
import glob
import snana
import os
import copy
_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_CIDS.LIST',unpack=True,dtype=str),
                            #np.loadtxt('output/goodcids/CfA_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_CIDS.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_CIDS.LIST',unpack=True,dtype=str)))

from raisin_cosmo import ovdatamc

# rsync -avz output/sim/PS1_RAISIN_SIM_NIR david@turanga.ucsc.edu:/Users/David/Dropbox/research/RAISIN/cosmo/output/sim/

class biascor:

    def __init__(self):

        self.nirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN_NIR_SIM',
                                 #'$RAISIN_ROOT/cosmo/output/fit_nir/CfA_RAISIN_NIR_SIM',
                                 '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN_NIR_SIM',
                                 '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN_NIR_SIM']
        self.nirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.nirsimfitreslist]
        
        self.opticalnirsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPTNIR_SIM',
                                        #'$RAISIN_ROOT/cosmo/output/fit_all/CfA_RAISIN_OPTNIR_SIM',
                                        '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPTNIR_SIM',
                                        '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPTNIR_SIM']
        self.opticalnirsimfitreslist = [os.path.expandvars(filepath) for filepath in self.opticalnirsimfitreslist]

        self.opticalsimfitreslist = ['$RAISIN_ROOT/cosmo/output/fit_all/CSP_RAISIN_OPT_SIM',
                                     #'$RAISIN_ROOT/cosmo/output/fit_all/CfA_RAISIN_OPT_SIM',
                                     '$RAISIN_ROOT/cosmo/output/fit_all/PS1_RAISIN_OPT_SIM',
                                     '$RAISIN_ROOT/cosmo/output/fit_all/DES_RAISIN_OPT_SIM']
        self.opticalsimfitreslist = [os.path.expandvars(filepath) for filepath in self.opticalsimfitreslist]

        
        self.nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                                  #'$RAISIN_ROOT/cosmo/output/fit_nir/CfA_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
        self.nirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.nirdatafitreslist]

        self.opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
                                         #'$RAISIN_ROOT/cosmo/output/fit_optical/CfA_RAISIN_optnir.FITRES.TEXT',
                                         '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
                                         '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']
        self.opticalnirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticalnirdatafitreslist]

        self.opticaldatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT',
                                      #'$RAISIN_ROOT/cosmo/output/fit_optical/CfA_RAISIN_optical.FITRES.TEXT',
                                      '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',
                                      '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optical.FITRES.TEXT']
        self.opticaldatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticaldatafitreslist]

        
        self.outfitres = 'output/fitres/RAISIN_stat.fitres'
        self.figdir = 'figs'
        
    def add_options(self):
        pass

    def apply_salt2_cuts(self,fr,salt2alpha=0.147,salt2beta=3.13,zmin=None,zmax=None,fitprobmin=0.001,trestmax=5):

        if not zmin: zmin = np.min(fr.zHD)
        if not zmax: zmax = np.max(fr.zHD)

        sf = -2.5/(fr.x0*np.log(10.0))
        invvars = 1./(fr.mBERR**2.+ salt2alpha**2. * fr.x1ERR**2. + \
                          salt2beta**2. * fr.cERR**2. +  2.0 * salt2alpha * (fr.COV_x1_x0*sf) - \
                          2.0 * salt2beta * (fr.COV_c_x0*sf) - \
                          2.0 * salt2alpha*salt2beta * (fr.COV_x1_c) )

        try:
            cols = np.where((fr.x1 > -3.0) & (fr.x1 < 3.0) &
                (fr.c > -0.3) & (fr.c < 0.3) &
                            (fr.x1ERR < 1) & (fr.PKMJDERR < 2*(1+fr.zHD)) &
                (fr.FITPROB >= fitprobmin) &
                            (invvars > 0) & (fr.zHD >= zmin) &
                (fr.zHD <= zmax) & (fr.TrestMAX > trestmax))
        except:
            print('Warning : Keyword TrestMAX not found!!!')
            cols = np.where((fr.x1 > -3.0) & (fr.x1 < 3.0) &
                            (fr.c > -0.3) & (fr.c < 0.3) &
                            (fr.x1ERR < 1) & (fr.PKMJDERR < 2*(1+fr.zHD)) &
                            (fr.FITPROB >= fitprobmin) &
                            (invvars > 0) & (fr.zHD >= zmin) &
                            (fr.zHD <= zmax))

        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][cols]

        return(fr)
        
    def apply_all_cuts(self,fr,fropt,restrict_to_good_list=False):

        # AV
        iGoodAV = np.zeros(len(fr.CID),dtype=bool)
        for j,i in enumerate(fr.CID):
            if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
                iGoodAV[j] = True

        # reasonable stretch
        iGoodSt = np.zeros(len(fr.CID),dtype=bool)
        for j,i in enumerate(fr.CID):
            if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.8 and fropt.STRETCH[fropt.CID == i] < 1.175:
                iGoodSt[j] = True

        #iGoodSt = (fropt.STRETCH > 0.8) & (fropt.STRETCH < 1.3)

        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]

        #restrict_to_good_list = False
        if restrict_to_good_list:
            # then case-by-case removal of bad stuff
            # because some RAISIN things have bad photometry
            iGood = np.array([],dtype=int)
            for j,i in enumerate(fr.CID):
                if i in _goodcids: iGood = np.append(iGood,j)
            for k in fr.__dict__.keys():
                fr.__dict__[k] = fr.__dict__[k][iGood]
            
        return fr

    def apply_biascor(self):
        # can interpolate for each SN individually because samples are so small
        for nirdatafitres,opticalnirdatafitres,nirsimfitres,opticalnirsimfitres,name,idx in zip(
                self.nirdatafitreslist,self.opticalnirdatafitreslist,
                self.nirsimfitreslist,self.opticalnirsimfitreslist,['CSP','PS1','DES'],range(4)):

            frdata = txtobj(nirdatafitres,fitresheader=True)
            fropt = txtobj(opticalnirdatafitres,fitresheader=True)
            frsim = txtobj(glob.glob('%s/*/FITOPT000.FITRES'%nirsimfitres)[0],
                           fitresheader=True)
            froptsim = txtobj(glob.glob('%s/*/FITOPT000.FITRES'%opticalnirsimfitres)[0],
                              fitresheader=True)

            
            frdata = self.apply_all_cuts(frdata,fropt,restrict_to_good_list=True)
            frsim = self.apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)
            
            frdata.DLMAG_biascor = np.array([-99.]*len(frdata.CID))
            for j,i in enumerate(frdata.CID):
                if name not in ['PS1','DES']: iBias = np.where(np.abs(frsim.zCMB - frdata.zCMB[j]) < 0.015)[0]
                else: iBias = np.where(np.abs(frsim.zCMB - frdata.zCMB[j]) < 0.03)[0]
                if len(iBias) < 95:
                    import pdb; pdb.set_trace()
                    raise RuntimeError('not enough biascor events!')
                bias_j = np.average(frsim.DLMAG[iBias]-frsim.SIM_DLMAG[iBias],weights=1/(frsim.DLMAGERR[iBias]))
                frdata.DLMAG_biascor[j] = bias_j
                frdata.DLMAG[j] -= bias_j

                if frdata.HOST_LOGMASS[j] < 10:
                    frdata.DLMAG[j] -= 0.04
                    #frdata.DLMAGERR[j] = np.sqrt(frdata.DLMAGERR[j]**2. + 0.05**2.)

                # intrinsic dispersion, approx for now
                frdata.DLMAGERR[j] = np.sqrt(frdata.DLMAGERR[j]**2. + 0.15**2.)

            frdata.writefitres(f"output/fitres_cosmo/{name}.FITRES",clobber=True)
            if idx == 0:
                frdata_combined = copy.deepcopy(frdata)
            else:
                for k in frdata_combined.__dict__.keys():
                    frdata_combined.__dict__[k] = np.append(frdata_combined.__dict__[k],frdata.__dict__[k])
            #import pdb; pdb.set_trace()
        frdata_combined.writefitres('output/cosmo_fitres/RAISIN_combined.FITRES',clobber=True)

        cosmomcfile = '$RAISIN_ROOT/cosmo/output/cosmo_fitres/RAISIN_combined_stat.cosmomc.txt'
        cosmomcini = '$RAISIN_ROOT/cosmo/cosmomc/RAISIN_combined_stat.ini'
        cosmomcds = '$RAISIN_ROOT/cosmo/cosmomc/RAISIN_combined_stat.dataset'
        cosmomcbatch = '$RAISIN_ROOT/cosmo/cosmomc/RAISIN_combined_stat.job'
        self.write_cosmomc_snoopy(frdata_combined,cosmomcfile)
        self.submit_batch(cosmomcfile,cosmomcds,cosmomcini,cosmomcbatch)

    def submit_batch(self,cosmomcfile,cosmomcds,cosmomcini,cosmomcbatch):

        initext = f"""#DEFAULT(batch3/BAO.ini)
DEFAULT(batch3/JLA.ini)
# DEFAULT(batch3/GAUSS.ini)
#high-L plik likelihood
DEFAULT(batch3/plik_rd12_HM_v22_TTTEEE.ini)

#low-L temperature
DEFAULT(batch3/lowl.ini)

#low-L EE polarization
DEFAULT(batch3/simall_EE.ini)

>>>>>>> 684288af150e443a89355c7a0c467976df096784
DEFAULT(/project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/CosmoMC/batch3/common.ini)
INCLUDE(/project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/cosmo/cosmomc/base.ini)

MPI_Converge_Stop = 0.01
MPI_Limit_Converge = 0.01
MPI_Limit_Converge_Err = 0.185

#propose_matrix= /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/planck_covmats/base_TT_lowTEB_plik.covmat

param[wa]=0
param[w]=-0.995 -2. 0. 0.001 0.001
param[omegak]=0
param[omegam]=0.3 0.2 0.4 0.002 0.002
compute_tensors=F
param[r]=0.0
param[calPlanck]=1
action = 0

file_root=RAISIN_stat
jla_dataset={os.path.expandvars(cosmomcds)}
root_dir = /scratch/midway2/rkessler/djones/cosmomc/chains/
"""

        batchtext = f"""#!/bin/bash
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

INI_FILES=({os.path.expandvars(cosmomcini)} {os.path.expandvars(cosmomcini)})
DONE_FILES=(done_0.txt done_1.txt)

cd /scratch/midway2/rkessler/djones/cosmomc/chains/
#mpirun /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/cosmomc ${{INI_FILES[$PARAMS]}}
mpirun $RAISIN_ROOT/CosmoMC/cosmomc ${{INI_FILES[$PARAMS]}}

if [ $? -eq 0 ]; then
    echo "SUCCESS" > ${{DONE_FILES[$PARAMS]}}
else
    echo "FAILURE" > ${{DONE_FILES[$PARAMS]}}
fi

"""

        datasettext = f"""#Settings for the joint SNLS/SDSS data analysis
name = JLA
data_file = {os.path.expandvars(cosmomcfile)}
pecz = 0
intrinsicdisp = 0
intrinsicdisp0 = 0
intrinsicdisp1 = 0
intrinsicdisp2 = 0
intrinsicdisp3 = 0
twoscriptmfit = F
#scriptmcut = 10.0
has_mag_covmat = F
mag_covmat_file = cosmomc_data/foundps1_all.covmat
has_stretch_covmat = F
has_colour_covmat = F
has_mag_stretch_covmat = F
has_mag_colour_covmat = F
has_stretch_colour_covmat = F"""
                

        with open(os.path.expandvars(cosmomcds),'w') as fout:
            print(datasettext,file=fout)
        
        with open(os.path.expandvars(cosmomcini),'w') as fout:
            print(initext,file=fout)

        with open(os.path.expandvars(cosmomcbatch),'w') as fout:
            print(batchtext,file=fout)


        # submit the job
        print(f"sbatch {cosmomcbatch}")
        os.system(f"sbatch {cosmomcbatch}")
            
    def write_cosmomc_snoopy(self,fr,outfile):
        with open(os.path.expandvars(outfile),'w') as fout:
            print('# name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor snana',file=fout)
            for i in range(len(fr.CID)):
                print(f'0.0 {fr.zHD[i]:.6f} {fr.zHD[i]:.6f} 0.000000 {fr.DLMAG[i]-19.3:.6f} {fr.DLMAGERR[i]:.6f} 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.0000e+00 0.0000e+00 0.0000e+00 0.0000 0.0000000 0.0000000 0.0000 0',file=fout)

            
    def mk_sim_validplots(self):
        # make sure sim/data line up for all three sims
        # will have to worry about systematics down the road
        plt.rcParams['figure.figsize'] = (16,4)
        for simfitreslist,datafitreslist,opticalnirsimfitreslist,opticalnirdatafitreslist,label in \
            zip([self.nirsimfitreslist,self.opticalnirsimfitreslist],
                [self.nirdatafitreslist,self.opticalnirdatafitreslist],
                [self.opticalnirsimfitreslist,self.opticalnirsimfitreslist],
                [self.opticalnirdatafitreslist,self.opticalnirdatafitreslist],
                ['NIR','Optical']):
            if label == 'NIR': continue
            for i,survey in enumerate(['CSP','PS1','DES']):
                plt.clf()
                hist = ovdatamc.ovhist()
                parser = hist.add_options(usage='')
                options,  args = parser.parse_args()
                hist.options = options
                hist.options.journal = True
                hist.options.cutwin = [('MURES',-2,2),('STRETCH',0.8,1.175),('AV',-0.5,0.93)]
                hist.options.nbins = 10
                hist.options.clobber = True
                hist.options.outfile = 'figs/sim_%s_%s.png'%(survey,label)

                # 'SNRMAX1',
                hist.options.histvar = ['MURES','AV','STRETCH']
                print(datafitreslist[i],simfitreslist[i])
                fr = txtobj(datafitreslist[i],fitresheader=True)
                fropt = txtobj(self.opticaldatafitreslist[i],fitresheader=True)
                print(simfitreslist[i])
                frsim = txtobj(glob.glob('%s/*/FITOPT000.FITRES'%simfitreslist[i])[0],
                               fitresheader=True)
                froptsim = txtobj(glob.glob('%s/*/FITOPT000.FITRES'%self.opticalsimfitreslist[i])[0],
                                  fitresheader=True)
                
                fr = self.apply_all_cuts(fr,fropt,restrict_to_good_list=True)
                frsim = self.apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)
                fropt = self.apply_all_cuts(fropt,fropt,restrict_to_good_list=True)
                froptsim = self.apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)
                
                #hist.main(datafitreslist[i],glob.glob('%s/*/FITOPT000.FITRES'%simfitreslist[i])[0])
                hist.main(data=fr,sim=frsim)
                
    def mk_biascor_validplots(self):
        # make biascor plots
        plt.rcParams['figure.figsize'] = (16,5)
        plt.subplots_adjust(left=0.1,right=0.9,wspace=0.3)
        for scatmod in ['_G10','']:
            for simfitreslist,datafitreslist,opticalnirsimfitreslist,opticalnirdatafitreslist,label in \
                zip([self.nirsimfitreslist,self.opticalnirsimfitreslist],
                    [self.nirdatafitreslist,self.opticalnirdatafitreslist],
                    [self.opticalnirsimfitreslist,self.opticalnirsimfitreslist],
                    [self.opticalnirdatafitreslist,self.opticalnirdatafitreslist],
                    ['NIR','OpticalNIR']):
                if label == 'NIR': continue
                if scatmod == '_G10': continue
                
                plt.clf()
                if 'G10' not in scatmod:
                    ax1,ax2,ax3 = plt.subplot(131),plt.subplot(132),plt.subplot(133)
                    for ax in [ax1,ax2,ax3]: ax.set_xlabel('$z_{CMB}$',fontsize=15)
                    ax1.set_ylabel('$\mu-\mu_{sim}$')
                    ax2.set_ylabel('$s-s_{sim}$')
                    ax3.set_ylabel('$AV-AV_{sim}$')
                else:
                    ax1 = plt.subplot(111)
                for i,survey in enumerate(['CSP','PS1','DES']):
                    if 'G10' in scatmod and 'fit_nir' in simfitreslist[i]: continue
                    #fr = self.apply_all_cuts(fr,fropt,restrict_to_good_list=True)
                    print(f'{simfitreslist[i]}{scatmod}/*/FITOPT000.FITRES')

                    frsim = txtobj(glob.glob(f'{simfitreslist[i]}{scatmod}/*/FITOPT000.FITRES')[0],
                                   fitresheader=True)
                    froptsim = txtobj(glob.glob(f'{opticalnirsimfitreslist[i]}{scatmod}/*/FITOPT000.FITRES')[0],
                                      fitresheader=True)

                    if 'G10' in scatmod:
                        frsim = self.apply_salt2_cuts(frsim)
                        frsim.DLMAG = frsim.mB+0.14*frsim.x1 - 3.1*frsim.c
                        frsim.SIM_DLMAG = frsim.SIM_mB+0.14*frsim.SIM_x1 - 3.1*frsim.SIM_c
                    else:
                        frsim = self.apply_all_cuts(frsim,froptsim,restrict_to_good_list=False)

                    #fropt = self.apply_all_cuts(fropt,fropt,restrict_to_good_list=True)
                    #froptsim = self.apply_all_cuts(froptsim,froptsim,restrict_to_good_list=False)


                    if np.max(frsim.zHD) < 0.15: zbins = np.linspace(0.01,0.1,10)
                    else: zbins = np.linspace(0.1,0.8,10)

                    mu_bin = binned_statistic(
                        frsim.zCMB,frsim.DLMAG-frsim.SIM_DLMAG,
                        statistic='median',bins=zbins).statistic
                    mu_errbin = binned_statistic(
                        frsim.zCMB,frsim.DLMAG-frsim.SIM_DLMAG,
                        statistic=errfnc,bins=zbins).statistic
                    ax1.errorbar((zbins[1:]+zbins[:-1])/2.,mu_bin,yerr=mu_errbin,
                                 fmt='o-',color=f'C{i}',capsize=0,lw=2,label=survey)
                    if 'G10' not in scatmod:
                        st_bin = binned_statistic(
                            frsim.zCMB,frsim.STRETCH-frsim.SIM_STRETCH,
                            statistic='median',bins=zbins).statistic
                        st_errbin = binned_statistic(
                            frsim.zCMB,frsim.STRETCH-frsim.SIM_STRETCH,
                            statistic=errfnc,bins=zbins).statistic
                        av_bin = binned_statistic(
                            frsim.zCMB,frsim.AV-frsim.SIM_AV,
                            statistic='median',bins=zbins).statistic
                        av_errbin = binned_statistic(
                            frsim.zCMB,frsim.AV-frsim.SIM_AV,
                            statistic=errfnc,bins=zbins).statistic

                        
                        ax2.errorbar((zbins[1:]+zbins[:-1])/2.,st_bin,yerr=st_errbin,
                                     fmt='o-',color=f'C{i}',capsize=0,lw=2,label=survey)
                        ax3.errorbar((zbins[1:]+zbins[:-1])/2.,av_bin,yerr=av_errbin,
                                     fmt='o-',color=f'C{i}',capsize=0,lw=2,label=survey)

                    if label == 'OpticalNIR' and survey == 'CSP': #if label =='OpticalNIR' and survey == 'CSP':
                        import pdb; pdb.set_trace()
                ax1.legend()
                ax1.set_xlabel(r'$z_{CMB}$',fontsize=15)
                ax1.set_ylabel(r'$\mu - \mu_{\mathrm{sim}}$',fontsize=15)
                plt.savefig(f'figs/biascor_{label}{scatmod}.png')


        def run_cosmomc(self):

            raisin_root = 'RAISIN_STAT'
            raisin_dataset = 'cosmomc_data/RAISIN_combined_stat.cosmomc.txt'

            cosmomc_batch = """
#!/bin/bash
#SBATCH --time=35:30:00
#SBATCH --partition=sandyb
#SBATCH --account=pi-rkessler
#SBATCH --job-name=RAISIN_stat
#SBATCH --output=logs/RAISIN_stat
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

timeout 18000 mpirun -np=4 /project/kicp/viniciusvmb/wfirst/cosmomc_wfirst_gaussian/cosmomc_rebekah cosmoini/BBC_FoundPS1_spec1d_SURVCAL.ini
"""
            
            cosmomc_ini = f"""DEFAULT(batch2/JLA.ini)
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
file_root={raisin_root}
jla_dataset={raisin_dataset}"""

            with open('RAISIN_stat.batch','w') as batch_fout:
                print(cosmomc_batch,file=batch_fout)
            with open('RAISIN_stat.ini','w') as batch_ini:
                print(cosmomc_ini,file=batch_ini)
            
class sim:
    def __init__(self):
        pass

    def runsim(self):
        # LOW-Z
        os.system('sim_SNmix.pl $RAISIN_ROOT/cosmo/sim/inputs/lowz/SIMGEN_MASTER_LOWZ.INPUT')

        # PS1
        os.system('sim_SNmix.pl $RAISIN_ROOT/cosmo/sim/inputs/PS1/SIMGEN_MASTER_PS1SPEC.INPUT')

        # DES
        os.system('sim_SNmix.pl $RAISIN_ROOT/cosmo/sim/inputs/DES/SIMGEN_MASTER_DESSPEC.INPUT')
        
    def runfit_sim(self):
        # optical
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CSP_RAISIN_optnir.nml')
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CfA_RAISIN_optnir.nml')
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/PS1_RAISIN_optnir.nml')
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/DES_RAISIN_optnir.nml')
        
        # NIR
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CSP_RAISIN_NIR.nml')
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/CfA_RAISIN_NIR.nml')
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/PS1_RAISIN_NIR.nml')
        os.system('split_and_fit.pl $RAISIN_ROOT/cosmo/fit/sim/DES_RAISIN_NIR.nml')

class lcfit:
    def __init__(self):
        self.CSP_optical_fitresfile = 'output/fit_nir/CSP_RAISIN_optical.FITRES.TEXT'
        self.CfA_optical_fitresfile = 'output/fit_nir/CfA_RAISIN_optical.FITRES.TEXT'
        self.CSP_optnir_fitresfile = os.path.expandvars('$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT')
        
        self.nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/CfA_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                                  '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
        self.nirdatafitreslist = [os.path.expandvars(filepath) for filepath in self.nirdatafitreslist]

        self.opticaldatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optical.FITRES.TEXT',
                                      '$RAISIN_ROOT/cosmo/output/fit_optical/CfA_RAISIN_optical.FITRES.TEXT',
                                      '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT',
                                      '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optical.FITRES.TEXT']
        self.opticaldatafitreslist = [os.path.expandvars(filepath) for filepath in self.opticaldatafitreslist]

        
    def runfit_data(self):
        # optical
        os.system('snlc_fit.exe fit/CSP_RAISIN_optnir.nml')
        os.system('snlc_fit.exe fit/CfA_RAISIN_optnir.nml')
        os.system('snlc_fit.exe fit/PS1_RAISIN_optnir.nml')
        os.system('snlc_fit.exe fit/DES_RAISIN_optnir.nml')

        # nir
        os.system('snlc_fit.exe fit/CSP_RAISIN.nml')
        os.system('snlc_fit.exe fit/CfA_RAISIN.nml')
        os.system('snlc_fit.exe fit/PS1_RAISIN.nml')
        os.system('snlc_fit.exe fit/DES_RAISIN.nml')

    def plot_data(self):
        os.system('python raisin_cosmo/plot_snana.py -v CSPDR3 -o figs -a -20 -f fit/plot/CSP_RAISIN_optnir.nml --plotAll -F output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT')
        os.system('python raisin_cosmo/plot_snana.py -v CfAIR2 -o figs -a -20 -f fit/plot/CfA_RAISIN_optnir.nml --plotAll -F output/fit_optical/CfA_RAISIN_optnir.FITRES.TEXT')
        os.system('python raisin_cosmo/plot_snana.py -v PS1_RAISIN -o figs -a -20 -f fit/plot/PS1_RAISIN_optnir.nml --plotAll -F output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT')
        os.system('python raisin_cosmo/plot_snana.py -v DES_RAISIN -o figs -a -20 -f fit/plot/DES_RAISIN_optnir.nml --plotAll -F output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT')
        
    def add_pkmjd(self):
        CSPDR3_files = glob.glob('data/Photometry/CSPDR3_RAISIN/*.DAT')
        self.csp_opt_fr = txtobj(self.CSP_optnir_fitresfile,fitresheader=True)

        for c in CSPDR3_files:
            if 'PKMJD' in c: continue
            with open(c) as fin:
                with open(c.replace('.DAT','.PKMJD.DAT'),'w') as fout:
                    for line in fin:
                        if line.startswith('PEAKMJD'):
                            try:
                                print('PEAKMJD:  %.2f                # from optical SNooPy fit'%(
                                    self.csp_opt_fr.PKMJD[self.csp_opt_fr.CID == SNID][0]),file=fout)
                            except:
                                import pdb; pdb.set_trace()
                                sn = snana.SuperNova(c)
                                print(SNID,sn.REDSHIFT_FINAL)
                                print(line.replace('\n',''),file=fout)
                        elif line.startswith('SNID'):
                            SNID = line.split()[1].replace('\n','')
                            print(line.replace('\n',''),file=fout)
                        else:
                            print(line.replace('\n',''),file=fout)

        # need to do the same for CfA


    def hubbleplot():
        pass
        
def errfnc(x):
    return(np.std(x)/np.sqrt(len(x)))

def main(frsim_ps1='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
         frsim_ps1_oir='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR/FITOPT000.FITRES',
         #frsim_ps1='output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
         #output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
         frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

    plt.clf()
    ax = plt.axes() #subplot(121)
    
    frp = txtobj(frsim_ps1,fitresheader=True)
    frpo = txtobj(frsim_ps1_oir,fitresheader=True)
    frd = txtobj(frsim_des,fitresheader=True)

    zbins = np.linspace(0.01,0.7,20)
    #iSNR = frp.SNRMAX1 > 50
    #ps1_bin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic='median',bins=zbins).statistic
    #ps1_errbin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic=errfnc,bins=zbins).statistic
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ps1_bin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic='median',bins=zbins).statistic
    ps1_bin_oir = binned_statistic(
        frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic='median',bins=zbins).statistic

    #ps1_bin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic
    ps1_errbin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic=errfnc,bins=zbins).statistic
    ps1_errbin_oir = binned_statistic(
        frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic=errfnc,bins=zbins).statistic

    #ps1_errbin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic=errfnc,bins=zbins).statistic

    
    #des_bin = binned_statistic(
    #   frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic

    #ax2 = plt.subplot(122)
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
                fmt='o-',color='r',capsize=0,lw=2)
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='b',capsize=0,lw=2)

    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ax.set_xlabel('$z_{CMB}$',fontsize=15)
    ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
    import pdb; pdb.set_trace()
    
    return

def smearmodelcomp(
        frsim_ps1='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
        frsim_ps1_oir='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR_NEWSIMLIB_z/FITOPT000.FITRES',
        #frsim_ps1='output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
        #output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
        frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

    plt.clf()
    ax = plt.axes() #subplot(121)
    
    frp = txtobj(frsim_ps1,fitresheader=True)
    frpo = txtobj(frsim_ps1_oir,fitresheader=True)
    frd = txtobj(frsim_des,fitresheader=True)

    zbins = np.linspace(0.01,0.7,20)
    #iSNR = frp.SNRMAX1 > 50
    #ps1_bin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic='median',bins=zbins).statistic
    #ps1_errbin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic=errfnc,bins=zbins).statistic
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ps1_bin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic='median',bins=zbins).statistic
    ps1_bin_oir = binned_statistic(
        frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic='median',bins=zbins).statistic

    #ps1_bin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic
    ps1_errbin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic=errfnc,bins=zbins).statistic
    ps1_errbin_oir = binned_statistic(
        frpo.zCMB,(frpo.mB-frpo.SIM_mB)+0.14*(frpo.x1-frpo.SIM_x1) - 3.1*(frpo.c-frpo.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic=errfnc,bins=zbins).statistic
    
    
    #ps1_errbin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic=errfnc,bins=zbins).statistic

    
    #des_bin = binned_statistic(
    #   frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic

    #ax2 = plt.subplot(122)
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
                fmt='o-',color='r',capsize=0,lw=2,label='G10 model')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='b',capsize=0,lw=2,label='SNooPy OIR model')
    ax.legend()
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ax.set_xlabel('$z_{CMB}$',fontsize=15)
    ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
    import pdb; pdb.set_trace()
    
    return

def saltvsnoopy(
        frsim_ps1='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR_newsimlib_z/FITOPT000.FITRES',
        frsim_ps1_oir='output/sim/PS1_RAISIN_SIM/PS1_RAISIN_SNOOPY_OIR_real_nironly/FITOPT000.FITRES',
        #frsim_ps1='output/sim/PS1_RAISIN_SIM_NIR.FITRES.TEXT',
        #output/sim/PS1_RAISIN_SIM_NIR/PS1_RAISIN_SNOOPY/FITOPT000.FITRES',
        frsim_des='output/sim/DES_RAISIN_SIM.FITRES.TEXT'): #'output/sim/DES_RAISIN_SIM_NIR.FITRES.TEXT'):

    plt.clf()
    ax = plt.axes() #subplot(121)
    
    frp = txtobj(frsim_ps1,fitresheader=True)
    frpo = txtobj(frsim_ps1_oir,fitresheader=True)
    frd = txtobj(frsim_des,fitresheader=True)
    #igood = frpo.AVERR < 0.05
    #for k in frpo.__dict__.keys():
    #   frpo.__dict__[k] = frpo.__dict__[k][igood]
    
    zbins = np.linspace(0.01,0.7,20)
    #iSNR = frp.SNRMAX1 > 50
    #ps1_bin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic='median',bins=zbins).statistic
    #ps1_errbin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic=errfnc,bins=zbins).statistic
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ps1_bin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic='median',bins=zbins).statistic
    ps1_bin_oir = binned_statistic(
        frpo.zCMB,frpo.DLMAG-frpo.SIM_DLMAG,
        statistic='median',bins=zbins).statistic
    ps1_stretch_oir = binned_statistic(
        frpo.zCMB,frpo.STRETCH-frpo.SIM_STRETCH,
        statistic='median',bins=zbins).statistic
    ps1_av_oir = binned_statistic(
        frpo.zCMB,frpo.AV-frpo.SIM_AV,
        statistic='median',bins=zbins).statistic

    
    #ps1_bin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic
    ps1_errbin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic=errfnc,bins=zbins).statistic
    ps1_errbin_oir = binned_statistic(frpo.zCMB,
        frpo.DLMAG - frpo.SIM_DLMAG,
        statistic=errfnc,bins=zbins).statistic
    
    
    #ps1_errbin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic=errfnc,bins=zbins).statistic

    
    #des_bin = binned_statistic(
    #   frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic

    #ax2 = plt.subplot(122)
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
                fmt='o-',color='r',capsize=0,lw=2,label='G10 model')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='b',capsize=0,lw=2,label='SNooPy OIR model')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_stretch_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='C0',capsize=0,lw=2,label='SNooPy OIR model - $s$')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_av_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='C1',capsize=0,lw=2,label='SNooPy OIR model - A_V')

    ax.legend()
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ax.set_xlabel('$z_{CMB}$',fontsize=15)
    ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
    import pdb; pdb.set_trace()
    
    return

def saltvsnoopy_des(
        frsim_des='output/sim/DES_RAISIN_SIM/DES_RAISIN_SALT2/FITOPT000.FITRES',
        frsim_des_oir='output/sim/DES_RAISIN_SIM_NIR/DES_RAISIN_SNOOPY_OIR/FITOPT000.FITRES'):
    plt.clf()
    ax = plt.axes() #subplot(121)
    
    frp = txtobj(frsim_des,fitresheader=True)
    frpo = txtobj(frsim_des_oir,fitresheader=True)
    frd = txtobj(frsim_des,fitresheader=True)
    #igood = frpo.AVERR < 0.05
    #for k in frpo.__dict__.keys():
    #   frpo.__dict__[k] = frpo.__dict__[k][igood]
    
    zbins = np.linspace(0.01,0.7,20)
    #iSNR = frp.SNRMAX1 > 50
    #ps1_bin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic='median',bins=zbins).statistic
    #ps1_errbin = binned_statistic(
    #   frp.zCMB[iSNR],frp.DLMAG[iSNR]-frp.SIM_DLMAG[iSNR],
    #   statistic=errfnc,bins=zbins).statistic
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ps1_bin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic='median',bins=zbins).statistic
    ps1_bin_oir = binned_statistic(
        frpo.zCMB,frpo.DLMAG-frpo.SIM_DLMAG,
        statistic='median',bins=zbins).statistic
    ps1_stretch_oir = binned_statistic(
        frpo.zCMB,frpo.STRETCH-frpo.SIM_STRETCH,
        statistic='median',bins=zbins).statistic
    ps1_av_oir = binned_statistic(
        frpo.zCMB,frpo.AV-frpo.SIM_AV,
        statistic='median',bins=zbins).statistic

    
    #ps1_bin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic
    ps1_errbin = binned_statistic(
        frp.zCMB,(frp.mB-frp.SIM_mB)+0.14*(frp.x1-frp.SIM_x1) - 3.1*(frp.c-frp.SIM_c), #frp.DLMAG-(frp.SIM_mB + 0.14*frp.SIM_x1 - 3.1*frp.SIM_c + 19.3),
        statistic=errfnc,bins=zbins).statistic
    ps1_errbin_oir = binned_statistic(frpo.zCMB,
        frpo.DLMAG - frpo.SIM_DLMAG,
        statistic=errfnc,bins=zbins).statistic
    
    
    #ps1_errbin = binned_statistic(
    #   frp.zCMB,frp.DLMAG-frp.SIM_DLMAG,
    #   statistic=errfnc,bins=zbins).statistic

    
    #des_bin = binned_statistic(
    #   frd.zCMB,frd.DLMAG-frd.SIM_DLMAG,
    #   statistic='median',bins=zbins).statistic

    #ax2 = plt.subplot(122)
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin,yerr=ps1_errbin,
                fmt='o-',color='r',capsize=0,lw=2,label='G10 model')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_bin_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='b',capsize=0,lw=2,label='SNooPy OIR model')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_stretch_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='C0',capsize=0,lw=2,label='SNooPy OIR model - $s$')
    ax.errorbar((zbins[1:]+zbins[:-1])/2.,ps1_av_oir,yerr=ps1_errbin_oir,
                fmt='o-',color='C1',capsize=0,lw=2,label='SNooPy OIR model - A_V')

    ax.legend()
    #ax.errorbar((zbins[1:]+zbins[:-1])/2.,des_bin,
    #           fmt='o-',color='b',capsize=0,lw=2)
    ax.set_xlabel('$z_{CMB}$',fontsize=15)
    ax.set_ylabel('$\mu - \mu_{sim}$',fontsize=15)
    import pdb; pdb.set_trace()
    
    return

if __name__ == "__main__":
    #main()
    #smearmodelcomp()
    #saltvsnoopy()
    #saltvsnoopy_des()

    #sm = sim()
    #sm.runsim()
    #sm.runfit()
    
    #lcf = lcfit()
    #lcf.add_pkmjd()
    
    bc = biascor()
    bc.mk_sim_validplots()
    bc.mk_biascor_validplots()
    
    #bc.mkcuts()
    #bc.mk_biascor_files()
    #bc.apply_biascor()
    #bc.mk_cosmomc_files()
    #bc.submit_batch('/project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/cosmo/output/cosmo_fitres/RAISIN_combined_stat.cosmomc.txt',
    #               '/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_combined_stat.dataset',
    #               '/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_combined_stat.ini',
    #               '/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_combined_stat.job')
