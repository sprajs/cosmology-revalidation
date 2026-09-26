#!/usr/bin/env python
# D. Jones - 11/21/19
"""
run:
DES/PS1 RAISIN simulations

check:
SIMLIB files vs. data
sim parameters vs. data
distance bias
"""

import os

def runsim():
	# PS1
	os.system('snlc_sim.exe sim/inputs/PS1/SIMGEN_PS1SPEC.INPUT')
	os.system('rm -r /Users/David/Dropbox/research/RAISIN/cosmo/output/sim/SIMDIR/PS1_RAISIN_SNOOPY')
	os.system('mv /usr/local/SNDATA_ROOT/SIM/PS1_RAISIN_SNOOPY /Users/David/Dropbox/research/RAISIN/cosmo/output/sim/SIMDIR')
	os.system('snlc_fit.exe fit/sim/PS1_RAISIN_SIM.nml')
	
	# DES
	os.system('snlc_sim.exe sim/inputs/DES/sim_DES_SNOOPY.input')
	os.system('rm -r /Users/David/Dropbox/research/RAISIN/cosmo/output/sim/SIMDIR/DES_RAISIN_SNOOPY')
	os.system('mv /usr/local/SNDATA_ROOT/SIM/DES_RAISIN_SNOOPY /Users/David/Dropbox/research/RAISIN/cosmo/output/sim/SIMDIR')

	# low-z
	os.system('snlc_sim.exe sim/inputs/lowz/sim_DES_SNOOPY.input')
	os.system('rm -r /Users/David/Dropbox/research/RAISIN/cosmo/output/sim/SIMDIR/DES_RAISIN_SNOOPY')
	os.system('mv /usr/local/SNDATA_ROOT/SIM/DES_RAISIN_SNOOPY /Users/David/Dropbox/research/RAISIN/cosmo/output/sim/SIMDIR')
	
def check_simlib():
	pass

def check_simpar():
	# PS1
	os.system('python code/ovdatamc.py output/fit_optical/PS1_RAISIN_optical.FITRES.TEXT output/sim/PS1_RAISIN_SIM.FITRES.TEXT MURES:SNRMAX1:AV:STRETCH --journal --cutwin MURES -2 2 --nbins 10 --clobber -o figs/sim_PS1.png')

	# DES
	os.system('python code/ovdatamc.py output/fit_optical/DES_RAISIN_optical.FITRES.TEXT output/sim/DES_RAISIN_SIM.FITRES.TEXT MURES:SNRMAX1:AV:STRETCH --journal --cutwin MURES -2 2 --nbins 10 --clobber -o figs/sim_DES.png')

	# low-z
	
def bias():
	pass

if __name__ == "__main__":
	runsim()

