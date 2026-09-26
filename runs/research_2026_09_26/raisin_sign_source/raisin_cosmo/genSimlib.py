#!/usr/bin/env python
# D. Jones - 7/26/19
"""generate simlib files for Pan-STARRS and DES RAISIN SNe"""

ps1lc = '/Users/David/Dropbox/research/raisin/SALT3_NIR/raisindata/snana_optical/PS*.snana.dat'
deslc = '/Users/David/Dropbox/research/raisin/SALT3_NIR/raisindata/snana_optical/DE*.snana.dat'
lowzlc = '/Users/David/Dropbox/research/raisin/SALT3_NIR/data/snana/sn*snana.dat'

import os
import numpy as np
import glob
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import snana
from txtobj import txtobj

ps1simlib = os.path.expandvars('$SNDATA_ROOT/simlib/PS1/Pantheon/PS1MD_FULL_fluxcorr.simlib')
dessimlib = os.path.expandvars('/Users/david/Dropbox/research/raisin/sim/DES_DIFFIMG.SIMLIB') #$SNDATA_ROOT/simlib/DES/DES_Y1reproc2.SIMLIB')
#lowzsimlib = os.path.expandvars('$SNDATA_ROOT/simlib/PS1/Pantheon/PS1_LOWZ_COMBINED.SIMLIB')
lowzsimlib = 'LOWZ_CSP.SIMLIB'
snid_mjd,t0_mjd = np.loadtxt('data/raisin_t0.txt',unpack=True,dtype=str,usecols=[0,1])
t0_mjd = t0_mjd.astype(float)
fr = txtobj('data/pantheon_plus_foundation.fitres',fitresheader=True)
lowz_snid_mjd,lowz_t0_mjd = fr.CID,fr.PKMJD
#lowz_snid_mjd,lowz_t0_mjd = np.loadtxt('data/lowz_t0.txt',unpack=True,dtype=str,usecols=[0,1])
#lowz_t0_mjd = lowz_t0_mjd.astype(float)

PS1header = """TELESCOPE: PS1 
SURVEY: PS1MD    FILTERS: grizJH
USER: djones     HOST: data-scope

              FILTs  LOG10(SNR)   ERROR/ERROR(SNANA)
FLUXERR_COR:  griz  -0.90      0.9810  0.8429  1.1278  1.0286
FLUXERR_COR:  griz  -0.70      0.9696  0.9162  1.0842  1.2051
FLUXERR_COR:  griz  -0.50      0.9928  0.9223  1.3334  1.1051
FLUXERR_COR:  griz  -0.30      1.0289  0.9763  1.1678  1.1607
FLUXERR_COR:  griz  -0.10      1.0136  0.9631  1.2165  1.1329
FLUXERR_COR:  griz   0.10      0.9858  0.8963  1.1406  1.1376
FLUXERR_COR:  griz   0.30      0.9530  0.9588  1.1088  1.1068
FLUXERR_COR:  griz   0.50      0.9728  0.9221  1.0763  1.0244
FLUXERR_COR:  griz   0.70      0.9366  0.8935  1.0134  0.9605
FLUXERR_COR:  griz   0.90      0.9825  0.8681  0.9181  0.9193
FLUXERR_COR:  griz   1.10      0.9150  0.8335  0.8808  0.9036
FLUXERR_COR:  griz   1.30      0.9322  0.8153  0.8303  0.9268
FLUXERR_COR:  griz   1.50      0.9272  0.8807  0.8610  0.7724
FLUXERR_COR:  griz   1.70      1.0148  0.8360  1.0080  0.7076
FLUXERR_COR:  griz   1.90      1.0414  0.8594  1.2145  1.1777
FLUXERR_COR:  griz   2.10      1.0000  1.0000  1.0000  1.0000
FLUXERR_COR:  griz   2.30      1.0000  1.0000  1.0000  1.0000
FLUXERR_COR:  griz   2.50      1.0000  1.0000  1.0000  1.0000
FLUXERR_COR:  griz   2.70      1.0000  1.0000  1.0000  1.0000
FLUXERR_COR:  griz   2.90      1.0000  1.0000  1.0000  1.0000

BEGIN LIBGEN  2016-3-28 
# -------------------------------------------- 
"""

DESheader = """SURVEY: DES      FILTERS: griz      TELESCOPE: CTIO
USER: rkessler    HOST: des20.fnal.gov
PIXSIZE: 0.27
PSF_UNIT: ARCSEC_FWHM

COMMENT: 'Based on 2nd reprocessing of Y1 DES (Fall 2013-Feb 2014)'
COMMENT: 'This simlib was used to generate SNIa for the red'
COMMENT: 'histograms in Figures 11-13 of arXiv:1507.05137
COMEMNT: 'http://adsabs.harvard.edu/abs/2015arXiv150705137K'

#--------------------------------------------------------------------------


BEGIN LIBGEN
"""

lowzheader = """SUBSURVEY_LIST: CFA3S,CFA3K,CFA4p1,CFA4p2,CSP,CFA1,CFA2
SURVEY:         PS1_LOWZ_COMBINED
FILTERS:        BVgriYJH
# abcdefghijklmnopqrstuvwxyzABCDEFGHIJK
PSF_UNIT:       ARCSEC_FWHM
SKYSIG_UNIT:    ADU_PER_SQARCSEC

# Assume instrument parameters for GROUND
# Assume SKYMAG(  2700.) = 23.80 mag/asec^2
# Assume SKYMAG(  3714.) = 22.70 mag/asec^2
# Assume SKYMAG(  4790.) = 21.00 mag/asec^2
# Assume SKYMAG(  6220.) = 20.40 mag/asec^2
# Assume SKYMAG(  7544.) = 19.40 mag/asec^2
# Assume SKYMAG(  8679.) = 18.10 mag/asec^2
# Assume SKYMAG( 10095.) = 17.90 mag/asec^2

BEGIN LIBGEN"""

hstdict = {'J':(0.137/0.25,26.2303),
		   'H':(0.151/0.25,25.9463)}
altnamesfinaldict = {'PScK450082':'PS1-12CBD',
					 'PScK450339':'PS1-12CGY',
					 'PScB480464':'PS1-13LO',
					 'PScB480794':'PS1-13QQ',
					 'PScA470041':'PS1-A470041',
					 'PScA470110':'PS1-A470110',
					 'PScA470240':'PS1-A470240',
					 'PScC490037':'PS1-C490037',
					 'PScC490521':'PS1-C490521',
					 'PScD500100':'PS1-D500100',
					 'PScD500301':'PS1-D500300',
					 'PScF510457':'PS1-E510457',
					 'PScF520062':'PS1-F520062',
					 'PScF520188':'PS1-F520188',
					 'PScH540087':'PS1-H540087',
					 'PScF520107':'PS1-F520107',
					 'PScG530251':'PS1-G530251',
					 'PScH540118':'PS1-H540118',
					 'PScJ550202':'PS1-J550202',
					 'PScJ560027':'PS1-J560027',
					 'PScJ560054':'PS1-J560054',
					 'PScJ440005':'SN01',
					 'PScJ440236':'SN02',
                     'DES15C1nhv':'DES15C1NHV',
                     'DES15C3odz':'DES15C3ODZ',
                     'DES15E2mhy':'DES15E2MHY',
                     'DES15E2nlz':'DES15E2NLZ',
                     'DES15E2uc':'DES15E2UC',
                     'DES15X2kvt':'DES15X2KVT',
                     'DES15X2mey':'DES15X2MEY',
                     'DES15X2nkz':'DES15X2NKZ',
                     'DES16C1cim':'DES16C1CIM',
                     'DES16C2cva':'DES16C2CVA',
                     'DES16C3cmy':'DES16C3CMY',
                     'DES16E1dcx':'DES16E1DCX',
                     'DES16E2clk':'DES16E2CLK',
                     'DES16E2cqq':'DES16E2CQQ',
                     'DES16E2cxw':'DES16E2CXW',
                     'DES16E2rd':'DES16E2RD',
                     'DES16S1agd':'DES16S1AGD',
                     'DES16S1bno':'DES16S1BNO',
                     'DES16S2afz':'DES16S2AFZ',
                     'DES16X1cpf':'DES16X1CPF',
                     'DES16X2crr':'DES16X2CRR',
                     'DES16X3cry':'DES16X3CRY',
                     'DES16X3zd':'DES16X3ZD',
                     'SNABELL370':'SNABELL370'}
filtdict = {'u':'B','v':'V','y':'g','z':'r','A':'i','Y':'Y','J':'J','H':'H',
			'B':'B','V':'V','g':'g','r':'r','i':'i','t':'u','w':'w','x':'x','l':'l','f':'f'}

def mkps1simlib():

	with open(ps1simlib) as fin:
		lines = fin.readlines()
	fout = open('PS1_RAISIN.simlib','w')
	print(PS1header,file=fout)
	lcfiles = glob.glob(ps1lc)
	for lci,lcf in enumerate(lcfiles):
		sn = snana.SuperNova(lcf)

		libidtext = """LIBID: %i
RA: 0.0    DECL: 0.0   NOBS: %i
MWEBV: %.3f   PIXSIZE: 0.250
REDSHIFT: %.3f   PEAKMJD: %.1f
FIELD: %s"""%(lci,len(sn.FLT),float(sn.MWEBV.split()[0]),sn.z,t0_mjd[sn.SNID == snid_mjd][0],sn.SNID)
		print(libidtext,file=fout)
		
		for m,f,i in zip(sn.MJD,sn.FLT,range(len(sn.FLT))):
			if f == 'H' or f == 'J':
				image = glob.glob('/Users/David/research2/hstdata/*%s*.e01/*_sub_masked.fits'%altnamesfinaldict[sn.SNID])[0]
				data = fits.getdata(image)
				header = fits.getheader(image)
				mn,md,skysig = sigma_clipped_stats(data)
				skysig *= np.sqrt(0.13/0.25)*header['EXPTIME']
				skysig = np.sqrt(skysig**2. + 0.01**2.)
				
				simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
					m,i,f,skysig,hstdict[f][0],hstdict[f][1]+2.5*np.log10(header['EXPTIME']),0.0)
				print(simlibline,file=fout)
				#import pdb; pdb.set_trace()
				continue
				
			for line in lines:
				if not line.startswith('S:'): continue
				mjd = float(line.split()[1])
				if np.abs(mjd - m) < 0.01:
					skysig = float(line.split()[6])
					psf = float(line.split()[7])
					zpt = float(line.split()[10])
					zpterr = float(line.split()[11])
					
					simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
						m,i,f,skysig,psf,zpt,zpterr)
					print(simlibline,file=fout)
					break

		print("""END_LIBID: %i
    
# --------------------------------------------"""%lci,file=fout)
	print('END_OF_SIMLIB:',file=fout)

def mkdessimlib():

	with open(dessimlib) as fin:
		lines = fin.readlines()
	fout = open('DES_RAISIN.simlib','w')
	print(DESheader,file=fout)
	lcfiles = glob.glob(deslc)
	for lci,lcf in enumerate(lcfiles):
		sn = snana.SuperNova(lcf)

		libidtext = """# --------------------------------------------
LIBID: %i
RA: 54.520306    DECL: -29.391666   NOBS: %i
MWEBV: %.3f   PIXSIZE: 0.270
REDSHIFT: %.3f   PEAKMJD: %.1f
FIELD: %s  # CCDS: [42]"""%(
	lci,len(sn.FLT),float(sn.MWEBV.split()[0]),sn.z,t0_mjd[sn.SNID == snid_mjd][0],sn.SNID)
		print(libidtext,file=fout)
		
		for m,f,i in zip(sn.MJD,sn.FLT,range(len(sn.FLT))):
			if f == 'H' or f == 'J':
				image = glob.glob('/Users/david/research2/raisin2/*%s*.e01/*_sub_masked.fits'%altnamesfinaldict[sn.SNID])[0]
				data = fits.getdata(image)
				header = fits.getheader(image)
				mn,md,skysig = sigma_clipped_stats(data)
				skysig *= np.sqrt(0.13/0.25)*header['EXPTIME']
				skysig = np.sqrt(skysig**2. + 0.01**2.)
				
				simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
					m,i,f,skysig,hstdict[f][0],hstdict[f][1]+2.5*np.log10(header['EXPTIME']),0.0)
				print(simlibline,file=fout)
				continue
				
			for line in lines:
				if not line.startswith('S:'): continue
				mjd = float(line.split()[1])
				if np.abs(mjd - m) < 0.01:
					skysig = float(line.split()[6])
					psf = float(line.split()[7])
					zpt = float(line.split()[10])
					zpterr = float(line.split()[11])
					
					simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
						m,i,f,skysig,psf,zpt,zpterr)
					print(simlibline,file=fout)
					break

		print("""END_LIBID: %i
    
# --------------------------------------------"""%lci,file=fout)
	print('END_OF_SIMLIB:',file=fout)

def mklowzsimlib():

	with open(lowzsimlib) as fin:
		lines = fin.readlines()
	fout = open('LOWZ_RAISIN.simlib','w')
	print(lowzheader,file=fout)
	lcfiles = glob.glob(lowzlc)
	count = 0
	for lci,lcf in enumerate(lcfiles):
		sn = snana.SuperNova(lcf)
		if sn.SNID[2:] not in lowz_snid_mjd:
			print('SNID %s will not be used!'%sn.SNID)
			continue
		count += 1

		simliblines = []
		nobs = 0
		for m,f,i in zip(sn.MJD,sn.FLT,range(len(sn.FLT))):
			#if f not in 'uvyzAYJHBVgri':
			#	print(f)
			#	continue
			if f in 'YJH': #['J','H','Y','a','b','c','d','e','f','g','l','m','n']:
				skysig = 5
				zpt = 25
				
				simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
					m,i,filtdict[f],skysig,1.0,zpt,0.0)
				simliblines += [simlibline]
				nobs += 1
				#print(simlibline,file=fout)
				continue

			for line in lines:
				
				if not line.startswith('S:') and not line.startswith('SUBSURVEY'): continue
				elif line.startswith('SUBSURVEY'):
					subsurvey = line.split()[1].replace(' ','')
				else:
					mjd = float(line.split()[1])
					if np.abs(mjd - m) < 0.01:
						#filt = line.split()[3]
						#if filt not in 'BVgri': continue
						skysig = np.sqrt(float(line.split()[6])**2.+0.01**2.)
						psf = float(line.split()[7])
						zpt = float(line.split()[10])
						zpterr = float(line.split()[11])

						if f in 'uvwxyZAYJH':
							simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
								m,i,filtdict[f],skysig,psf,zpt,zpterr)
							simliblines += [simlibline]
							nobs += 1
						break
				if subsurvey != 'CSP': continue

		libidtext = """# --------------------------------------------
LIBID:        %i
RA:   167.624320        DECL:   55.160840     NOBS:   %i
MWEBV:  %.3f  PIXSIZE:  0.500     REDSHIFT: %.5f     PEAKMJD: %.1f
FIELD: %s"""%(
	count,nobs,float(sn.MWEBV.split()[0]),sn.z,lowz_t0_mjd[sn.SNID[2:] == lowz_snid_mjd][0],sn.SNID)
		print(libidtext,file=fout)
				
		for simlibline in simliblines:
			print(simlibline,file=fout)

		print("""END_LIBID: %i
    
# --------------------------------------------"""%count,file=fout)
	print('END_OF_SIMLIB:',file=fout)

def modcspsimlib(filename='LOWZ_CSP.SIMLIB'):
	from random import sample
	lcfiles = glob.glob(lowzlc)
	
	with open(filename,'r') as fin:
		with open('LOWZ_CSP_MODFILT.SIMLIB','w') as fout:
			for line in fin:
				line = line.replace('\n','')
				print(line)
				if line.startswith('S:'):
					if has_sn and len(simliblines):
						filt = line.split()[3]
						lineparts = line.split()
						lineparts[3] = filtdict[filt]
						newline = ' '.join(lineparts)
						print(newline,file=fout)
				elif line.startswith('LIBID:') and 'cadence' in line:
					simlibsnid = line.split()[-1]
					print(line,file=fout)
				elif line.startswith('NOBS:'):
					lcfile = glob.glob('/Users/David/Dropbox/research/raisin/SALT3_NIR/data/snana/sn*%s*snana.dat'%simlibsnid)
					if not len(lcfile):
						has_sn = False
						continue
					else: lcfile = lcfile[0]
					has_sn = True
					sn = snana.SuperNova(lcfile)

					simliblines = []
					lineparts = line.split()
					nobs = int(line.split()[1])
					for m,f,i in zip(sn.MJD,sn.FLT,range(len(sn.FLT))):
						if f in 'YJH':
							skysig = 5
							zpt = 25
				
							simlibline = "S: %.3f     %i %s  1.0  0.00    %.2f  %.2f 0.00 0.000  %.3f  %.3f  99.000"%(
								m,i,filtdict[f],skysig,1.0,zpt,0.0)
							simliblines += [simlibline]
							nobs += 1
					lineparts[1] = str(nobs)
					print(' '.join(lineparts),file=fout)
				elif line.startswith('#  MJD'):
					print(line,file=fout)
					if has_sn:
						for s in simliblines:
							print(s,file=fout)
				else:
					print(line,file=fout)

if __name__ == "__main__":
	mkdessimlib()
	#mkps1simlib()
	#mklowzsimlib()
	#modcspsimlib()
