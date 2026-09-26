This is the readme file for the NIR spectra of SNe Ia that are published in Lu et al. 2023, titled as 
"Carnegie Supernova Project-II: Near-infrared spectral diversity and template of Type Ia Supernovae".
The link to the publication: https://ui.adsabs.harvard.edu/abs/2022arXiv221105998L/abstract

This folder contains 339 spectra of 98 individual SNe obtained as part of the Carnegie Supernova Project-II. 
These spectra were obtained with the FIRE spectrograph on the 6.5m Magellan Baade telescope and have a 
spectral range of 0.8–2.5um. 

*** Note the the wavelength in the files are all in um ***

---------------------------------------------------------------------------------------------------
## File: Lu2023_TableA1_extension.csv ###
List of SNe Ia and the NIR spectra published in this work, include the information that 
were listed in Table A1 plus few more columns. Here are the columns information:

- SNe information:
	name: SN name
	untargeted?: whether the SN were among "targeted" or "untargted" search
	zhel: heliocentric redshift of the SN
	Tmax(MJD), e(Tmax): the time of the rest frame B-band maximum in MJD, and its error
	sBV,e(sBV): The color stretch sBV fitted with SNooPy, and its error
	EBV_MW: The Milkyway extinction E(B-V)
	spec_counts: The number of spectra of this SN

- Spectra information:
	filename: the file name of the spectrum 
	Date: the UT data of when the spectrum was observed
	MJD: the MJD of when the spectrum was observed
	epoch: the rest-frame phase of the spectrum
	EXPTOT: on-target exposure time of the spectrum excluding overhead
	SNRY,SNRJ,SNRH,SNRH,SNRK: the median signal-to-noise ratio of the spectrum in Y,J,H,K band
	host_contamination: if the spectrum is contaminated by host galaxy

---------------------------------------------------------------------------------------------------
## Folder: observed_spectra ###
The observed NIR spectra in this work, the file name (without the '.txt' or '.fits extension) is 
matched with the "filename" column in fileLu2023_TableA1_extension.csv.

- For txt files, each one contains 3 columns:
	0: the wavelength in observer frame (in unit of um)
	1: flux
	2: flux error

- For fits files of FIRE:
	from astropy.io import fits
	hdu = fits.open(file)
        wv,fx,fx_err = hdu[0].data[0],hdu[0].data[1],hdu[0].data[2]

---------------------------------------------------------------------------------------------------