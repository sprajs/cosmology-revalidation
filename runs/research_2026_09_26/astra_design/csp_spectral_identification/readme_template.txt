This is the readme file for the NIR spectral template of SNe Ia that are published in Lu et al. 2023, titled as 
"Carnegie Supernova Project-II: Near-infrared spectral diversity and template of Type Ia Supernovae".
The link to the publication: https://ui.adsabs.harvard.edu/abs/2022arXiv221105998L/abstract

*** Note the the wavelength in um ***

---------------------------------------------------------------------------------------------------
## File: NIR_Ia_template_buildingblocks.pkl ###
The file that contains the PCA+GPR that were performed on the 7 wavelength wavelength regions
You can access the template which is a function of phase and sBV, using the python package BYOST. 

First need to "pip install BYOST" in terminal. 
For more instuctions of how to use the package please check
see instructions in https://github.com/DeerWhale/BYOST. 
...............................................................................................
import BYOST
df_buildingblock = pd.read_pickle('NIR_Ia_template_buildingblocks.pkl')
phase = 10  # input phase of the template spectra
sBV = 0.9  # input sBV of the template spectra
template = BYOST.template.get_template(df_buildingblock,phase,sBV,return_template_error=True)
...............................................................................................
## set return_template_error=False if you don't need template_flux_err,it takes some time to compute
## template[0] is the template wavelength, and template[1] is the template flux
## template[2] is the template flux error if return_template_error=True 
---------------------------------------------------------------------------------------------------
