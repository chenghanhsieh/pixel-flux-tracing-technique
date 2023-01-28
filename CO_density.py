import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import math
from astropy.io import fits


def Ray_Jeans(I,nu,bmin,bmaj):
    """This function change the flux (Jy) to brightness temperautre (K)
    using the Rayleigh Jeans limit
    I: flux in Jy/beam
    nu: frequency in Hz
    bmin: beam minor axis (arcsec) 
    bmax: beam major axis (arcsec)"""
    # See https://science.nrao.edu/facilities/vla/proposing/TBconv
    #h  = 6.62607004*10**-34 
    #kb = 1.38064852*10**-23  
    return 1.222*10**6*I/(nu/10**9)**2/bmaj/bmin #K/beam

def N(I,TEX,chanw): #Dunham et al. 2014 C.10.
    '''for J3-2
    I: Integrated Intensity map in K kms^-1
    TEX: Excitation temperature ~50 K
    XCO: H2 to CO abundance factor ~10^-4 #Frerking et al. 1982;Lacy et al. 1994; Hatchell et al. 2007a)
    chanw: channel width in m/s
    '''
    return 2.8488138*10**13*np.abs(chanw)/1000.0*I


Name = 'calibrated.ms.image.line.HH212.spw3.CO.rotated-22deg.fits'
path = './'
image_file = path + Name
image_data = fits.getdata(image_file, ext=0)
image_header = fits.getheader(image_file, ext=0)
nu = 345.7959899*10**9
bmaj = image_header['BMAJ']*3600
print(bmaj)
bmin = image_header['BMIN']*3600
I = image_data
TMB = Ray_Jeans(I,nu,bmin,bmaj)
TEX = 50 #K
chanw = image_header['CDELT3']

Den_co = N(TMB,TEX,chanw)

#Write file to fits
OUT = Den_co[:,:,:] #'image0219_flux'  #input \n",
hm = './Den_HH212_CO32.fits'  #output\n",

header_d = image_header #header imformation\n", #fits.getheader(im+'.fits')

header_d['BUNIT']='cm^-2' #
header_d['BTYPE']='CO column_density'

fits.writeto(hm,OUT,header_d,overwrite=True)
