import numpy as np
import matplotlib.pyplot as plt
import math
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.utils import data
from astropy.io import fits
import  aplpy
from astropy.utils.data import download_file
from astropy.io import fits
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
from astropy import log
log.setLevel('ERROR')
from pvextractor import PathFromCenter
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS

def Mass(NH2,image_header,X,d):
    '''X: abundace ratio between molecule and H2
    d: distance to the source in pc'''
    mu = 2.8 #Kauffmann et al. 2008
    mH = 1.674*10**-24 #g mass. of hydrogen
    au2cm = 1.49597871*10**13
    pixArea= abs(image_header['CDELT1'])*3600*abs(image_header['CDELT2'])*3600*(d*au2cm)**2
    #print(pixArea)
    #print(image_header)
    Msun = 1.9891*10**33
    #print(mu*mH*pixArea/Msun*(1000/image_header['CDELT3']))
    return mu*mH*pixArea/Msun*NH2/X*(1000/abs(image_header['CDELT3'])) #Msun /(kms^-1)

def generate_spec(image_file_co,d,vsyst,inc,vstart,vend,Name):
    '''image_file_co: CO file path (fits file)
    d:distance to the source in pc
    vsyst: system velocity in km/s
    inc: inclination angle in degrees (0 is pole on; 90 is face on)
    Name: save fits name (string)
    '''
    #CO
    image_data_co = fits.getdata(image_file_co, ext=0)
    image_header3 = fits.getheader(image_file_co, ext=0)
    bmaj = image_header3['BMAJ']*3600
    bmin = image_header3['BMIN']*3600
    velo = np.linspace(image_header3['CRVAL3']/1000.0,(image_header3['CRVAL3']+image_header3['CDELT3']*image_header3['NAXIS3'])/1000.0,image_header3['NAXIS3'])
    spec_co = np.nan_to_num(image_data_co)

    #print(velo)
    M_CO = Mass(spec_co,image_header3,10**-4,d)

    #plt.imshow(M_CO[243,:,:])
    #plt.show()
    #Crossing time estimates
    vout = abs(velo-vsyst)/np.cos(inc*np.pi/180)
    au2km = 149597871.
    dr = image_header3['CDELT2']*3600*d*au2km/np.sin(inc*np.pi/180) #1 pixel ring
    print('dr',dr)
    print('pixel size',image_header3['CDELT2'])
    #print('vout',vout)
    tau = dr/vout/60/60/24/365.25/10**3
     
    MSpect_final = np.full((M_CO.shape[0],M_CO.shape[1],M_CO.shape[2],len(vstart)),np.nan)
    for count in range(len(vstart)):
        for i in range(M_CO.shape[1]):
            for j in range(M_CO.shape[2]):
                MSpect = np.where(velo >vstart[count], M_CO[:,i,j], np.nan)
                MSpect_final[:,i,j,count] = np.where(velo > vend[count], np.nan, MSpect)/tau[:]*1000
    MSpect_final2 = np.nansum(MSpect_final,axis=3) 
    print(MSpect_final2.shape)
    dM_2d = np.nansum(MSpect_final2[:,:,:],axis=0)
    dM_2d[dM_2d ==0] = np.nan
    #write to fits
    header_d = fits.getheader(image_file_co)
    del header_d['NAXIS3']
    del header_d['PC3_1']
    del header_d['PC3_2']
    del header_d['PC3_3']
    del header_d['PC1_3']
    del header_d['PC2_3']
    del header_d['CTYPE3']
    del header_d['CRVAL3']
    del header_d['CDELT3']
    del header_d['CRPIX3']
    del header_d['CUNIT3']
    
    OUT = dM_2d[:,:] #'image0219_flux'  #input \n",
    hm = str(Name)+'_dM_map_2d.fits'  #output\n",   
    header_d['BUNIT']='Msun/Myr' #
    header_d['BTYPE']='Mass rate'
    fits.writeto(hm,OUT,header_d,overwrite=True)

    plt.imshow(dM_2d,origin='lower')
    plt.title(str(Name)+'_dM_map')
    plt.xlabel('RA [pixel]')
    plt.ylabel('DEC [pixel]')
    cbar = plt.colorbar()
    cbar.set_label('$M_\odot Myr$')
    plt.savefig(str(Name)+'_dM_map_2d.pdf',dpi=300)
    #plt.show()
    plt.close()
    #Momentum rate dP 
    dP_3dmap = MSpect_final2*vout.reshape(MSpect_final2.shape[0],1,1)
    dv = abs(velo[1]-velo[0])
    dP_2d = np.nansum(dP_3dmap[:,:,:]*dv,axis=0)
    dP_2d[dP_2d ==0] = np.nan

    OUT = dP_2d[:,:] #'image0219_flux'  #input \n",
    hm = str(Name)+'_dP_map_2d.fits'  #output\n",   
    
    header_d = fits.getheader(image_file_co) 
    header_d['BUNIT']='Msun km/s/Myr' #
    header_d['BTYPE']='Momentum rate'
    del header_d['NAXIS3']
    del header_d['PC3_1']
    del header_d['PC3_2']
    del header_d['PC3_3']
    del header_d['PC1_3']
    del header_d['PC2_3']
    del header_d['CTYPE3']
    del header_d['CRVAL3']
    del header_d['CDELT3']
    del header_d['CRPIX3']
    del header_d['CUNIT3']
    fits.writeto(hm,OUT,header_d,overwrite=True)
    plt.imshow(dP_2d,origin='lower')
    plt.title(str(Name)+'_dP_map')
    plt.xlabel('RA [pixel]')
    plt.ylabel('DEC [pixel]')
    cbar = plt.colorbar()
    cbar.set_label('$M_\odot km/s/Myr$')
    plt.savefig(str(Name)+'_dP_map_2d.pdf',dpi=300)
    plt.close()

    #dE Energy rate maps
    Msun = 1.989E30 #kg
    fac = Msun*(10**3)**2*10**7/1e43 # Msun (km/s)^2 to erg
    dE_3dmap = 1/2*MSpect_final2*vout.reshape(MSpect_final2.shape[0],1,1)**2*dv*fac
    dE_2d = np.nansum(dE_3dmap[:,:,:],axis=0)
    dE_2d[dE_2d ==0] = np.nan
    OUT = dE_2d[:,:] #'image0219_flux'  #input \n",
    hm = str(Name)+'_dE_map_2d.fits'  #output\n",

    header_d = fits.getheader(image_file_co)

    header_d['BUNIT']='1e43 erg/Myr' #
    header_d['BTYPE']='Energy rate'
    del header_d['NAXIS3']
    del header_d['PC3_1']
    del header_d['PC3_2']
    del header_d['PC3_3']
    del header_d['PC1_3']
    del header_d['PC2_3']
    del header_d['CTYPE3']
    del header_d['CRVAL3']
    del header_d['CDELT3']
    del header_d['CRPIX3']
    del header_d['CUNIT3']
    fits.writeto(hm,OUT,header_d,overwrite=True)

    plt.imshow(dE_2d,origin='lower')
    plt.title(str(Name)+'_dE_map')
    plt.xlabel('RA [pixel]')
    plt.ylabel('DEC [pixel]')
    cbar = plt.colorbar()
    cbar.set_label('$10^{43}erg/Myr$')
    plt.savefig(str(Name)+'_dE_map_2d.pdf',dpi=300)
    plt.close()

    #Mass map
    MSpect_final = np.full((M_CO.shape[0],M_CO.shape[1],M_CO.shape[2],len(vstart)),np.nan)
    for count in range(len(vstart)):
        for i in range(M_CO.shape[1]):
            for j in range(M_CO.shape[2]):
                MSpect = np.where(velo >vstart[count], M_CO[:,i,j], np.nan)
                MSpect_final[:,i,j,count] = np.where(velo > vend[count], np.nan, MSpect)
    MSpect_final2 = np.nansum(MSpect_final,axis=3)
    
    print(MSpect_final2.shape)
    M_2d = np.nansum(MSpect_final2[:,:,:],axis=0)
    M_2d[M_2d ==0] = np.nan
    #write to fits
    header_d = fits.getheader(image_file_co)
    del header_d['NAXIS3']
    del header_d['PC3_1']
    del header_d['PC3_2']
    del header_d['PC3_3']
    del header_d['PC1_3']
    del header_d['PC2_3']
    del header_d['CTYPE3']
    del header_d['CRVAL3']
    del header_d['CDELT3']
    del header_d['CRPIX3']
    del header_d['CUNIT3']

    OUT = M_2d[:,:] #'image0219_flux'  #input \n",
    hm = str(Name)+'_M_map_2d.fits'  #output\n",
    header_d['BUNIT']='Msun' #
    header_d['BTYPE']='Mass'
    fits.writeto(hm,OUT,header_d,overwrite=True)

    #Momentum Map
    P_3dmap = MSpect_final2*vout.reshape(MSpect_final2.shape[0],1,1)
    dv = abs(velo[1]-velo[0])
    P_2d = np.nansum(P_3dmap[:,:,:]*dv,axis=0)
    P_2d[P_2d ==0] = np.nan

    OUT = P_2d[:,:] #'image0219_flux'  #input \n",
    hm = str(Name)+'_P_map_2d.fits'  #output\n",   

    header_d = fits.getheader(image_file_co)
    header_d['BUNIT']='Msun km/s' #
    header_d['BTYPE']='Momentum rate'
    del header_d['NAXIS3']
    del header_d['PC3_1']
    del header_d['PC3_2']
    del header_d['PC3_3']
    del header_d['PC1_3']
    del header_d['PC2_3']
    del header_d['CTYPE3']
    del header_d['CRVAL3']
    del header_d['CDELT3']
    del header_d['CRPIX3']
    del header_d['CUNIT3']
    fits.writeto(hm,OUT,header_d,overwrite=True)
    plt.imshow(P_2d,origin='lower')
    plt.title(str(Name)+'_P_map')
    plt.xlabel('RA [pixel]')
    plt.ylabel('DEC [pixel]')
    cbar = plt.colorbar()
    cbar.set_label('$M_\odot km/s$')
    plt.savefig(str(Name)+'_P_map_2d.pdf',dpi=300)
    plt.close()

    #E Energy rate maps
    Msun = 1.989E30 #kg
    fac = Msun*(10**3)**2*10**7/1e43 # Msun (km/s)^2 to erg
    E_3dmap = 1/2*MSpect_final2*vout.reshape(MSpect_final2.shape[0],1,1)**2*dv*fac
    E_2d = np.nansum(E_3dmap[:,:,:],axis=0)
    E_2d[E_2d ==0] = np.nan
    OUT = E_2d[:,:] #'image0219_flux'  #input \n",
    hm = str(Name)+'_E_map_2d.fits'  #output\n",

    header_d = fits.getheader(image_file_co)

    header_d['BUNIT']='1e43 erg' #
    header_d['BTYPE']='Energy rate'
    del header_d['NAXIS3']
    del header_d['PC3_1']
    del header_d['PC3_2']
    del header_d['PC3_3']
    del header_d['PC1_3']
    del header_d['PC2_3']
    del header_d['CTYPE3']
    del header_d['CRVAL3']
    del header_d['CDELT3']
    del header_d['CRPIX3']
    del header_d['CUNIT3']
    fits.writeto(hm,OUT,header_d,overwrite=True)

    plt.imshow(E_2d,origin='lower')
    plt.title(str(Name)+'_E_map')
    plt.xlabel('RA [pixel]')
    plt.ylabel('DEC [pixel]')
    cbar = plt.colorbar()
    cbar.set_label('$10^{43}erg$')
    plt.savefig(str(Name)+'_E_map_2d.pdf',dpi=300)
    plt.close()


    return velo,tau, M_2d, P_2d, E_2d, dM_2d, dP_2d, dE_2d

#System velocity:
vsyst = [1.3] # Codella et al. 2014

d = 450 #pc Zinnecker et al. 1998; Codella et al. 2014
#inc = 57.3 this is for random

inc = [86.0] #4 degrees to the plane of sky (Claussen et al. 1998) 

path = './'
image_file_co = 'Den_HH212_CO32.fits'
Name= 'HH212'
vsyst_in =  vsyst[0]
inc_in = inc[0] # outflinw inclination angle

vstart = [-23.73,8.01] #in km/s
vend   = [-3.42,17.96]

generate_spec(image_file_co,d,vsyst_in,inc_in,vstart,vend,Name)
