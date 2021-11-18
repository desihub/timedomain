from abc import ABC
import numpy as np
import numpy.ma as ma

import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.constants import c
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

from astropy.table import Table

import os
from glob import glob

from astropy.table import Table, join, vstack, hstack, unique

from desispec.io import read_spectra
from desispec.spectra import Spectra
from desispec.coaddition import coadd_cameras

import redrock.templates

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings

from scipy.signal import butter,filtfilt

from desispec.interpolation import resample_flux
from desispec.resolution import Resolution

####

def singlegaus(x,A,mu,sigma):
    y = A*np.exp(-(x-mu)**2/(2*sigma**2))
    return y
def doublegaus(x,A1,mu1,sigma1,A2,mu2,sigma2):
    y = A1*np.exp(-(x-mu1)**2/(2*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2*sigma2**2)) 
    return y


def triplegaus(x,A1,mu1,sigma1,A2,mu2,sigma2,A3,mu3,sigma3):
    y = A1*np.exp(-(x-mu1)**2/(2*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2*sigma2**2)) + \
        A3*np.exp(-(x-mu3)**2/(2*sigma3**2))
    return y

def Combine_multifilt(wave,flux, mask,ivar):
        difwave_single = []
        difflux_single = []
        difmask_single = []
        difivar_single = []
        for band in wave:
            difwave_single += list(wave[band]) 
            difflux_single += list(flux[band])
            difivar_single += list(ivar[band])
            difmask_single += list(mask[band])
                
        diftable = Table([difwave_single, difflux_single,difmask_single,difivar_single], names = ('wave','flux', 'mask','ivar'))

        diftable.sort('wave')
        difwave_single = np.array(diftable['wave'])
        difflux_single = np.array(diftable['flux'])
        difivar_single = np.array(diftable['ivar'])
        difmask_single = np.array(diftable['mask'])
        return difwave_single,difflux_single,difmask_single,difivar_single
       

def line_finder(wave, flux,ivar,mask,z):
    c = 2.99e5 #km/s

    lines = ['Halpha','Hbeta', 'Hgamma','HeII4686','OIII','NIII']
    restwavelengths = [6562,4861,4340,4686,5007,4100]
    
    wave,flux,mask,ivar = Combine_multifilt(wave,flux,mask,ivar)
    wave = wave/(1+z)

    HB_center = list(abs(wave-4861.4)).index(min(abs(wave - 4861.4)))
    HBmask = mask[HB_center - 500:HB_center + 500]
    HBroi = ma.masked_array(wave[HB_center - 500:HB_center + 500], HBmask)
    HBflux = ma.masked_array(flux[HB_center - 500:HB_center + 500], HBmask)
    HBsigma = np.sqrt(1/ivar[HB_center - 500:HB_center + 500])
    
    Ha_center = list(abs(wave-6562.79)).index(min(abs(wave - 6562.79)))
    Hamask = mask[Ha_center - 300:Ha_center + 300]
    Haroi = ma.masked_array(wave[Ha_center - 300:Ha_center + 300], Hamask)
    Haflux = ma.masked_array(flux[Ha_center - 300:Ha_center + 300], Hamask)
    Hasigma = np.sqrt(1/ivar[Ha_center - 300:Ha_center + 300])
    
    NIII_center = list(abs(wave-4200)).index(min(abs(wave - 4200)))
    NIIImask = mask[NIII_center - 200:NIII_center + 200]
    NIIIroi = ma.masked_array(wave[NIII_center- 200:NIII_center+200], NIIImask)
    NIIIflux = ma.masked_array(flux[NIII_center- 200:NIII_center+200], NIIImask)    
    NIIIsigma = np.sqrt(1/ivar[NIII_center - 200:NIII_center + 200])
    
    
    HBopt, HBcov = curve_fit(triplegaus, HBroi[HBroi.mask], HBflux[HBflux.mask], \
                             p0 = [1,4686,3,1,4861,5,5,5007,0.125],sigma = HBsigma, bounds = (0,np.inf))
    
    Haopt, Hacov = curve_fit(singlegaus, Haroi[Haroi.mask], Haflux[Haflux.mask], \
                         p0 = [3,6562,5],sigma = Hasigma, bounds = (0,np.inf))
   
    NIIIopt, NIIIcov = curve_fit(doublegaus, NIIIroi[NIIIroi.mask], NIIIflux[NIIIflux.mask], \
                     p0 = [1,4100,2,1,4340,3],sigma = NIIIsigma, bounds = (0,np.inf))
    
    
    heights = []
    heights_err = []
    means = []
    means_err = []
    sigmas = []
    sigmas_err = []
    vs = []
    variance = []
    
    heights.append(Haopt[0]) #Ha data
    means.append(Haopt[1])
    sigmas.append(Haopt[2])
    vs.append((Haopt[2]*2.355/Haopt[1])*c)
    heights_err.append((Hacov[0][0])**2)
    means_err.append((Hacov[1][1])**2)
    sigmas_err.append((Hacov[2][2])**2)
    variance.append(np.mean(1/ivar[Ha_center - 300:Ha_center + 300]))
    
    
    heights.append(HBopt[3]) #HB data
    means.append(HBopt[4])
    sigmas.append(HBopt[5])
    vs.append((HBopt[5]*2.355/HBopt[4])*c)
    heights_err.append((HBcov[3][3])**2)
    means_err.append((HBcov[4][4])**2)
    sigmas_err.append((HBcov[5][5])**2)
    variance.append(np.mean(1/ivar[HB_center - 500:HB_center + 500]))
    
    heights.append(NIIIopt[3]) #Hgamma data
    means.append(NIIIopt[4])
    sigmas.append(NIIIopt[5])
    vs.append((NIIIopt[5]*2.355/NIIIopt[4])*c)
    heights_err.append((NIIIcov[3][3])**2)
    means_err.append((NIIIcov[4][4])**2)
    sigmas_err.append((NIIIcov[5][5])**2)
    variance.append(np.mean(1/ivar[NIII_center - 200:NIII_center + 200]))     
         
    heights.append(HBopt[0]) #HeII data
    means.append(HBopt[1])
    sigmas.append(HBopt[2])
    vs.append((HBopt[2]*2.355/HBopt[1])*c)
    heights_err.append((HBcov[0][0])**2)
    means_err.append((HBcov[1][1])**2)
    sigmas_err.append((HBcov[2][2])**2)
    variance.append(np.mean(1/ivar[HB_center - 500:HB_center + 500]))
    
    heights.append(HBopt[6]) #OIII data
    means.append(HBopt[7])
    sigmas.append(HBopt[8])
    vs.append((HBopt[8]*2.355/HBopt[7])*c)
    heights_err.append((HBcov[6][6])**2)
    means_err.append((HBcov[7][7])**2)
    sigmas_err.append((HBcov[8][8])**2)   
    variance.append(np.mean(1/ivar[HB_center - 500:HB_center + 500]))

    heights.append(NIIIopt[0])#NIII data
    means.append(NIIIopt[1])
    sigmas.append(NIIIopt[2])
    vs.append((NIIIopt[2]*2.355/NIIIopt[1])*c)
    heights_err.append((NIIIcov[0][0])**2)
    means_err.append((NIIIcov[1][1])**2)
    sigmas_err.append((NIIIcov[2][2])**2)
    variance.append(np.mean(1/ivar[NIII_center - 200:NIII_center + 200]))
    

    
    
    linetable = Table([lines,restwavelengths,heights,heights_err,means,means_err,sigmas,sigmas_err,\
                      vs,variance],names = ('Line','Wavelength','Height','e_Height','Mean','e_Mean','Sigma',\
                                   'e_Sigma','Velocity','Variance'))
    linetable = linetable.to_pandas()
    return linetable


    
    
    
def TDE_filter(linetable, flux):
    filter_pass = []
    score = 0
    lines = list(linetable['Line'])
    #flux input should be prior to subtraction of continuum
    blue_end = np.nanmean(flux[1000:3000])
    red_end = np.nanmean(flux[5000:7000])
    
    #Halpha
    i = lines.index('Halpha')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 3*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 800\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('Halpha')
    
    #HBeta
    i = lines.index('Hbeta')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 3*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 650\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('Hbeta')
        
    #Hgamma
    i = lines.index('Hgamma')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 3*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 450\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('Hgamma')
    #HeII4686
    i = lines.index('HeII4686')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 3*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 500\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('HeII4686')
    #NIII
    i = lines.index('NIII')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 3*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 350\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('NIII')

    
    #OIII
    i = lines.index('OIII')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 3*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 10\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score -= 1.1
   
    #Blueness
    if blue_end/red_end >= 2:
        score += 1
        filter_pass.append('Blue')
    
    return(score)
#####
def Hline_filter(linetable):
    filter_pass = []
    score = 0
    lines = list(linetable['Line'])
    #flux input should be prior to subtraction of continuum

    #Halpha
    i = lines.index('Halpha')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 5*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 75\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('Halpha')
    
    #HBeta
    i = lines.index('Hbeta')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 5*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 75\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('Hbeta')
        
    #Hgamma
    i = lines.index('Hgamma')
    if linetable['e_Height'][i]/linetable['Height'][i] < 0.05 and \
            linetable['Height'][i]> 5*linetable['Variance'][i] and \
            linetable['e_Sigma'][i]/linetable['Sigma'][i] < 0.05 and linetable['Velocity'][i] > 50\
            and abs(linetable['Wavelength'][i] - linetable['Mean'][i]) < 20 and \
            linetable['e_Mean'][i] < 30:
        score +=1
        filter_pass.append('Hgamma')
  
    return(score)