from abc import ABC
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.constants import c
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

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
warnings.filterwarnings("ignore", category=RuntimeWarning)
#warnings.filterwarnings("ignore", category=OptimizeWarning)
mpl.rc('font', size=14)

class SpectralLine(ABC):

    def __init__(self, name, wavelength):

        self.name = name
        self.line = wavelength

    def get_linewidth(self, wave, flux):
        pass
#Just a naive gaussian definition for later fitting
def gaus(x, sigma, line, A, b):
    flux = A*np.exp(-(x - line)**2/(2*sigma)**2) + b
    return flux

class SpectralLine_pEW(SpectralLine):

    def __init__(self, name, wavelength, cont_blue=None, cont_red=None):
        """Initialize a spectral line used for a pEW calculation.
        
        Parameters
        ----------
        line : str
            Name of line.
        wavelength : float
            restframe wavelength of line.
        cont_blue : list of two floats, or None
            Lower and upper wavelength for computing continuum on 'blue' end of line.
        cont_red : list of two floats, or None
            Lower and upper wavelength for computing continuum on 'red' end of line.
        """
        super().__init__(name, wavelength)

        if cont_blue is not None:
            self.cont_blue = cont_blue
        else:
            self.cont_blue = [self.line - 205., self.line - 195.]

        if cont_red is not None:
            self.cont_red = cont_red
        else:
            self.cont_red = [self.line + 195., self.line + 205.]
            
        self.vshift = None
        self.roi_wave = None
        self.roi_flux = None
        self.roi_flux_scaled = None
        self.roi_continuum = None

    def get_roi(self, wave, flux, dminwave=200., dmaxwave=200.):
        """Produce a region of interest for estimating the line pEW.
        
        Parameters
        ----------
        wave : ndarray
            Input wavelength from a spectrum.
        flux : ndarray
            Input flux from a spectrum.
        dminwave : float
            Distance below the line center (in Angstrom) for computing the ROI.
        dmaxwave : float
            Distance above the line center (in Angstrom) for computing the ROI.
        
        Returns
        -------
        vshift : ndarray
            Array of velocities in the ROI w.r.t. to the line rest frame.
        roi_wave : ndarray
            Selected wavelengths in the ROI.
        roi_flux : ndarray
            Selected flux values in the ROI.
        """
        if wave[-1] < self.line - dminwave or wave[1] > self.line + dmaxwave:
            raise ValueError('Analysis region not covered for {} line'.format(self.name))
        
        # Get velocity shift
        vshift = c * (wave / self.line - 1.)
    
        # Select the wavelengths in the region of interest around the line.
        select = (wave > self.line - dminwave) & (wave < self.line + dmaxwave)
        
        return vshift[select], wave[select], flux[select]
    
    def estimate_continuum(self, roi_wave, roi_flux):
        """Perform a linear fit of the continuum using 'blue' and 'red' sidebands of the spectra line.
        
        Parameters
        ----------
        roi_wave : ndarray
            Wavelengths in the region of interest of this spectral line.
        roi_flux : ndarray
            Fluxes in the region of interest of the spectral line.
        
        Returns
        -------
        a : float
            Intercept of the linear fit to the continuum.
        b : float
            Slope of the linear fit to the continuum.
        continuum : ndarray
            Computed continuum under the spectral line using a linear model.
        scaled_roi_flux :
            Flux in the ROI, scaled by the continuum to remove its slope.
        """
        blue_min, blue_max = self.cont_blue
        blue_avg = 0.5*(blue_min + blue_max)
        select_blue = (roi_wave > blue_min) & (roi_wave < blue_max)
        flux_blue = roi_flux[select_blue]
        
        red_min, red_max = self.cont_red
        red_avg = 0.5*(red_min + red_max)
        select_red = (roi_wave > red_min) & (roi_wave < red_max)
        flux_red = roi_flux[select_red]
        
        # Perform a linear fit to the mean flux vs mean wavelength in the blue and red sidebands of the line.
        # a = intercept, b = slope.
        b = (np.nanmean(flux_red) - np.nanmean(flux_blue)) / (red_avg - blue_avg)
        a = np.nanmean(flux_blue) - (b * blue_avg)

        # Calculate the continuum using the linear fit.
        continuum = a + b*roi_wave
        
        # Scale the flux by the linear continuum.
        scaled_roi_flux = roi_flux / continuum
        return a, b, continuum, scaled_roi_flux
    
    def estimate_peqw(self, roi_wave, scaled_roi_flux):
        """Compute the pseudo-equivalent width of the spectral line.
        
        Parameters
        ----------
        roi_wave : ndarray
            Wavelength in the fit region of interest of the line.
        scaled_roi_flux : ndarray
            Flux in the ROI scaled by the continuum to remove its slope.
            
        Returns
        -------
        """
        continuum = np.ones_like(scaled_roi_flux)
        _sum = np.sum((scaled_roi_flux - continuum) / continuum)
        _delta = np.diff(roi_wave)[0]
        return _sum * _delta

    def get_linewidth(self, wave, flux):
        """Return the line width (pseudo-equivalent width).
        
        Parameters
        ----------
        wave : ndarray
            Input wavelength from a spectrum.
        flux : ndarray
            Input flux from a spectrum.
        
        Returns
        -------
        peqw : float
            Pseudo-equivalent width of the spectral line.
        """
        vs, roi_wave, roi_flux = self.get_roi(wave, flux,self.line - np.nanmean(self.cont_blue),\
                                              np.nanmean(self.cont_red)-self.line)
        a, b, continuum, scaled_roi_flux = self.estimate_continuum(roi_wave, roi_flux)
        peqw = self.estimate_peqw(roi_wave, scaled_roi_flux)
        
        # Store the details of the continuum fit for debugging.
        self.vshift = vs
        self.roi_wave = roi_wave
        self.roi_flux = roi_flux
        self.roi_flux_scaled = scaled_roi_flux
        self.roi_continuum = continuum
        
        return peqw

    
    def __str__(self):
        s = '{}: {}\nBlue continuum fit: {} - {}\nRed continuum fit:  {} - {}'.format(
            self.name, self.line,
            self.cont_blue[0], self.cont_blue[1],
            self.cont_red[0], self.cont_red[1])
        return s
    

def is_TDE(Ha,HB,HeII,OIII,NIII, flux):
    '''
    Ha,HB, HeII, and OIII should be pEW for each of those lines respectively.
    '''
    TDE_score = 0 
    featurelist = []
    blue_end = np.nanmean(flux[1000:3000])
    red_end = np.nanmean(flux[5000:7000])
    mean_flux = np.nanmean(flux)
    try:
        if OIII/Ha > 0.05 or OIII > 1:
            TDE_score -= 1.1 #if OIII is present, will be marked by an index with a noninteger part
    except ZeroDivisionError:
        if OIII > 0.5:
            TDE_score -= 1.1
    if Ha > 25: #Just looks for generally wide Ha lines
        TDE_score += 1
        featurelist.append('Ha')
    if HB > 20:
        TDE_score += 1
        featurelist.append('HB')
                        
    try:
        if HeII/Ha >0.75 and HeII > 15: #looks for comparable He and Ha lines
            TDE_score += 1
            featurelist.append('HeII')
    except ZeroDivisionError:
        pass
    if NIII > 15: #looks for NIII sometimes seen in TDE-bowen 
        TDE_score += 1
        featurelist.append('NIII')
    if blue_end/red_end >= 1.5:
        TDE_score += 1
        featurelist.append('Blue')
    if mean_flux < 1.5:
        TDE_score -= 1.1
    
    return TDE_score, featurelist

def Combine_multifilt(wave,flux, mask = False):
        difwave_single = []
        difflux_single = []
        difmask_single = []
        for band in wave:
            difwave_single += list(wave[band]) 
            difflux_single += list(flux[band])
            if mask:
                difmask_single += list(mask[band])
                
        if mask:
            diftable = Table([difwave_single, difflux_single,difmask_single], names = ('wave','flux', 'mask'))
        else:
            diftable = Table([difwave_single, difflux_single], names = ('wave','flux'))
        diftable.sort('wave')
        difwave_single = list(diftable['wave'])
        difflux_single = list(diftable['flux'])
        if mask:
            difmask_single = list(diftable['mask'])
            return difwave_single,difflux_single,difmask_single
        else:
            return difwave_single,difflux_single

    
def Get_Spectra(tilepath):
    specdict = {}
    specs= sorted(glob('{}/spectra*.fits'.format(tilepath)))
    for spec in specs:
        spectra = read_spectra(spec)

        zname = spec.replace('spectra','redrock')
        rrdata = Table.read(zname)['TARGETID','Z','ZWARN','SPECTYPE','DELTACHI2']
        select = (rrdata['SPECTYPE'] == 'GALAXY') & (rrdata['ZWARN'] == 0) & \
        (rrdata['DELTACHI2'] >= 25) & (rrdata['TARGETID'] > 0  & (rrdata['Z'] > 0))
        spectra = coadd_cameras(spectra)
        rrdata = rrdata[select]
#         select2 = (spectra.fibermap['TARGETID'] in rrdata['TARGETID']) & (spectra.fibermap['OBJTYPE'] == 'TGT')
#         spectra = spectra[select2]
        for i in range(len(spectra.fibermap)):
            tid = spectra.exp_fibermap[i]['TARGETID']
            if  tid in rrdata['TARGETID'] and spectra.fibermap[i]['OBJTYPE'] == 'TGT':
                z_ind = list(rrdata['TARGETID']).index(tid)
                specdict[tid] = {'Z':rrdata['Z'][z_ind]}
                for band in spectra.bands:
                    specdict[tid]['wave'] = spectra.wave[band]
                    specdict[tid]['flux'] = spectra.flux[band][i]
                specdict[tid]['FILEINFO'] = spec #Format= petal-tile-date
                specdict[tid]['RA'] = spectra.fibermap[i]['TARGET_RA']
                specdict[tid]['DEC'] = spectra.fibermap[i]['TARGET_DEC']
                specdict[tid]['PETAL'] = spectra.exp_fibermap[i]['PETAL_LOC']
                specdict[tid]['TILE'] = spectra.exp_fibermap[i]['TILEID']
    return(specdict)
def TDE_Check(wave,flux, redshift, scale = 3.,cont_width = 10,multifilt = False, spectype = 'GALAXY', mask = False):
    '''
    Creates a scaled flux (assumed to be roughly gaussian) for each line of interest. Then fits this scaled flux to a
    Gaussian, pulling the $\sigma$. Recalculates pEW using an ROI of the line $\pm$ scale*sigma. Appends this pEW to
    a Dictionary (lines), each line should have a list in format [wavelength, pEW]. pEW used in logic to determine whether 
    or not it's a TDE.
    Parameters:
    wave, flux: wavelength and flux of spectrum
    scale: # of standard deviations to fit continuum
    cont_width: width of region to sample cont_blue and cont_red, default 5 grabs area 2.5 angstroms on either side of 
        bound determined by scale*sigma
    '''
    if spectype == 'GALAXY':

        if multifilt:
            if mask:
                wave,flux, mask = Combine_multifilt(wave, flux,mask)
            else:
                wave,flux= Combine_multifilt(wave, flux)
        wave = wave/(1+redshift)
        flux = gaussian_filter1d(flux, sigma = 3)
        half_cont = cont_width/2
        lines = {'Halpha':[6562.79],'Hbeta':[4861.4], 'HeII4686':[4686],'OIII':[5007],'NIII':[4100]}
        for line in lines:
            pass1 = SpectralLine_pEW(line,lines[line][0])
            try:
                peqw_pass1 = pass1.get_linewidth(wave, flux)
                #plt.plot(pass1.roi_wave, pass1.roi_flux)
                #plt.show()
                #print(line)
                #print(peqw_pass1)

                try:
                    popt = curve_fit(gaus, pass1.roi_wave, pass1.roi_flux_scaled, p0 = [10.,pass1.line,max(pass1.roi_flux_scaled),1],check_finite = False)[0]
                    std = popt[0]
                    #print(std)
                    if scale*std <= half_cont: #Fit failed, line probably not present enough for measurement. 
                        peqw_specline = 0
                        #print('Failed')

                    else:
                        cont_blue=[pass1.line - scale*std-half_cont, pass1.line - scale*std+half_cont]
                        cont_red=[pass1.line + scale*std-half_cont, pass1.line + scale*std+half_cont]
                        specline = SpectralLine_pEW(line, lines[line][0], cont_blue=cont_blue, \
                                    cont_red=cont_red)
                        if mask:
                            roistart_ind = wave.index(specline.roi_wave[0])
                            roiend_ind = wave.index(specline.roi_wave[-1]) + 1
                            if 1 in mask[roistart_ind:roiend_ind]:
                                peqw_specline = 0
                                
                        else:
                            try:
                                peqw_specline = specline.get_linewidth(wave, flux)
                                #print(cont_blue,cont_red)
                                #print(peqw_specline)
                                #plt.plot(specline.roi_wave, specline.roi_flux)
                                #plt.show()
                            except IndexError:
                                peqw_specline = 0
                #print('###########')
                except RuntimeError:
                    peqw_specline = 0
            except:
                peqw_specline = 0
            lines[line].append(peqw_specline)
        return is_TDE(lines['Halpha'][-1],lines['Hbeta'][-1],lines['HeII4686'][-1],lines['OIII'][-1],\
                      lines['NIII'][-1], flux)
    else:
        return 0