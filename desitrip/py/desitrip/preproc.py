"""Spectrum preprocessing routines before object classification. Assumes a
single observation of a spectrum.
"""

import numpy as np
from desispec.interpolation import resample_flux

def _rebin_logwave(wave, flux, ivar, targetids, minwave=3500., maxwave=10000., dlogwave=1e-3):
    """
    Parameters:
    -----------
    wave : ndarray
        Input wavelength binning, assumed in Angstroms.
    flux : ndarray
        Input spectrum.
    ivar : ndarray
        Spectrum weight (inverse variance).
    targetids : ndarray
        List of TARGETIDs for spectra.
    minwave : float
        Minimum output wavelength, in Angstroms.
    maxwave : float
        Maximum output wavelength, in Angstroms.
    dlogwave : float
        Log wavelength bin size (per Angstrom).

    Returns:
    --------
    fl : ndarray
        Rebinned spectrum.
    iv : ndarray
        Rebinned inverse variance.
    """
    # Copied from QuasarNET; for checks.
    logmin = np.log10(minwave)
    logmax = np.log10(maxwave)
    logwave = np.log10(wave)
    logbins = np.floor((logwave-logmin)/dlogwave).astype(int)

    # Select wavelengths with data. Zero-pad everything else.
    w = (logwave>logmin) & (logwave<logmax)
    fl_aux = flux[:,w]
    iv_aux = ivar[:,w]
    ivfl_aux = fl_aux*iv_aux

    # Compute rebinned spectrum.
    utids = np.unique(targetids)
    nbins = int((logmax-logmin)/dlogwave)
    nspec = len(utids)
    fl = np.zeros((nspec, nbins))
    iv = np.zeros((nspec, nbins))

    for i, t in enumerate(targetids):
        j = np.argwhere(utids == t)[0]
        c = np.bincount(logbins, weights=ivfl_aux[i])
        fl[j,:len(c)] += c

        c = np.bincount(logbins, weights=iv_aux[i])
        iv[j,:len(c)] += c

    w = iv > 0
    fl[w] /= iv[w]

    return fl, iv

def rebin_flux(wave, flux, ivar=None, minwave=3600., maxwave=9900., nbins=600,  log=False):
    """Rebin differential flux vs wavelength using desispec resample_flux.

    Parameters
    ----------
    wave : ndarray
        Input wavelength; assume units of Angstroms.
    flux : ndarray
        Input differential spectra as a function of wavelength.
    ivar : None or ndarray
        Inverse variance (weight) of spectra vs wavelength.
    minwave : float
        Minimum output wavelength, in units of Angstroms.
    maxwave : float
        Maximum output wavelength, in units of Angstroms.
    nbins : int
        Number of output wavelength bins.
    log : bool
        If true, use logarithmic bins between minwave and maxwave.

    Returns
    -------
    basewave : ndarray
        Output wavelength, in units of Angstroms.
    fl : ndarray
        Rebinned spectra.
    iv : ndarray
        Rebinned inverse variance.
    """
    if log:
        basewave = np.logspace(np.log10(minwave), np.log10(maxwave), nbins)
    else:
        basewave = np.linspace(minwave, maxwave, nbins)

    if flux.ndim > 1:
        nspec = len(flux)
        fl = np.zeros((nspec, nbins))
        iv = np.ones((nspec, nbins))
        for i in range(nspec):
            if ivar is not None:
                iv[i], fl[i] = resample_flux(basewave, wave, flux[i], ivar[i])
            else:
                fl[i] = resample_flux(basewave, wave, flux[i])
    else:
        resampled = resample_flux(basewave, wave, flux, ivar)
        if ivar is not None:
            fl, iv = resampled
        else:
            fl, iv = resampled, None

    return basewave, fl, iv

def rescale_flux(flux):
    """Rescale flux so that it ranges from 0 to 1.

    Parameters
    ----------
    flux : ndarray
        Input flux array.

    Returns
    -------
    rsfl : ndarray
        Flux rescaled to range between 0 and 1.
    """
    if flux.ndim > 1:
        a, b = np.min(flux,axis=1)[:,None], np.max(flux,axis=1)[:,None]
    else:
        a, b = np.min(flux), np.max(flux)
    
    return (flux - a) / (b - a)
