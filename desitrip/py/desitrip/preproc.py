"""Spectrum preprocessing routines before object classification. Assumes a
single observation of a spectrum.
"""

import numpy as np

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
    """
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

    return fl
