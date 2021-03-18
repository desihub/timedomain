"""
Calculate the difference between spectral and fiber magnitudes
"""

import numpy as np
import speclite.filters
from functools import reduce

def delta_mag(cspectra, fibermap, select, nsigma):
    """
    Parameters:
    -----------
    cspectra : astropy table
        Coadded spectra table to access wavelength and flux arrays
    fibermap : astropy table
        Fibermap table to access imaging flux information
    select : ndarray
        Boolean array to isolate valid targets and fibers
    nsigma : int
        n sigma deviation from mean deltamag to evaluate transient candidates

    Returns:
    --------
    fibermap : astropy table
        Updated fibermap table containing delta mag per filter
        and boolean array of potential transient candidates
    """
    # Set up deltamag and candidate output arrarys
    for i in 'GRZ':
        fibermap['DELTAMAG_{}'.format(i)] = np.empty(len(fibermap))
        fibermap['DELTAMAG_{}'.format(i)][:] = np.NaN
    fibermap['CANDIDATE'] = np.zeros(len(fibermap), dtype='bool')

    outliers = []
    arm_avail = [key for key in cspectra.wave.keys()][0]
    # Loop through spectrograph arms with available information
    for arm in arm_avail:
        # Imaging data uses 'g' not 'b'
        if arm == 'b':
            arm = 'g'

        # Calculate fiber magnitudes from imaging data
        fiberflux = fibermap[select]['FIBERFLUX_{}'.format(arm.upper())]
        fibermag = 22.5 - 2.5*np.log10(fiberflux)

        # Get flux and wavelength information per filter
        wave = cspectra.wave[arm_avail]
        flux = 1e-17 * cspectra.flux[arm_avail][select]
        band = speclite.filters.load_filter('decam2014-{}'.format(arm))
        flux, wave = band.pad_spectrum(flux, wave)

        # Calculate magnitude difference and add to fibermap
        specmag = band.get_ab_magnitude(flux, wave)
        deltamag = specmag - fibermag
        fibermap['DELTAMAG_{}'.format(arm.upper())][select] = deltamag

        # Define deltamag cutoff and find transient candidates
        mean = np.median(deltamag)
        std = nsigma * np.std(deltamag)
        thresh = mean - std
        outliers.append(np.where(deltamag <= thresh))

    # Add transient candidates to fibermap
    outliers = reduce(np.union1d, (outliers))
    fibermap['CANDIDATE'][outliers] = True

    return fibermap

