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
    # Calculate fiber magnitudes from imaging data
    flux_g = fibermap[select]['FIBERFLUX_G']
    flux_r = fibermap[select]['FIBERFLUX_R']
    flux_z = fibermap[select]['FIBERFLUX_Z']
    fiber_gmag = 22.5 - 2.5*np.log10(flux_g)
    fiber_rmag = 22.5 - 2.5*np.log10(flux_r)
    fiber_zmag = 22.5 - 2.5*np.log10(flux_z)

    # Get flux and wavelength information
    wave = cspectra.wave['brz']
    flux = 1e-17 * cspectra.flux['brz'][select]

    # Load filter response information
    gband = speclite.filters.load_filter('decam2014-g')
    rband = speclite.filters.load_filter('decam2014-r')
    zband = speclite.filters.load_filter('decam2014-z')

    # Pad wavelength coverage to match filters
    gflux, gwave = gband.pad_spectrum(flux, wave)
    rflux, rwave = rband.pad_spectrum(flux, wave)
    zflux, zwave = zband.pad_spectrum(flux, wave)

    # Calculate spectral (AB) magnitudes
    spec_gmag = gband.get_ab_magnitude(gflux, gwave)
    spec_rmag = rband.get_ab_magnitude(rflux, rwave)
    spec_zmag = zband.get_ab_magnitude(zflux, zwave)

    # Calculate difference between spectral and fiber magnitudes
    delta_gmag = spec_gmag - fiber_gmag
    delta_rmag = spec_rmag - fiber_rmag
    delta_zmag = spec_zmag - fiber_zmag

    # Add magnitude difference to fibermap
    fibermap['DELTAMAG_G'] = np.empty(len(fibermap))
    fibermap['DELTAMAG_R'] = np.empty(len(fibermap))
    fibermap['DELTAMAG_Z'] = np.empty(len(fibermap))
    fibermap['DELTAMAG_G'][select] = delta_gmag
    fibermap['DELTAMAG_R'][select] = delta_rmag
    fibermap['DELTAMAG_Z'][select] = delta_zmag

    # Define deltamag cutoff for transient candidates
    nsigma = nsigma
    gmean = np.median(delta_gmag)
    rmean = np.median(delta_rmag)
    zmean = np.median(delta_zmag)
    gstd = nsigma * np.std(delta_gmag)
    rstd = nsigma * np.std(delta_rmag)
    zstd = nsigma * np.std(delta_zmag)
    dg_min = gmean - gstd
    dr_min = rmean - rstd
    dz_min = zmean - zstd

    # Find target ID of transient candients
    g_outliers = np.where(delta_gmag <= dg_min)
    r_outliers = np.where(delta_rmag <= dr_min)
    z_outliers = np.where(delta_zmag <= dz_min)
    outliers = reduce(np.union1d, (g_outliers,r_outliers,z_outliers))

    # Add transient candidates to fibermap
    fibermap['CANDIDATE'] = np.zeros(len(fibermap), dtype='bool')
    fibermap['CANDIDATE'][outliers] = True

    return fibermap

