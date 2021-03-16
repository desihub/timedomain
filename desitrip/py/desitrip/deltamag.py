"""
Calculate the difference between spectral and fiber magnitudes
"""

import numpy as np
import speclite.filters

def delta_mag(cspectra, fibermap, select):
    """
    Parameters:
    -----------
    cspectra : astropy table
        Coadded spectra table to access wavelength and flux arrays
    fibermap : astropy table
        Fibermap table to access imaging flux information
    select : ndarray
        Boolean array to isolate valid targets and fibers

    Returns:
    --------
    fibermap : astropy table
        Updated fibermap table containing delta mag per filter
    """
    # Calculate fiber magnitudes from imaging data
    flux_g = fibermap[select]['FLUX_G']
    flux_r = fibermap[select]['FLUX_R']
    flux_z = fibermap[select]['FLUX_Z']
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

    return fibermap

