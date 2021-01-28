"""
Calculate the difference between fiber and spectral magnitudes
"""

import numpy as np
from pkg_resources import resource_filename

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

    # Get filter response information from speclite
    response_file_g = resource_filename('speclite', 'data/filters/{}.ecsv'.format('decam2014-g'))
    response_file_r = resource_filename('speclite', 'data/filters/{}.ecsv'.format('decam2014-r'))
    response_file_z = resource_filename('speclite', 'data/filters/{}.ecsv'.format('decam2014-z'))

    # Get wavelength (convert from nm to angstroms) and response information from files
    gfile = np.loadtxt(response_file_g, skiprows=15, unpack=True)
    gwave = 10. * gfile[0]
    gresponse = gfile[1]
    rfile = np.loadtxt(response_file_r, skiprows=15, unpack=True)
    rwave = 10. * rfile[0]
    rresponse = rfile[1]
    zfile = np.loadtxt(response_file_z, skiprows=15, unpack=True)
    zwave = 10. * zfile[0]
    zresponse = zfile[1]

    # Extrapolate to cframe wavelength grid as filter data have per nm resolution
    wave = cspectra.wave['brz']
    flux = 1e-17 * cspectra.flux['brz'][select]
    # g band
    gres = np.zeros(wave.shape)
    for wg in range(gresponse.shape[0]):
        if wg >= 1 and wg <= gresponse.shape[0]-2:
            lo = (gwave[wg] - gwave[wg-1]) / 2
            wlo = gwave[wg] - lo
            indlo = np.abs(wave - wlo).argmin()
            hi = (gwave[wg+1] - gwave[wg]) / 2
            whi = gwave[wg] + hi
            indhi = np.abs(wave - whi).argmin()
            gres[indlo:indhi] = gresponse[wg]
    # r band
    rres = np.zeros(wave.shape)
    for wr in range(rresponse.shape[0]):
        if wr >= 1 and wr <= rresponse.shape[0]-2:
            lo = (rwave[wr] - rwave[wr-1]) / 2
            wlo = rwave[wr] - lo
            indlo = np.abs(wave - wlo).argmin()
            hi = (rwave[wr+1] - rwave[wr]) / 2
            whi = rwave[wr] + hi
            indhi = np.abs(wave - whi).argmin()
            rres[indlo:indhi] = rresponse[wr]
    # z band
    zres = np.zeros(wave.shape)
    for wz in range(zresponse.shape[0]):
        if wz >= 1 and wz <= zresponse.shape[0]-2:
            lo = (zwave[wz] - zwave[wz-1]) / 2
            wlo = zwave[wz] - lo
            indlo = np.abs(wave - wlo).argmin()
            hi = (zwave[wz+1] - zwave[wz]) / 2
            whi = zwave[wz] + hi
            indhi = np.abs(wave - whi).argmin()
            zres[indlo:indhi] = zresponse[wz]

    # Convolve/integrate flux to get average flux density
    gflux = np.trapz(wave*flux*gres, wave) / np.trapz(wave*gres, wave)
    rflux = np.trapz(wave*flux*rres, wave) / np.trapz(wave*rres, wave)
    zflux = np.trapz(wave*flux*zres, wave) / np.trapz(wave*zres, wave)

    # Convert flux units to nanomaggies using effective wavelengths
    eff_gwave = np.trapz(wave*flux*gres, wave) / np.trapz(flux*gres, wave)
    eff_rwave = np.trapz(wave*flux*rres, wave) / np.trapz(flux*rres, wave)
    eff_zwave = np.trapz(wave*flux*zres, wave) / np.trapz(flux*zres, wave)
    specflux_g = (3.34e4 * eff_gwave**2 * gflux) / 3.631e-6
    specflux_r = (3.34e4 * eff_rwave**2 * rflux) / 3.631e-6
    specflux_z = (3.34e4 * eff_zwave**2 * zflux) / 3.631e-6

    # Convert flux to spectral magnitudes
    spec_gmag = 22.5 - 2.5*np.log10(specflux_g)
    spec_rmag = 22.5 - 2.5*np.log10(specflux_r)
    spec_zmag = 22.5 - 2.5*np.log10(specflux_z)

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


