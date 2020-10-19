"""Example scipt to plot a LIGO skymap with angular error contours
and confidence level containing LIGO PDF.

Based on plotting script by Raamis Hussain, UW-Madison, 21 May 2019.
Updated SYB, Dec 2020.
"""

import os
import logging
import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    import meander  # pip install this to draw contours around skymap
except ImportError as e:
    log.error(e)
    raise SystemExit

from glob import glob
from urllib.request import urlretrieve
from matplotlib import cm
from astropy.io import fits

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def compute_quantiles(proportions, samples):
    """Get quantiles of HEALPix probability map.

    Parameters:
    -----------
    proportions: list
        list of containment level to make contours for.
        E.g [0.68,0.9]
    samples: array
        array of values read in from healpix map
        E.g samples = hp.read_map(file)

    Returns:
    --------
    levels: list
        List of map values corresponding to a given containment level.
    """
    levels = []
    sorted_samples = list(reversed(list(sorted(samples))))
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > proportion).tolist().index(True)
        level = (sorted_samples[level_index] + (sorted_samples[level_index+1] if level_index+1 < len(samples) else 0)) / 2.0
        levels.append(level)
    return levels


def compute_contours(proportions, samples):
    """Plot containment contour around desired level. E.g 90% containment of a
    PDF on a healpix map.

    Parameters:
    -----------
    proportions: list
        list of containment level to make contours for.
        E.g [0.68,0.9]
    samples: array
        array of values read in from healpix map
        E.g samples = hp.read_map(file)

    Returns:
    --------
    ra_list: list
        List of arrays containing RA values for desired contours [deg].
    dec_list: list
        List of arrays containing Dec values for desired contours [deg].
    """

    levels = []
    sorted_samples = list(reversed(list(sorted(samples))))
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > proportion).tolist().index(True)
        level = (sorted_samples[level_index] + (sorted_samples[level_index+1] if level_index+1 < len(samples) else 0)) / 2.0
        levels.append(level)
    contours_by_level = meander.spherical_contours(sample_points, samples, levels)

    ra_list = []; dec_list=[]
    for contours in contours_by_level:
        for contour in contours:
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            dec_list.append(90 - np.degrees(theta))
            ra_list.append(np.degrees(phi))
    
    return ra_list, dec_list


def plot_gwmap(lvc_healpix_file, levels=[0.5, 0.9], rot=255):
    """Plot the GW map with the DESI footprint in a Mollweide projection.
    
    Parameters
    ----------
    lvc_healpix_file : str
        Relative or absolute path to LIGO/Virgo HEALPix angular reconstruction file.
    levels : list
        List of credible interval thresholds, e.g., 0.5, 0.9, etc.
    rot : float
        Rotation angle (around z) for map projection, in degrees.
    
    Returns
    -------
    fig : matplotlib.Figure
        Figure object for accessing or saving a plot.
    """
    # Read metadata from FITS.
    hdus = fits.open(lvc_healpix_file)
    header = hdus[1].header

    # instruments = header['INSTRUME']
    distmean = header['DISTMEAN']
    diststd = header['DISTSTD']
    origin = header['ORIGIN']
    date = header['DATE-OBS']
    
    # Read HEALPix map.
    gwmap = hp.read_map(lvc_healpix_file)

    # Compute GW contours.
    prob64 = hp.pixelfunc.ud_grade(gwmap, 64) #reduce nside to make it faster
    prob64 = prob64/np.sum(prob64)
    pixels = np.arange(prob64.size)
    ra_contour, dec_contour = compute_contours(levels, prob64)

    cmap = mpl.cm.OrRd
    cmap.set_under('w')

    # Access DESI contours.
    nside = hp.pixelfunc.get_nside(gwmap)
    desi_mask = hp.read_map('desi_mask_nside{:04d}.fits'.format(nside))
    probs = np.copy(gwmap)
    probs[desi_mask == 0] = hp.UNSEEN

    fig = plt.figure(num=1, figsize=(10,6))
    hp.mollview(probs, fig=1, cbar=True, unit=r'$dp/d\Omega$ [deg$^{-2}$]',
                title='{} {}'.format(origin, date),
                min=0, max=2e-4, flip='astro', rot=rot, cmap=cmap)
    hp.graticule(ls=':', alpha=0.5, dpar=30, dmer=45)

    ramin, ramax = 1e99, -1e99
    decmin, decmax = 1e99, -1e99
    for i, (rc, dc) in enumerate(zip(ra_contour, dec_contour)):
        ramin = np.minimum(ramin, np.min(rc))
        ramax = np.maximum(ramax, np.max(rc))
        decmin = np.minimum(decmin, np.min(dc))
        decmax = np.maximum(decmax, np.max(dc))
        hp.projplot(rc, dc, lonlat=True, linewidth=1, c='k')
        
    logging.info('RA_min={:.1f}, RA_max={:.1f}, Dec_min={:.1f}, Dec_max={:.1f}'.format(ramin, ramax, decmin, decmax))

#    ax = plt.gca()

    # Label latitude lines.
    for _dec in [-60,-30,0,30,60]:
        vert = 'top' if _dec < 0 else 'bottom' if _dec > 0 else 'center'
        hp.projtext(180+rot, _dec, r'{0:+}$^\circ$'.format(_dec),
                    ha='left', va=vert, fontsize=12, lonlat=True)
        
    # Label longitude lines.
    for _ra in np.arange(0,360,45):
        hp.projtext(_ra, -40, r'{:d}$^\circ$'.format(_ra),
                    horizontalalignment='center', fontsize=12, lonlat=True)
    
#    fig = plt.gcf()
    return fig


if __name__ == '__main__':
    p = ArgumentParser(description='GW Event Plotter',
                       formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-d', '--display', dest='display',
                   action='store_true', default=False,
                   help='Display plot')

    sp = p.add_subparsers(title='subcommands', description='valid subcommands',
                          help='Additional help.')

    # Event mode: access GW event from GraceDB using ID number.
    pemode = sp.add_parser('eventmode', help='Access event from GraceDB.')
    pemode.add_argument('event_id', nargs=1,
                        type=str,
                        help='LVC event ID')
    pemode.add_argument('-t', '--type', dest='event_type',
                        default=None, type=str,
                        help='Event type [BBH, BNS, ..]')

    # File mode: access GW event in a downloaded file.
    pfmode = sp.add_parser('filemode', help='Access event from FITS file.')
    pfmode.add_argument('fitsfile', nargs=1,
                        type=str,
                        help='FITS file with GW maps.')

    args = p.parse_args()
    
    # Fits file containing Gravitational wave skymap
    if 'event_id' in args:
        event_id = args.event_id[0]
        fitsFile = 'bayestar_{}.fits.gz'.format(event_id)
        if not os.path.isfile(fitsFile):
            url = 'https://gracedb.ligo.org/api/superevents/{}/files/bayestar.fits.gz'.format(event_id)
            print('Downloading {}'.format(url))
            urlretrieve(url, fitsFile)
    elif 'fitsfile' in args:
        fitsFile = args.fitsfile[0]
    
    fig = plot_gwmap(fitsFile)

    # Output plot to a PDF file.
    if 'event_id' in args:
        fig.savefig('gw_desi_{}.pdf'.format(event_id))
    elif 'fitsfile' in args:
        basename = fitsFile.split('.')[0]
        fig.savefig('gw_desi_{}.pdf'.format(basename))
    
    # Display results.
    if args.display:
        plt.show()
