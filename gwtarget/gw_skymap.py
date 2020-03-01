"""Example scipt to plot a LIGO skymap with angular error contours
and confidence level containing LIGO PDF.

Based on plotting script by Raamis Hussain, UW-Madison, 21 May 2019.
Updated SYB, Dec 2020.
"""

import os
import numpy  as np
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

from argparse import ArgumentParser


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
    theta_list: list
        List of arrays containing theta values for desired contours
    phi_list: list
        List of arrays containing phi values for desired contours
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

    theta_list = []; phi_list=[]
    for contours in contours_by_level:
        for contour in contours:
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            theta_list.append(theta)
            phi_list.append(phi)

    return theta_list, phi_list


def plot_gwmap(lvc_healpix_file, levels=[0.5, 0.9]):
    """Plot the GW map with the DESI footprint.
    
    Parameters
    ----------
    lvc_healpix_file : str
        Relative or absolute path to LIGO/Virgo HEALPix angular reconstruction file.
    levels : list
        List of credible interval thresholds, e.g., 0.5, 0.9, etc.
    
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
    date = header['DATE']
    
    # Read HEALPix map.
    gwmap = hp.read_map(lvc_healpix_file)

    # Compute GW contours.
    prob64 = hp.pixelfunc.ud_grade(gwmap, 64) #reduce nside to make it faster
    prob64 = prob64/np.sum(prob64)
    pixels = np.arange(prob64.size)
    #sample_points = np.array(hp.pix2ang(nside,pixels)).T
    theta_contour, phi_contour = compute_contours(levels, prob64)

    cmap = mpl.cm.OrRd
    cmap.set_under('w')

    # Access DESI contours.
    nside = hp.pixelfunc.get_nside(gwmap)
    desi_mask = hp.read_map('desi_mask_nside{:04d}.fits'.format(nside))
    probs = np.copy(gwmap)
    probs[desi_mask == 0] = hp.UNSEEN

    hp.mollview(probs, cbar=True, unit=r'probability', title='{} {}'.format(origin, date),
                min=0, max=2e-4, flip='astro', rot=180, cmap=cmap)
    hp.graticule(ls=':', alpha=0.5, dpar=30, dmer=45)

    thmin, thmax = 1e99, -1e99
    phmin, phmax = 1e99, -1e99
    for i, (tc, pc) in enumerate(zip(theta_contour, phi_contour)):
        thmin = np.minimum(thmin, np.min(tc))
        thmax = np.maximum(thmax, np.max(tc))
        phmin = np.minimum(phmin, np.min(pc))
        phmax = np.maximum(phmax, np.max(pc))
        hp.projplot(tc, pc, linewidth=1, c='k')

    ramin, ramax = np.degrees(phmin), np.degrees(phmax)
    decmin, decmax = 90-np.degrees(thmax), 90-np.degrees(thmin)

    ax = plt.gca()

    # Label latitude lines.
    ax.text( 2.00,  0.10, r'$0^\circ$', horizontalalignment='left')
    ax.text( 1.80,  0.45, r'$30^\circ$', horizontalalignment='left')
    ax.text( 1.30,  0.80, r'$60^\circ$', horizontalalignment='left')
    ax.text( 1.83, -0.45, r'$-30^\circ$', horizontalalignment='left')
    ax.text( 1.33, -0.80, r'$-60^\circ$', horizontalalignment='left')
    ax.text(-2.00,  0.10, r'$0^\circ$', horizontalalignment='right')
    ax.text(-1.80,  0.45, r'$30^\circ$', horizontalalignment='right')
    ax.text(-1.30,  0.80, r'$60^\circ$', horizontalalignment='right')
    ax.text(-1.85, -0.45, r'$-30^\circ$', horizontalalignment='right')
    ax.text(-1.35, -0.80, r'$-60^\circ$', horizontalalignment='right')

    # Label longitude lines.
    ax.text( 2.0, -0.15, r'0$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 1.5, -0.15, r'3$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 1.0, -0.15, r'6$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 0.5, -0.15, r'9$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 0.0, -0.15, r'12$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-0.5, -0.15, r'15$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-1.0, -0.15, r'18$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-1.5, -0.15, r'21$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-2.0, -0.15, r'24$^\mathrm{h}$', horizontalalignment='center')
    
    fig = plt.gcf()
    return fig


if __name__ == '__main__':
    parser = ArgumentParser(description='GW Event Plotter')
    parser.add_argument('event_id', nargs=1,
                        type=str,
                        help='LVC event ID')
    parser.add_argument('-d', '--display', dest='display',
                        action='store_true', default=False,
                        help='Display plot')
    parser.add_argument('-t', '--type', dest='event_type',
                        default=None, type=str,
                        help='Event type [BBH, BNS, ..]')
    args = parser.parse_args()
    
    # Fits file containing Gravitational wave skymap
    event_id = args.event_id[0]
    fitsFile = 'bayestar_{}.fits.gz'.format(event_id)
    if not os.path.isfile(fitsFile):
        url = 'https://gracedb.ligo.org/api/superevents/{}/files/bayestar.fits.gz'.format(event_id)
        print('Downloading {}'.format(url))
        urlretrieve(url, fitsFile)
    
    # Read metadata.
    hdus = fits.open(fitsFile)
    header = hdus[1].header
    instruments = header['INSTRUME']
    distmean = header['DISTMEAN']
    diststd = header['DISTSTD']
    
    # Read map and get probabilities
    probs = hp.read_map(fitsFile)
    nside = hp.pixelfunc.get_nside(probs)
    
    # Choose color map and set background to white
    cmap = cm.OrRd
    cmap.set_under("w")
    
    # Compute GW contours.
    prob64 = hp.pixelfunc.ud_grade(probs, 64) #reduce nside to make it faster
    prob64 = prob64/np.sum(prob64)
    pixels = np.arange(prob64.size)
    #sample_points = np.array(hp.pix2ang(nside,pixels)).T
    levels = [0.50, 0.90]
    theta_contour, phi_contour = compute_contours(levels, prob64)
    
    # Access DESI contours.
    desi_mask = hp.read_map('desi_mask_nside{:04d}.fits'.format(nside))
    probs[desi_mask == 0] = hp.UNSEEN
    
    # Plot GW skymap in Mollweide projection
    hp.mollview(probs, cbar=True, unit=r'probability', min=0, max=3e-5, rot=180, cmap=cmap)
    hp.graticule(ls=':', alpha=0.5, dpar=30, dmer=45) # Set grid lines
    
    # Draw containment contour around GW skymap
    nregion = len(theta_contour) // len(levels)
    # ls = ['-', '--', '-.']
    # label = ''
    for i, (tc, pc) in enumerate(zip(theta_contour, phi_contour)):
        hp.projplot(tc, pc, linewidth=1, c='k')
    #     j = i // nregion
    #     print(len(theta_contour), nregion, i, j)
    #     newlabel = '{:g}% credible region'.format(100*levels[j])
    #     if newlabel == label:
    #         hp.projplot(tc, pc, linewidth=1, c='k', linestyle=ls[j])
    #     else:
    #         hp.projplot(tc, pc, linewidth=1, c='k', linestyle=ls[j], label=newlabel)
    #         label = newlabel
    
    ax = plt.gca()
    
    # Label latitude lines.
    ax.text( 2.00,  0.10, r'$0^\circ$', horizontalalignment='left')
    ax.text( 1.80,  0.45, r'$30^\circ$', horizontalalignment='left')
    ax.text( 1.30,  0.80, r'$60^\circ$', horizontalalignment='left')
    ax.text( 1.83, -0.45, r'$-30^\circ$', horizontalalignment='left')
    ax.text( 1.33, -0.80, r'$-60^\circ$', horizontalalignment='left')
    ax.text(-2.00,  0.10, r'$0^\circ$', horizontalalignment='right')
    ax.text(-1.80,  0.45, r'$30^\circ$', horizontalalignment='right')
    ax.text(-1.30,  0.80, r'$60^\circ$', horizontalalignment='right')
    ax.text(-1.85, -0.45, r'$-30^\circ$', horizontalalignment='right')
    ax.text(-1.35, -0.80, r'$-60^\circ$', horizontalalignment='right')
    
    # Label longitude lines.
    ax.text( 2.0, -0.15, r'0$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 1.5, -0.15, r'3$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 1.0, -0.15, r'6$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 0.5, -0.15, r'9$^\mathrm{h}$', horizontalalignment='center')
    ax.text( 0.0, -0.15, r'12$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-0.5, -0.15, r'15$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-1.0, -0.15, r'18$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-1.5, -0.15, r'21$^\mathrm{h}$', horizontalalignment='center')
    ax.text(-2.0, -0.15, r'24$^\mathrm{h}$', horizontalalignment='center')
    
    # Label participating detectors.
    ax.text(-1.75, -1.1, instruments, horizontalalignment='left', fontsize=12)
    
    if args.event_type is not None:
        plt.title('{} : {}, $D_L={:.1f}\\pm{:.1f}$ Mpc'.format(event_id, args.event_type, distmean, diststd), fontsize=12)
    else:
        plt.title('{} : $D_L={:.1f}\\pm{:.1f}$ Mpc'.format(event_id, distmean, diststd), fontsize=12)
    #plt.legend(loc='lower right', bbox_to_anchor=(1.02, -0.075), fontsize=9,
    #           frameon=False, facecolor=None)
    plt.savefig('gw_desi_{}.pdf'.format(event_id))
    
    if args.display:
        plt.show()
