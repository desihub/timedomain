#!/usr/bin/env python
# coding: utf-8

# # GW_transient_compare
# ### This now does many things but all are related to handling Gravitational Wave (GW) localization maps and DESI observations.
# ### Most generally, this finds those observations in the area of a specified confidence interval (CI) of the GW map.
# ### It also does the following:

# 1. Finds *not-previously-observed* ALERCE alerts in the area and creates a ToO ledger for those.
# 2. Identifies Bright and Dark targets in the CI contour.
# 3. Compares previous observations with dr9 targets in the CI contour to find the dr9 targets not already observed. Can perform this comparison by matching the targetids or by matching the RA's and DEC's which we feel is slightly more robust.
# 4. A bit of plotting but it takes awhile to produce a proper cartesian map, be forewarned.

import sys
import os
user_home = os.environ['HOME']
user_name = os.getlogin()
# For importing useful functions
sys.path.append(f"/global/homes/{user_name[0]}/{user_name}/timedomain/")
_ = [sys.path.append(f"{user_home}/desi/{x}/py/") for x in os.listdir(f"/global/homes/{user_name[0]}/{user_name}/desi/")]

from astropy.io import fits
from astropy.table import Table, Column, join, hstack, vstack, unique, setdiff
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, Angle
from astropy.time import Time

try:
    import astropy_healpix as ah
except:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "astropy-healpix"])
    import astropy_healpix as ah

import requests
from alerce.core import Alerce
from alerce.exceptions import APIError

from gw_skymap import plot_mollmap, compute_quantiles, compute_contours #, plot_cartmap
from desitarget import io, cuts
from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

import matplotlib.pyplot as plt
import healpy as hp
import pandas as pd
import numpy as np

from copy import deepcopy
from glob import glob
import argparse
import logging # to disable output when reading in FITS info
import subprocess

import psycopg2
import requests
import sqlite3

import pickle

# # Some handy, frequently used things
global db_filename
db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'

global exposure_path
# Hardcoding this for now until I can get the environment set up right
exposure_path = "/global/cfs/cdirs/desi/spectro/redux" #os.environ["DESI_SPECTRO_REDUX"]

global color_band
color_band = "r"

global today
today = Time.now()

global tile_rad
tile_rad = 1.6*u.deg

# ## Some function definitions which are hopefully self-explanatory

# From decam_TAMU_ledgermaker.ipynb - https://github.com/desihub/timedomain/blob/master/too_ledgers/decam_TAMU_ledgermaker.ipynb
# Thanks Antonella!
def write_too_ledger(filename, too_table, checker, overwrite=False, verbose=False, tabformat='TAMU'):
    """Write ToO ledger in the ECSV format specified by Adam Meyers.
    These can be passed to fiberassign for secondary targeting.
    
    Parameters
    ----------
    filename : str
        Output filename of the ledger (can be an absolute path).
    too_table : pandas.DataFrame
        Table of ToO data, using DECam format.
    checker : str
        Initials of individual(s) who have verified the ToO list.
    overwrite : bool
        If True, overwrite the output file.
    verbose : bool
        If True, 
    """
    today = Time.now()
    
    mode = 'w' if overwrite else 'a'
    if verbose:
        mode = mode + '+'
    
    if tabformat == 'ND':
        last_id = "TILEID"
    else:
        last_id = "TOOID"
        
    with open(filename, mode) as outf:
        if overwrite:
            outf.write("""# %ECSV 0.9
    # ---
    # datatype:
    # - {name: RA, unit: deg, datatype: float64}
    # - {name: DEC, unit: deg, datatype: float64}
    # - {name: PMRA, unit: mas / yr, datatype: float32}
    # - {name: PMDEC, unit: mas / yr, datatype: float32}
    # - {name: REF_EPOCH, unit: yr, datatype: float32}
    # - {name: CHECKER, datatype: string}
    # - {name: TOO_TYPE, datatype: string}
    # - {name: TOO_PRIO, datatype: string}
    # - {name: OCLAYER, datatype: string}
    # - {name: MJD_BEGIN, unit: d, datatype: float64}
    # - {name: MJD_END, unit: d, datatype: float64}
    # - {name: TOOID, datatype: int32}
    # - {name: PROB_COVERED, datatype: float64} if applicable
    # meta: {DEPNAM00: desitarget, DEPNAM01: desitarget-git, DEPVER00: 0.53.0.dev4635, DEPVER01: 0.53.0-24-g58c9a719, EXTNAME: TOO, RELEASE: 9999}
    # schema: astropy-2.0\n""")
            outf.write(f"RA DEC PMRA PMDEC REF_EPOCH CHECKER TOO_TYPE TOO_PRIO OCLAYER MJD_BEGIN MJD_END {last_id} PROB_COVERED\n")
            
        datedict = {}
        reporting = ['DESIRT','DDF','ALERCE', 'LEGACY']
        
        if tabformat=='TAMU':
            for i in range(too_table.shape[0]):

                row=too_table.iloc[i]           
                coord = SkyCoord(ra=row['RA-OBJECT'], dec=row['DEC-OBJECT'], unit=(u.degree, u.degree), frame='icrs')
                ra, dec = coord.ra.to('deg').value, coord.dec.to('deg').value
                t_disc = Time(row['Discovery-Time'], scale='utc')

                mag  = row['Discovery-Magnitude']
                too_type = 'FIBER'
                too_prog = 'BRIGHT' #if mag < 21 else 'DARK'
                too_prio = 'HI'

                # Encode the ToO ID as: MJD + ID + NNN.
                mjd_disc = int(t_disc.mjd) 
                reporter = 'DESIRT'
                if reporter not in reporting:
                    reporting.append(reporter)
                mjd_exp = 100*mjd_disc + reporting.index(reporter)
                if mjd_exp in datedict:
                    datedict[mjd_exp] += 1
                else:
                    datedict[mjd_exp] = 1
                too_id = 100*mjd_exp + datedict[mjd_exp]

                epoch = 2000.0

                outf.write('{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {}\n'.format(
                        ra, dec, 0, 0, epoch, checker, too_type, too_prio, too_prog, t_disc.mjd, today.mjd+14, too_id))
                
        if tabformat=='ddf':
            for i in range(too_table.shape[0]):

                row=too_table[i]           
                coord = SkyCoord(ra=row['RA'], dec=row['DEC'], unit=(u.degree, u.degree), frame='icrs')
                ra, dec = coord.ra.to('deg').value, coord.dec.to('deg').value
                t_disc = today.mjd #Do not have discovery time so using today

                too_type = 'FIBER'
                too_prog = 'BRIGHT' #if mag < 21 else 'DARK'
                too_prio = 'HI'

                # Encode the ToO ID as: MJD + ID + NNN.
                mjd_disc = int(t_disc) 
                reporter = 'DDF'
                reporting_id=2
                mjd_exp = 100*mjd_disc + reporting_id
                if mjd_exp in datedict:
                    datedict[mjd_exp] += 1
                else:
                    datedict[mjd_exp] = 1
                too_id = 100*mjd_exp + datedict[mjd_exp]

                epoch = 2000.0

                outf.write('{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {}\n'.format(
                        ra, dec, 0, 0, epoch, checker, too_type, too_prio, too_prog, t_disc, t_disc+14, too_id))
        
        # **************** My addition ****************
        if tabformat=='ALERCE':
            for i in range(too_table.shape[0]):

                row = too_table.iloc[i]           
                coord = SkyCoord(ra=row['meanra'], dec=row['meandec'], unit=(u.degree, u.degree))
                ra, dec = coord.ra.to('deg').value, coord.dec.to('deg').value
                t_disc = Time(row['lastmjd'], format = 'mjd')

                #mag  = row['Discovery-Magnitude']
                too_type = 'FIELD' # FIBER
                too_prog = 'BRIGHT' #if mag < 21 else 'DARK'
                too_prio = 'HI'

                # Encode the ToO ID as: MJD + ID + NNN.
                mjd_disc = int(t_disc.mjd) 
                reporter = 'ALERCE'
                if reporter not in reporting:
                    reporting.append(reporter)
                mjd_exp = 100*mjd_disc + reporting.index(reporter)
                
                if mjd_exp in datedict:
                    datedict[mjd_exp] += 1
                else:
                    datedict[mjd_exp] = 1
                    
                too_id = 100*mjd_exp + datedict[mjd_exp]

                epoch = 2000.0

                outf.write('{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {}\n'.format(
                        ra, dec, 0, 0, epoch, checker, too_type, too_prio, too_prog, t_disc.mjd, today.mjd+14, too_id))
    
        if tabformat == 'LEGACY': #DR9, Legacy Survey
            counter = 0
            for i in range(too_table.shape[0]):

                row = too_table.iloc[i]           
                coord = SkyCoord(ra=row['RA'], dec=row['DEC'], unit=(u.degree, u.degree))
                ra, dec = coord.ra.to('deg').value, coord.dec.to('deg').value

                #mag  = row['Discovery-Magnitude']
                too_type = 'FIELD' 
                
                # Should this be the original program?
                too_prog = 'BRIGHT' #if mag < 21 else 'DARK'
                too_prio = 'LOW'

                # Encode the ToO ID as: ID + NNNNN.
                #mjd_disc = int(t_disc.mjd) 
                reporter = 'LEGACY'
                
                '''
                if reporter not in reporting:
                    reporting.append(reporter)
                mjd_exp = 100*mjd_disc + reporting.index(reporter)
                
                 if mjd_exp in datedict:
                     datedict[mjd_exp] += 1
                 else:
                     datedict[mjd_exp] = 1
                     
                too_id = 100*mjd_exp + datedict[mjd_exp]
                ''' 
                too_id = 99 + counter
                counter += 1

                epoch = 2000.0

                outf.write('{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {}\n'.format(
                        ra, dec, 0, 0, epoch, checker, too_type, too_prio, too_prog, today.mjd, today.mjd+30, too_id))
    
        if tabformat == 'ND': #ND Mode
            counter = 0
            for i in range(too_table.shape[0]):

                row = too_table.iloc[i]           
                coord = SkyCoord(ra=row['RA'], dec=row['DEC'], unit=(u.degree, u.degree))
                ra, dec = coord.ra.to('deg').value, coord.dec.to('deg').value

                #mag  = row['Discovery-Magnitude']
                too_type = 'FIELD' 
                too_prog = 'BRIGHT' #if mag < 21 else 'DARK'
                too_prio = 'LOW'

                # Encode the ToO ID as: ID + NNNNN.
                # mjd_disc = int(t_disc.mjd) 
                # reporter = 'LEGACY'
                
                '''
                if reporter not in reporting:
                    reporting.append(reporter)
                mjd_exp = 100*mjd_disc + reporting.index(reporter)
                
                 if mjd_exp in datedict:
                     datedict[mjd_exp] += 1
                 else:
                     datedict[mjd_exp] = 1
                     
                too_id = 100*mjd_exp + datedict[mjd_exp]
                ''' 
                #too_id = 99 + counter
                #counter += 1
                
                tile_id = row['TILEID']
                prob = row['PROB_COVERED']

                epoch = 2000.0

                outf.write('{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {} {:>.3f}%\n'.format(
                        ra, dec, 0, 0, epoch, checker, too_type, too_prio, too_prog, today.mjd, today.mjd+30, tile_id, prob))
    
        if verbose:
            outf.seek(0)
            for line in outf:
                print(line.strip())
                
    return None

# Borrowed from gw_skymap.py with minor modifications

def plot_cartmap(lvc_healpix_file, levels=[0.5, 0.9], angsize=3., tile_ra=None, tile_dec=None, targ_ra=None, targ_dec=None):
    """Plot the GW map with the DESI footprint in a Cartesian projection.
    
    Parameters
    ----------
    lvc_healpix_file : str
        Relative or absolute path to LIGO/Virgo HEALPix angular reconstruction file.
    levels : list
        List of credible interval thresholds, e.g., 0.5, 0.9, etc.
    angsize : float
        Size of plot (-angsize, +angsize) in degrees about the center.
    tile_ra : list or ndarray
        List of RAs for DESI tiles (in deg).
    tile_dec : list or ndarray
        List of declinations for DESI tiles (in deg).
    targ_ra : list or ndarray
        List of RAs for DESI targets (in deg).
    targ_dec : list or ndarray
        List of declinations for DESI targets (in deg).
    
    Returns
    -------
    fig : matplotlib.Figure
        Figure object for accessing or saving a plot.
    """
    gwmap = hp.read_map(lvc_healpix_file)
    npix = len(gwmap)
    nside = hp.npix2nside(npix)

    # Compute contours.
    if nside > 256:
        _gwmap = hp.pixelfunc.ud_grade(gwmap, 256)
        _gwmap = _gwmap / np.sum(_gwmap)
    else:
        _gwmap = gwmap
    ra_contour, dec_contour = compute_contours(levels, _gwmap)

    # Create a temporary plot to produce a nice image array.
    # This code sets the size of the map around the maximum value.
    maxpix = np.argmax(gwmap)
    ra_c, dec_c = hp.pix2ang(nside, maxpix, lonlat=True)

    xmin = np.round(ra_c - angsize)
    xmax = np.round(ra_c + angsize)
    if xmax < xmin:
        xmin, xmax = xmax, xmin
    cxmin, cxmax = xmin, xmax
    frot = 0.
    if xmax > 90 and xmax < -90:
        frot, cxmin, cmax = 180., xmax-180., xmax+180.
    ymin = np.round(dec_c - angsize)
    ymax = np.round(dec_c + angsize)

    faspect = np.abs(cxmax - cxmin)/np.abs(ymax-ymin)
    fysize = 4
    figsize = (fysize*faspect+1, fysize+2)

    # Open and close the temporary plot.
    tfig   = plt.figure(num=2,figsize=figsize)
    rotimg = hp.cartview(gwmap, fig=2,coord='C', title="", cbar=False, flip='astro',
                         lonra=[cxmin,cxmax], latra=[ymin,ymax], rot=frot,
                         notext=True, xsize=1000,
                         return_projected_map=True)
    plt.close(tfig)

    # Now make the real plot with the desired angular contours.
    fig, ax = plt.subplots(1,1, num=1, figsize=figsize)
    img = ax.imshow(rotimg, extent=[cxmax, cxmin, ymin, ymax],
                    origin='lower', cmap='OrRd')

    for i, (rc, dc, lstyle, clev) in enumerate(zip(ra_contour, dec_contour, ['--', '-'], ['50', '90'])):
        p = ax.plot(rc, dc, 'g-', ls=lstyle, lw=2, label='{}% CI'.format(clev))

    ax.set(xlim=(cxmax, cxmin),
       xlabel='RA [deg]',
       ylabel='Dec [deg]')
    ax.grid(ls=':')

    _h, _l = ax.get_legend_handles_labels()

#     # Add DESI tile drawings, specified by central RA, Dec.
#     if tile_ra is not None and tile_dec is not None:
#         for _ra_c, _dec_c in zip(tile_ra, tile_dec):
#             circ = plt.Circle((_ra_c, _dec_c), radius=1.6, fc='None', ec='b', ls=':', lw=2)
#             ax.add_artist(circ)

#         _h.append(circ)
#         _l.append('DESI FOV')

    # Add DESI targets, specified by RA, Dec.
    if targ_ra is not None and targ_dec is not None:
        ax.plot(targ_ra, targ_dec, 'k.', alpha=0.5) # temp change, alpha = 0.1 -> alpha = 0.5 (maybe push command line arg for this)

    ax.legend(handles=_h, labels=_l, fontsize=10, ncol=2)

    cb = fig.colorbar(img, orientation='horizontal', shrink=0.95,
                      fraction=0.04, pad=0.2, ax=ax)
    cb.set_label(r'$dp/d\Omega$ [deg$^{-2}$]')

    return fig

# From ALeRCE_ledgermaker https://github.com/alercebroker/alerce_client
# I have had trouble importing this before so I copy, paste it, and modify it here.

# Choose cone_radius of diameter of tile so that, whatever coord I choose for ra_in, dec_in, we cover the whole tile
# KEPT FOR POSTERITY, MAY DELETE IN NEXT PUSH TO REPO
def access_alerts_old(radius = 3600*4, order_by = 'oid', order_mode = 'DESC', classifier = 'stamp_classifier', class_names=['SN', 'AGN'], **kwargs):
    alerce_client = Alerce()
    if type(class_names) is not list:
        raise TypeError('Argument `class_names` must be a list.')
        
    dataframes = []
#     if not lastmjd:
#         date_range = 60
#         lastmjd = [Time.now().mjd - date_range, Time.now().mjd]
#         print('Defaulting to a lastmjd range of', str(date_range), 'days before today.')
        
    #print("lastmjd:", lastmjd_in)
    for class_name in class_names:
        data = alerce_client.query_objects(classifier = classifier,
                                           class_name = class_name, 
                                           radius = radius, # in arcseconds
                                           page_size = 5000,
                                           order_by = order_by,
                                           order_mode = order_mode,                             
                                           format = 'pandas',
                                           **kwargs)
        
        #if lastmjd is not None:
        #    select = data['lastmjd'] >= lastmjd
        #    data = data[select]
            
        dataframes.append(data)
    
    #print(pd.concat(dataframes).columns)
    return pd.concat(dataframes)#.sort_values(by = 'lastmjd')

# Choose radius of diameter of tile (in degrees) so that, whatever coord I choose for ra_in, dec_in, we cover the whole tile
def access_alerts(ra = 0, dec = 0, radius = 4, order_by = 'oid', # radius = 4
                  order_mode = 'DESC', classifier = 'stamp_classifier', 
                  class_names = ['SN'], **kwargs):
    
    if not np.any(ra) or not np.any(dec):
        raise NameError('RAs and DECs must be fed in as key word arguments (ra =..., dec =...) to `access_alerts`.')
    else:
        objects_str = ",\n".join([f"({r}, {d})" for r,d in zip(ra, dec)])
    
    if not isinstance(class_names, (list,tuple)):
        raise TypeError('Argument `class_names` in `access_alerts` must be a list or a tuple.')
        
    # ***** Backwards compatibility *****
    # May remove in next push
    if ("firstmjd" in kwargs) and ("lastmjd" in kwargs):
        raise AttributeError('Can only accept argument `firstmjd` or `lastmjd` not both in `access_alerts`.')
        
    if "firstmjd" in kwargs:
        mjd_choose = "firstmjd"
        days_forward = kwargs.get("days_forward", 30)
        days_backward = 0
        mjd_date = kwargs.get("firstmjd", Time.now().mjd - days_forward)
        if isinstance(mjd_date, list):
            mjd_date = mjd_date[0]
        
    elif "lastmjd" in kwargs:
        mjd_choose = "lastmjd"
        days_forward = 0
        days_backward = kwargs.get("days_backward", 60)
        mjd_date = kwargs.get("lastmjd", Time.now().mjd + days_backward)
        if isinstance(mjd_date, list):
            mjd_date = mjd_date[1]
        
    else:
        raise KeyError('Please specify either `firstmjd` or `lastmjd` in `access_alerts`.')
        
    credentials_file = "https://raw.githubusercontent.com/alercebroker/usecases/master/alercereaduser_v4.json"
    params = requests.get(credentials_file).json()["params"]
    conn = psycopg2.connect(dbname=params["dbname"], user=params["user"], host=params["host"], password=params["password"])
    

    matches = pd.DataFrame()
    
    for class_name in class_names:
        prefix = class_name.lower()
        # Query example taken from 
        # https://github.com/alercebroker/usecases/blob/43e7775c00f6f949c5368aca3cd5dcc6bb64376c/notebooks/ALeRCE_Other_Watchlist.ipynb
        # Also, Python is too flexible... no error for slicing an empty string [2:]
        query = f"""
                WITH catalog ( ra, dec) AS (
                    VALUES
                        {objects_str}
                ),
                {prefix} (oid, classifier_name, class_name) AS (
                    SELECT
                        o.oid, p.classifier_name, p.class_name, p.probability, p.ranking
                    FROM
                        probability p
                    INNER JOIN 
                        object o
                    ON 
                        o.oid=p.oid
                    WHERE
                        p.classifier_name='{classifier}'
                        AND p.class_name IN ('{class_name}')
                        AND p.ranking=1 /* Ranking <==> Highest probability class name (AGN, SN, etc.) */
                        AND o.{mjd_choose} BETWEEN {mjd_date - days_backward} AND {mjd_date + days_forward}
                )

                SELECT 
                    o.oid, o.meanra, o.meandec, q3c_dist(c.ra,c.dec,o.meanra,o.meandec), 
                    o.{mjd_choose}, {prefix}.classifier_name, {prefix}.class_name, {prefix}.probability

                FROM object o INNER JOIN {prefix} ON {prefix}.oid=o.oid, catalog c
                    /*
                     * It is REALLY important to first use the catalog then the object ra,dec for speed. The radius is in degrees.
                     */
                WHERE
                    q3c_join(c.ra, c.dec, o.meanra, o.meandec, {radius})
                    AND o.{mjd_choose} BETWEEN {mjd_date - days_backward} AND {mjd_date + days_forward}
                ORDER BY
                    o.{order_by} {order_mode}
                """
        matches = matches.append(pd.read_sql(query, conn), ignore_index = True)
    
    conn.close()
    
    return matches.drop_duplicates(subset='oid')

# This function takes information in from the GW file and returns them in a dictionary
# I promise I'll make my functions properly documented... eventually
def read_gwfile(filepath: str, hdu_num = 1):
    
    properties = {}
    
    try:
        # Use astropy's unified file I/O
        hdu = Table.read(filepath, hdu = hdu_num)
        properties["mjd"] = hdu.meta["MJD-OBS"]

        properties["nest"] = True if hdu.meta["ORDERING"] == "NESTED" else False #save myself some time here
        properties["nside"] = hdu.meta["NSIDE"]
        properties["prob"] = hdu["PROB"].data
            
    except:
        filename = filepath.split("/")[-1]
        print("Could not open or use:", filename)
        print("In path:", filepath)
        #print("Trying the next file if it exists...")
        return properties 
    
    #if transient_candidate:
    #    targ_mjd = filepath.split("/")[-1].split("_")[-2] #to grab the date
    #    targ_mjd = targ_mjd[:4]+"-"+targ_mjd[4:6]+"-"+targ_mjd[6:] # Adding dashes for Time
    #    targ_mjd = Time(targ_mjd).mjd
    
    return properties #targ_ra, targ_dec, targ_mjd

# Borrowed and slightly modified from 'find_map_pixels.py', thanks Antonella!
def prob_pixel_locs(gw_in: dict, percentile = [0.9]):
    
    all_idx = []
    
    sort_percentile = sorted(percentile)
    max_percentile = sort_percentile[-1]
    #npix = len(gw_map)
    pb = gw_in["prob"] #hs.data['PROB']
    NSIDE = gw_in["nside"]
    idx_sort = np.argsort(pb)
    idx_sort_up = list(reversed(idx_sort))
    sum = 0.
    id = 0
    p_idx = 0
    
    while sum <= sort_percentile[-1]:
        
        this_idx = idx_sort_up[id]
        sum += pb[this_idx]
        id += 1
        
        if sum >= percentile[p_idx]:
            all_idx.append(idx_sort_up[:id])
            p_idx += 1
            
            if p_idx > len(percentile) - 1:
                break
        
    all_idx = dict(zip(sort_percentile, all_idx))
    total_area = {}
    for i in sort_percentile:
        area = hp.nside2pixarea(gw_in['nside'], degrees=True)*len(all_idx[i])
        print("The ", i*100. ,"% map is ",area," deg^2", sep = '')
        total_area[i] = area

    return all_idx, total_area


# ## Matching DESI observations to x% CI contour
# This merely matches skymap angles to tile pointings and indicates which tiles matched and their program (bright/dark).
# 
# *The code does not currently write anything to file*

# ## The two *major* matching functions, **initial_check** and **inner_matching**
# ## Initial check
# Performs the initial match with tile pointings. Grabs info from the exposures sql table, filters it, condenses by date, then sends all the good stuff to be matched via **inner_matching**. Finally outputs the results from that into a dictionary with the dates as keys and the sql table information as the values (sqlite3 row type, similar to a namedtuple/well-keyed dictionary). So there can be multiple elements as the values but only one key hence it's all organized by date. A convenient thing but also necessary for ALERCE efficiency.
# 
# e.g. {date:(exposure info)}

# From light_transient_matching, copy, pasted, and modified.

# Note, only matches 1-to-1. If you would like to find *all* matches in the area, set around = True and sit back.
# If you do set around = True, there are no data reduction methods in place yet so expect a lot of duplicates.

def initial_check(skymap_ra_in, skymap_dec_in, around = False) -> dict:
    
    # outputs a dictionary (desi_matches_dict) with keys as dates and values as sqlite3 rows containing the info grabbed from the exposures sql table
    # = {date:(exposure info)...}
    
    db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'

    # Per Antonella, no need to go further back
    query_date_start = "20201201"
    
    #today = Time.now()
    smushed_YMD = today.iso.split(" ")[0].replace("-","")
    
    query_date_end = smushed_YMD 

    # Handy queries for debugging/general information
    query2 = "PRAGMA table_info(exposures)"
    query3 = "PRAGMA table_info(tiles)"
    # Crossmatch across tiles and exposures to grab obsdate via tileid
    query_match = "SELECT distinct tilera, tiledec, obsdate, obsmjd, expid, program, exposures.tileid from exposures INNER JOIN tiles ON exposures.tileid = tiles.tileid where obsdate BETWEEN " +         query_date_start + " AND " + query_date_end + ";" #obsdate>20210228 

    # Querying sql and returning a data type called sqlite3 row, it's kind of like a namedtuple/dictionary
    conn = sqlite3.connect(db_filename)

    conn.row_factory = sqlite3.Row # https://docs.python.org/3/library/sqlite3.html#sqlite3.Row

    cur = conn.cursor()
    
    ''' ***** For handy purposes only *****
    
     cur.execute(query2)
     row2 = cur.fetchall()
     for i in row2:
         print(i[:])
    cur.execute(query)
    rows = cur.fetchall()
    
    *********************************** '''

    cur.execute(query_match)
    matches_list = cur.fetchall()
    cur.close()

    # I knew there was a way! THANK YOU!
    # https://stackoverflow.com/questions/11276473/append-to-a-dict-of-lists-with-a-dict-comprehension
    date_dict = {k['obsdate'] : list(filter(lambda x:x['obsdate'] == k['obsdate'], matches_list)) for k in matches_list}


    desi_matches_dict = {} # {i['obsdate']: [] for i in matches_list}
    all_confidence_matches = []
    
    # Uncomment the exp_ras/decs lines if retaining said coordinates is of interest (say for scatter plotting)
    # Or make a new function parameter... I quite like how clean it is now but up to you future user
    #all_exp_ras = []
    #all_exp_decs = []
    
    # Iterating day by day as a way to keep track
    for date, row in date_dict.items():
        #print(date)
        
        date_str = str(date)
        date_str = date_str[:4]+"-"+date_str[4:6]+"-"+date_str[6:] # Adding dashes for Time
        obs_mjd = Time(date_str).mjd

        # This method is *technically* safer than doing a double list comprehension with set albeit slower
        # The lists are small enough that speed shouldn't matter here
        # row has a separate element for all the tiles on that date
        unique_tileid = {i['tileid']:(i['tilera'], i['tiledec']) for i in row}
        exposure_ras, exposure_decs = zip(*unique_tileid.values())
        #all_exp_ras.extend(exposure_ras)
        #all_exp_decs.extend(exposure_decs)
        
        # **** Old method, deprecated ****
        #unique_ra_dec = list(set([(i[0], i[1]) for i in row])) # there's probably a way to do this in SQL... oh well
        #exposure_ras, exposure_decs = zip(*unique_ra_dec)
        # ********************************
        
        # Some renaming still to be done...
        #print(date)
        # Sends all the prepared information above to be matched and retains the coordinates of the matches
        desi_matches, confidence_matches = inner_matching(exposure_ras_in = exposure_ras, exposure_decs_in = exposure_decs, 
                                                              ra_in = skymap_ra_in, dec_in = skymap_dec_in, 
                                                              query_dates = query_date_start + query_date_end,
                                                              max_sep = 2, sep_units = 'deg', around = around)
        
        if desi_matches.size: # Interchangeable with confidence_matches
            # If any matches, print!
            print(date, '-', len(confidence_matches), 'match(es)')
            #print(len(confidence_matches), 'match(es)')
            if not around:
                all_confidence_matches.extend(confidence_matches) 
        else:
            continue

        desi_matches_dict[date] = []

        # Prepping output
        # Populating the dictionary by date (a common theme)
        # Each element in the dictionary thus contains the entire sqlite3 row (all info from sql tables with said headers)
        for tup in desi_matches:
            ra = tup.ra.deg
            dec = tup.dec.deg
            match_rows = [i for i in row if (i['tilera'], i['tiledec']) == (ra, dec)] # Just rebuilding for populating, this shouldn't change/exclude anything
            desi_matches_dict[date].extend(match_rows)
            
    return desi_matches_dict, all_confidence_matches #, all_exp_ras, all_exp_decs


# ## inner_matching
# #### AKA 'the bread & butter'
# 
# **inner_matching** does the actual matching via **match_coordinates_sky**. It outputs a SkyCoord array containing the 1-to-1 matches for DESI and whatever the second 'catalog' is. 
# 
# It handles a little bit of prep including removing NaNs but for the most part simply performs the match and spits it back.

# Perform the actual matching with a bit of pre-conditioning to ensure happiness for everyone involved (looking at you astropy)
def inner_matching(exposure_ras_in = np.array([]), exposure_decs_in = np.array([]), ra_in = np.array([]), dec_in = np.array([]), 
                   max_sep = 2, sep_units = 'deg', around = False, query_dates = ""):
    
    # hashtag units
    if sep_units == 'arcsec':
        max_sep *= u.arcsec
    elif sep_units == 'arcmin':
        max_sep *= u.arcmin
    elif sep_units == 'deg':
        max_sep *= u.deg
    else:
        print("Separation unit specified is invalid for matching. Defaulting to arcsecond.")
        max_sep *= u.arcsec
    
    # Just in case
    if not np.array(exposure_ras_in).size:
        print('Array of size 0 fed in. Returning no matches.')
        return np.array([]), np.array([])
    
    # Removing nans because they don't place nice in match_coordinates_sky
    nan_ra = np.isnan(exposure_ras_in)
    nan_dec = np.isnan(exposure_decs_in)
    
    if np.any(nan_ra) or np.any(nan_dec):
        print("NaNs found, removing them from array (not FITS) before match.")
        #print("Original length (ra, dec): ", len(target_ras), len(target_decs))
        nans = np.logical_not(np.logical_and(nan_ra, nan_dec))
        exposure_ras_in = exposure_ras_in[nans] # Logic masking, probably more efficient
        exposure_decs_in = exposure_decs_in[nans]
    
    # Sometimes this comes in handy but the code runs fast enough these days
    tree_name = "_".join(("kdtree", query_dates))
    
    desi_coords = SkyCoord(exposure_ras_in*u.deg, exposure_decs_in*u.deg)
    catalog_coords = SkyCoord(ra_in*u.deg, dec_in*u.deg)
    
    # When running against the targlist this goes crazy... not worth it but I'll keep it in for posterity
    # Lots of repeats and not yet filtered because it may not be necessary to ever use this
    if around: 
        idx_desi, idx_catalog, d2d, d3d = catalog_coords.search_around_sky(desi_coords, max_sep)
        
        desi_matches = desi_coords[idx_desi] 
        catalog_matches = catalog_coords[idx_catalog]
        
    else:
        idx_catalog, d2d, d3d = match_coordinates_sky(desi_coords, catalog_coords) #, storekdtree = tree_name) # store tree to speed up subsequent results

        sep_constraint = d2d < max_sep
        desi_matches = desi_coords[sep_constraint]
        catalog_matches = catalog_coords[idx_catalog[sep_constraint]]
    
    
    #if desi_matches.size:
    #    print(len(desi_matches), "match(es) found with separation -", max_sep)
        #print()

    return desi_matches, catalog_matches

# Function to read in fits table info, RA, DEC, MJD and targetid if so desired
# Uses control parameter tile to determine if opening tile exposure file or not since headers are different
def read_fits_info(filepath: str, tile = False):
    
    if not tile:
        hdu_num = 5
    else:
        hdu_num = 2 # For Zbest?
    
    try:
        with fits.open(filepath) as hdu1:
    
            data_table = Table(hdu1[hdu_num].data) #columns
        
            targ_ID = data_table['TARGETID']
            targ_ra = data_table['TARGET_RA'].data # Now it's a numpy array
            targ_dec = data_table['TARGET_DEC'].data
            
            # Not needed... kept for posterity
            # if not tile: 
            #     targ_mjd = hdu1[hdu_num].header['MJD-OBS'] # This is a string
            # else:
            #     targ_mjd = data_table['MJD'].data
            #     targ_mjd = Time(targ_mjd[0], format = 'mjd')
            
    except:
        filename = filepath.split("/")[-1]
        print("Could not open or use:", filename)
        #print("In path:", filepath)
        #print("Trying the next file...")
        return np.array([]), np.array([]), np.array([])
    
    return targ_ra, targ_dec, targ_ID

# Grabbing the frame fits files
def glob_frames(exp_d: str):   
    
    # This function grabs the names of all input files in the transient directory and does some python string manipulation
    # to grab the names of the input files with full path and the filenames themselves.
    try:
        filenames_read = glob(exp_d + "/cframe-" + color_band + "*.fits") # Only need one of b, r, z
        # sframes not flux calibrated
        # May want to use tiles... coadd (will need later, but not now)
    except:
        try:
            # Unsure if this is necessary but in looking in some directories, there are no cframes...
            # Put this here but the issue may have been somewhere else
            # TODO: Ask Antonella
            filenames_read = glob(exp_d + "/frame-" + color_band + "*.fits") # Only need one of b, r, z
        except:
            print("Could not grab/find any fits in the exposure directory:")
            print(exp_d)
            filenames_read = [] # Just in case
            #raise SystemExit("Exiting.")
            print("Continuing...")

    #else:
        #filenames_out = [s.split(".")[0] for s in filenames_read]
        #filenames_out = [s.split("/")[-1] for s in filenames_read]
        #filenames_out = [s.replace("in", "out") for s in filenames_out]
        
    return filenames_read #, filenames_out

#path_to_transient = "/global/cfs/cdirs/desi/science/td/daily-search/desitrip/out"
#print(all_candidate_filenames(path_to_transient)[1])

def pickle_dump(data_structure, filepath):
    
    with open(filepath, "wb") as f:    
        pickle.dump(data_structure, f)
    
    return None

def pickle_load(filepath):
    
    if not os.path.exists(filepath):
        print("Cannot find pickle file in")
        print(filepath)
        print("Continuing...")
        return None
    
    with open(filepath, "rb") as f:    
        data_structure = pickle.load(f)
    
    return data_structure

# ## Closer checking (handler) function
# In the algorithm, **initial_check** has already been run and has produced a dictionary of DESI matches and 1-to-1 targetlist matches. Using said dictionary, we only check the dates/tiles in which we found a match and then compare all 5000 fibers to the whole targetlist again to save on operation time and memory.
# 
# This function handles via **glob_frames** and **read_fits_info** FITS I/O and feeds that directly into **inner_matching** to perform the actual match. So technically it doesn't actually perform the matching but it's the handler. 
# 
# It ultimately outputs two lists, one with the DESI exposure matches and the other with matches to the second catalog. In this case, the first (DESI) isn't needed so that will be all but empty while the second has the necessary info.

def closer_check(matches_dict = {}, catalog2_ras = [], catalog2_decs = [], exclusion_list = []):
    
    # Initalizing output structures
    #all_exp_matches = {}
    all_targlist_matches = []
    
    # Just in case
    if not matches_dict:
        print("No far matches fed in for nearby matching. Returning none.")
        return {}
    
    # Iterate through dictionary from initial_check
    for date, row in matches_dict.items(): 
        print("\n", date)
        
        # If I have to run this multiple times, I probably already know what days there won't be a match so this is to save time
        if date in exclusion_list:
            print(f"{date} previously matched, continuing...")
            continue
        
        #all_exp_matches[date] = []
        file_indices = {}

        all_targ_ras = np.array([])
        all_targ_decs = np.array([])

        # Iterating tile by tile and grabbing all of the necessary info
        for i in row:
            exp_paths = '/'.join((exposure_path, "daily/exposures", str(i['obsdate']), f"{str(i['expid']):0>8}"))
            #print(exp_paths)
            
            # Now finally going file by file and putting all of those ras/decs into one giant numpy array
            # It should start making sense why we go date by date now... :)
            for path in glob_frames(exp_paths):
                #print(path)
                targ_ras, targ_decs, _ = read_fits_info(path, False)

                # For comparing to data tables in original file/sanity checking
                all_len = len(all_targ_ras)
                new_len = len(targ_ras)
                if all_len:
                    all_len -= 1
                    file_indices[path] = (all_len, all_len + new_len) # The start and end index, modulo number
                else:
                    file_indices[path] = (0, new_len) # The start and end index, modulo number

                if len(targ_ras) != len(targ_decs):
                    print("Length of all ras vs. all decs do not match.")
                    print("Something went wrong!")
                    print("Continuing but not adding those to match...")
                    continue

                all_targ_ras = np.append(all_targ_ras, targ_ras)
                all_targ_decs = np.append(all_targ_decs, targ_decs)

        # Here we finally perform the match, re-using inner_matching but with a much more strict matching radius
        # desi_dict_matches, targlist_fiber_matches 
        _ , targlist_fiber_matches = inner_matching(exposure_ras_in = all_targ_ras, exposure_decs_in = all_targ_decs,
                ra_in = catalog2_ras, dec_in = catalog2_decs, 
                max_sep = 1, sep_units = 'arcsec', around = False, query_dates = "")

        if targlist_fiber_matches:
            print(len(targlist_fiber_matches), "matches")
            all_targlist_matches.extend(targlist_fiber_matches)
        else:
            print("No matches found. Continuing...")
            continue
        
    #return all_exp_matches, all_targlist_matches
    return all_targlist_matches


# ## Building the DR9 targetlist 

# Borrowed from Segev's code, thanks Segev! 
# https://github.com/desihub/timedomain/blob/master/gwtarget/gw_dr9.ipynb

def build_targlist_table(nside, pix_map):
    hpdirnames = ['/global/project/projectdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright',
                  '/global/project/projectdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/dark']

    readcols = ['TARGETID', 'BRICKID', 'BRICKNAME', 'BRICK_OBJID',
                'RA', 'DEC', 'PMRA', 'PMDEC', 'REF_EPOCH',
                'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET',
                'FLUX_G', 'FLUX_R', 'FLUX_Z',
                'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z']

    targlist_threshold = None
    for hpdirname in hpdirnames:
        if targlist_threshold is None:
            targlist_threshold = Table(io.read_targets_in_hp(hpdirname, nside = nside, pixlist=pix_map, columns=readcols))
        else:
            targlist_threshold = vstack([targlist_threshold, Table(io.read_targets_in_hp(hpdirname, nside = nside, pixlist=pix_map, columns=readcols))])

    # targlist90.rename_column('BRICK_OBJID', 'OBJID')
    if not len(targlist_threshold): # i.e. 0 elements, I don't think Table has an 'empty' method
        return None
        
    targlist_threshold = unique(targlist_threshold)
    
    return targlist_threshold


def targlist_write(table_in, filename, overwrite = True):
    # building the table takes *forever*
    # if overwrite = True, then overwrites regardless
    # if false, checks if the tables are equivalent and does nothing if they are
    # otherwise appends a 1 to the filename and writes
    # So if this were to happen a bunch you'd have xyz1111111.ecsv ;)
    
    if not overwrite:
        try:
            table_read = Table.read(filename)
        except:
            print("Can't find table specified per", filename)
            print("Writing a new one to", filename)
            overwrite = True

    if overwrite:
        table_in.write(filename, overwrite = True)
    else:
        if setdiff(table_in, table_read): # If there are differences
            name, ext = filename.split('.')
            name += "1"
            table_in.write(".".join((name, ext)))
    
    return None

# Operates in a similar fashion to closer_check, just with targetids
# Also uses setdiff because I like it and I can iteratively reduce the original targetlist with little fuss and have it done
def grab_desi_targetid(matches_dict = {}, targlist_table = Table(), exclusion_list = []):
    
    if not targlist_table:
        print("No targetlist table fed in. Exiting and returning empty table.")
        return targlist_table
    
    if not matches_dict:
        print("No far matches fed in for targetid matching. Returning an empty table.")
        return Table()
    
    # For calling setdiff
    # I don't overwrite targlist_table... just in case... out of habit... because I'm scared
    table_reduced_targids = targlist_table
    
    for date, row in matches_dict.items(): 
        
        print(date)
        if date in exclusion_list:
            continue
            
        target_ids_date = []

        for i in row:
            exp_paths = '/'.join((exposure_path, "daily/exposures", str(i['obsdate']), f"{str(i['expid']):0>8}"))
            #print(exp_paths)
            #all_exp_fits[date].extend()
            for path in glob_frames(exp_paths):
                #print(path)
                _, _, target_ids = read_fits_info(path, False)

                target_ids_date.extend(target_ids)

        # Convert targetids to astropy table
        target_ids_date = np.array(target_ids_date, dtype = 'longlong')
        DESI_targ_id_table = Table([target_ids_date], names = ['TARGETID'])
        #targetid_matches = setdiff(targlist_table['TARGETID'], DESI_targ_id_table)
        
        # Perform setdiff
        reduced_targids = setdiff(table_reduced_targids, DESI_targ_id_table, keys = ['TARGETID'])
        if reduced_targids: # Since setdiff will return an empty table if there are no matches... which would overwrite everything
            table_reduced_targids = reduced_targids 
        else:
            print("No targetid matches found on:", date)
            
    return table_reduced_targids


# ## Matching CI interval pixel locations to DESI tile pointings

def plot_cartmap_tiles(lvc_healpix_file, levels=[0.5, 0.9], angsize=3., tile_ra=None, tile_dec=None, targ_ra=None, targ_dec=None, program_names = None):
    """Plot the GW map with the DESI footprint in a Cartesian projection.
    
    Parameters
    ----------
    lvc_healpix_file : str
        Relative or absolute path to LIGO/Virgo HEALPix angular reconstruction file.
    levels : list
        List of credible interval thresholds, e.g., 0.5, 0.9, etc.
    angsize : float
        Size of plot (-angsize, +angsize) in degrees about the center.
    tile_ra : list or ndarray
        List of RAs for DESI tiles (in deg).
    tile_dec : list or ndarray
        List of declinations for DESI tiles (in deg).
    targ_ra : list or ndarray
        List of RAs for DESI targets (in deg).
    targ_dec : list or ndarray
        List of declinations for DESI targets (in deg).
    
    Returns
    -------
    fig : matplotlib.Figure
        Figure object for accessing or saving a plot.
    """
    gwmap = hp.read_map(lvc_healpix_file)
    npix = len(gwmap)
    nside = hp.npix2nside(npix)

    # Compute contours.
    if nside > 256:
        _gwmap = hp.pixelfunc.ud_grade(gwmap, 256)
        _gwmap = _gwmap / np.sum(_gwmap)
    else:
        _gwmap = gwmap
    ra_contour, dec_contour = compute_contours(levels, _gwmap)

    # Create a temporary plot to produce a nice image array.
    # This code sets the size of the map around the maximum value.
    maxpix = np.argmax(gwmap)
    ra_c, dec_c = hp.pix2ang(nside, maxpix, lonlat=True)

    xmin = np.round(ra_c - angsize)
    xmax = np.round(ra_c + angsize)
    if xmax < xmin:
        xmin, xmax = xmax, xmin
    cxmin, cxmax = xmin, xmax
    frot = 0.
    if xmax > 90 and xmax < -90:
        frot, cxmin, cmax = 180., xmax-180., xmax+180.
    ymin = np.round(dec_c - angsize)
    ymax = np.round(dec_c + angsize)

    faspect = np.abs(cxmax - cxmin)/np.abs(ymax-ymin)
    fysize = 4
    figsize = (fysize*faspect+1, fysize+2)

    # Open and close the temporary plot.
    tfig   = plt.figure(num=2,figsize=figsize)
    rotimg = hp.cartview(gwmap, fig=2,coord='C', title="", cbar=False, flip='astro',
                         lonra=[cxmin,cxmax], latra=[ymin,ymax], rot=frot,
                         notext=True, xsize=1000,
                         return_projected_map=True)
    plt.close(tfig)

    # Now make the real plot with the desired angular contours.
    fig, ax = plt.subplots(1,1, num=1, figsize=figsize)
    img = ax.imshow(rotimg, extent=[cxmax, cxmin, ymin, ymax],
                    origin='lower', cmap='OrRd')

    for i, (rc, dc, lstyle, clev) in enumerate(zip(ra_contour, dec_contour, ['--', '-'], ['50', '90'])):
        p = ax.plot(rc, dc, 'g-', ls=lstyle, lw=2, label='{}% CI'.format(clev))

    ax.set(xlim=(cxmax, cxmin),
       xlabel='RA [deg]',
       ylabel='Dec [deg]')
    ax.grid(ls=':')
    
    bright_list = ['BRIGHT', 'BGS']
    dark_list = ['DARK', 'ELG', 'QSO', 'LRG']
    # Add DESI targets, specified by RA, Dec.
    if targ_ra is not None and targ_dec is not None:
        if program_names:
            
            # if bright in program_info or if '---'.lower() == bgs, then plot as bright
            # if dark in ... or if ''.lower() == elg or qso or lrg
            
            bright_pointing = np.array([(i,j) for i,j,k in zip(targ_ra, targ_dec, program_names) if any(substring in k.upper() for substring in bright_list)])
            dark_pointing = np.array([(i,j) for i,j,k in zip(targ_ra, targ_dec, program_names) if any(substring in k.upper() for substring in dark_list)])
            
            bright_dots = ax.plot(bright_pointing[:, 0], bright_pointing[:, 1], 'c.', alpha = 1, label = 'BRIGHT', marker = "+") 
            dark_dots = ax.plot(dark_pointing[:, 0], dark_pointing[:, 1], 'm.', alpha = 0.7, label = 'DARK', marker = "x") 
            
        else:
            ax.plot(targ_ra, targ_dec, 'k.', alpha=0.5, label = 'Matches') # temp change, alpha = 0.1 -> alpha = 0.5 (maybe push command line arg for this)

    _h, _l = ax.get_legend_handles_labels()

    # Add DESI tile drawings, specified by central RA, Dec.
    if tile_ra is not None and tile_dec is not None:
        #for _ra_c, _dec_c in zip(tile_ra, tile_dec):
        for _ra_c, _dec_c in zip(targ_ra, targ_dec):
            circ = plt.Circle((_ra_c, _dec_c), radius=1.6, fc='None', ec='b', ls=':', lw=2)
            ax.add_artist(circ)

        _h.append(circ)
        _l.append('DESI FOV')
    
    ax.legend(handles=_h, labels=_l, fontsize=10, ncol=2)

    cb = fig.colorbar(img, orientation='horizontal', shrink=0.95,
                      fraction=0.04, pad=0.2, ax=ax)
    cb.set_label(r'$dp/d\Omega$ [deg$^{-2}$]')

    return fig

def mainInjector_env(path = "./"):
    
    _ = [sys.path.append(i.path) for i in os.scandir(os.path.join(user_home, "timedomain/gwtarget/DESI_mainInjector/Main-Injector-master")) if i.is_dir()]
    sys.path.append(os.path.join(user_home, "timedomain/gwtarget/DESI_mainInjector"))
    
    os.environ["DESGW_DIR"] = os.path.join(path, "../Main-Injector-master/python/")
    os.environ["DESGW_DATA_DIR"] = os.path.join(path, "../Main-Injector-master/data/")
    if os.environ["DESGW_DIR"] not in os.environ["PYTHONPATH"]:
        os.environ["PYTHONPATH"] = os.environ["PYTHONPATH"] + ":" + os.environ["DESGW_DIR"]
     
    # My, my, what convention to use here?
    #return

def disruptive_mode(gwfile_path = "./", gw_name = ""):
    # TODO: What to do with maininjector as a whole? That's a lot of files to copy!
    # Running maininjector
    # user_home = "/global/u2/p/portmanm/"
    
    # TODO: Set this as an environment variable? Or grab it from the environment?
    path_to_MI = os.path.join(user_home, "timedomain/gwtarget/DESI_mainInjector/work_dark/")
    
    mainInjector_env(path_to_MI)
    
    import resimulator
    
    try: 
        from pyslalib import slalib
    except:
        _ = subprocess.run(sys.executable, '-m', 'pip', 'install', '--user', 'pyslalib', shell = True, check = True)
        
    if not gw_name:
        gw_name = os.path.basename(gwfile_path).split(".fits")[0]
        
    current_dir = os.path.join(os.getcwd(), "./") # To make sure we end in a /
        
    os.chdir(path_to_MI)
    #resimulator.recycle(gw_name, gwfile_path, "dark", path_to_MI, do_make_maps = True, do_make_hexes = True, do_make_jsons = False, do_make_gifs = False)
    # Default to keeping maps but for testing it's false
    # In future, delete after 1 month or something
    resimulator.recycle(gw_name, gwfile_path, "dark", current_dir, do_make_maps = False, do_make_hexes = True, do_make_jsons = False, do_make_gifs = False)
    os.chdir(current_dir)
    
    ra_dec_in = os.path.join(current_dir, gw_name + "-ra-dec-id-prob-mjd-slot-dist.txt")
    ra_dec_out = os.path.join(current_dir, gw_name + "_disruptive_ToO_ledger_pointings.txt")
    
    with open(ra_dec_in, "r") as i:
        file_text = i.readlines()
        
    with open(ra_dec_out, "w") as o:
        header = "name," + file_text[0].strip("#")
        new_header = header.split(" ")
        # It should be the same every time but just in case...
        new_header = " ".join(["mjd_start, mjd_end," if "mjd" in head else head for head in new_header])
        
        o.write(new_header)
        
        for line in file_text[1:]:
            # We use .replace For an errant extra space in the test output
            values = line.replace("  "," ").split(" ")
            mjd_end = str(float(values[4]) + 10.) #mjd
            values.insert(0, gw_name)
            values.insert(6, mjd_end)
            o.write(", ".join(values))
            
    # Can I split this off and run maininjector *while* running the rest of the code to generate the below?
    
    return ra_dec_out

def recursive_pix_filter(map_properties, best_tile_pix, tile_dict, r_num, best_ids, threshold):
    # This function performs a recursive search of the tiles
    # that cover the CI region and chooses tiles which cover
    # the most probability while excluding the pixels that overlap.
    # 
    # This is mostly deprecated but left in just in case an area
    # has a great deal of coverage in one of the programs.    
    # 
    # This method is no longer prioritized so the following TODO
    # may take awhile...
    # TODO: Make this work with full astropy TABLE
    
    if r_num >= len(tile_dict):
        return best_ids
    
    if len(best_ids) == threshold:
        return best_ids
    
    print('Recursion', r_num)
    
    #pix_in_tiles = np.zeros(np.shape(best_tile_pix))
    
    all_probs = np.zeros((r_num,2))
    
    # pass in prob, best pix
    for ID in np.setdiff1d(list(tile_dict.keys()), best_ids):
        remain = np.isin(tile_dict[ID], best_tile_pix, assume_unique = True, invert = True)
        pix_remain = tile_dict[ID][remain]
        
        if pix_remain.size == 0:
            continue
        
        all_probs = np.vstack((all_probs, (ID, np.sum(np.array(map_properties["prob"])[pix_remain]))))
        
    # Grabs argmax of prob column of all probs, spits back ID
    new_max_id = all_probs[np.argmax(np.array(all_probs)[:,1]),0]
    
    if new_max_id in best_ids:
        print(f"Something went wrong! ID {new_max_id} already exists in best.")
        _ = [print(int(i[0]), i[1]) for i in all_probs]
    elif not new_max_id:
       #print(all_probs)
        print("No non-overlapping tiles found in chosen CI threshold.")
        print("Returning rest of tiles in area (up to number of pointings) as ranked by probability.")
        #print("Either turn off CI restriction or allow overlap.")
        return best_ids
    
    r_num += 1
    
    old_len = len(best_tile_pix)
    best_tile_pix = np.union1d(best_tile_pix, tile_dict[new_max_id])
    
    if not len(best_tile_pix) - old_len:
        # Grab second to last
        new_max_id = all_probs[np.argsort(np.array(all_probs)[:,1]),0][-2]
    
    best_ids = np.append(best_ids, new_max_id)
    
    best_ids = recursive_pix_filter(map_properties, 
                                    best_tile_pix, 
                                    tile_dict,
                                    r_num, 
                                    best_ids,
                                    threshold)
    
    return best_ids

def nondisruptive_mode(map_properties: dict, degrade_map_properties: dict, pixmap, 
                       num_pointings:int = 10,  CI_level = 0.9,
                       restrict = False, overlap = True):
    
    # This ND mode operates by calculating the probability covered by tiles in the 
    # 90% CI region of the map. This is slightly slower than interpolating over the probability
    # of the pixels in each tile although the timing is fairly negligible. 
    # This way is less cool but the results seem better grouped -- Matt P.
    #
    # The interpolation method can be found in the jupyter notebook of the same name.
    # EXAMPLE CALL:
    #     result_dict = nondisruptive_mode(gw_properties, gw_degraded_properties, pixmap = pixmap, num_pointings = 15)
    # 
    # I DO NOT RECOMMEND RUNNING THIS WITH RESTRICT ON AND OVERLAP OFF
    #
    # Restrict limits to pixmap, which in here represents the 90% CI so there has to be 
    # *SOME* overlap with that exact boundary
    #
    # Overlap limits whether or not tiles can overlap. 
    # Defaults to True since the overlaps are likely to be 
    # found in different passes and programs
    
    pixmap = np.array(pixmap)
    
    # Grab ecsv with tile numbers, positions, and # of observations
    tiles_path = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv"
    tile_info = unique(Table.read(tiles_path), keys=['RA', 'DEC'])
    max_pass_num = np.max(tile_info['PASS'])
    # print(tile_info.columns)
    
    # PRIORITIZE LOWER PASS
    # MAX SHOULD BE PER PROGRAM, num_pointings *per* program (bright, dark, backup[may remove] -- we think)
    
    conditions = tile_info['IN_DESI'] & (tile_info['STATUS'] == 'unobs') #& (tile_info['PROGRAM'] != 'BACKUP')
    tile_info = tile_info[conditions] # filter out by IN_DESI
    programs = ['DARK', 'BRIGHT', 'BACKUP']
    
    best_tiles = {}
    
    for p in programs:
        
        in_DESI = tile_info[tile_info['PROGRAM'] == p]
    
        tile_ra = np.array(in_DESI['RA'])
        tile_dec =  np.array(in_DESI['DEC'])

        tile_skycoord = SkyCoord(tile_ra*u.deg, tile_dec*u.deg)
    
        # grab top 20
        prob_sum = 0
        count = 0
        idxs = []
        d_prob = degrade_map_properties["prob"]
        arg_prob = np.argsort(d_prob)[::-1]

        while prob_sum < CI_level:
            prob_sum += d_prob[arg_prob[count]]
            count += 1
            idxs.append(arg_prob[count])
        
        ra_d, dec_d = hp.pix2ang(degrade_map_properties["nside"], idxs, nest = degrade_map_properties["nest"], lonlat = True) 
        pixmap_skycoord = SkyCoord(ra_d*u.deg, dec_d*u.deg)

        # I bet this is faster than conesearch and checking matches of indices
        # Don't know the order but this should be better optimized than a loop
        # ... unless we parallelize but it's not worth it here
        _, d2d_tile, _ = match_coordinates_sky(tile_skycoord, pixmap_skycoord)
        # Filtering by maximum separation and closest match
        sep_constraint = d2d_tile < 4*u.deg # Extended, doubled tile radius
        tile_matches = tile_skycoord[sep_constraint]
        in_constraint = in_DESI[sep_constraint]
            
        # Initialize astropy HEALPix
        # SkyCoord defaults to ICRS so we hardcode it here

        astro_hp = ah.HEALPix(nside = map_properties["nside"], order = "nested" if map_properties["nest"] else "ring", frame = 'icrs')

        prob_vals = []
        # Calculate best match
        for coord in tile_matches:
            hp_idx = astro_hp.cone_search_skycoord(coord, tile_rad)
            both = np.intersect1d(hp_idx, pixmap, assume_unique = True)

            if both.size:
                #count += 1
                prob_vals.append(np.sum(map_properties["prob"][both]))
            else:
                prob_vals.append(0)

        prob_vals = np.array(prob_vals)
        #args_prob_sum = np.argsort(prob_vals)[-num_pointings:]
        
        # Named to match write_too_ledger
        in_constraint['PROB_COVERED'] = prob_vals*100
        
        # Sort by pass *and then* by probability
        in_constraint.sort(['PASS','PROB_COVERED'])
        
        # Reverse the order so that highest probability is at the top
        in_constraint.reverse()
        in_constraint = in_constraint[in_constraint['PROB_COVERED'] > 0]
        
        tile_table = in_constraint[in_constraint['PASS'] == 0]
        pass_num = 1
        
        while (len(tile_table) < num_pointings) and (pass_num < max_pass_num):
            in_pass = in_constraint[in_constraint['PASS'] == pass_num]
            tile_table = vstack([tile_table, in_pass])
            pass_num += 1
            
        # Limit to num_pointings in case/when we go over the amount needed
        best_tiles[p] = tile_table[:num_pointings]
        
        if len(tile_table) < num_pointings:
            print(f"Could not find requested number of tiles ({num_pointings}) in/near 90% CI for {p} program.")
            print(f"Found {len(tile_table)} tiles in/near 90% CI covering a total probability of {np.sum(tile_table['PROB_COVERED']):.4f}%.\n")
        
        # Go to next program... avoid a tab cluster on the way ;) 
        if overlap or not len(tile_table):
            continue
            
        args_interp = np.argsort(in_constraint['PROB_COVERED'])[::-1]

        tileid_red = in_constraint['TILEID'][args_interp]
        
        tile_ra_c = np.array(in_constraint['RA'])[args_interp]
        tile_dec_c =  np.array(in_constraint['DEC'])[args_interp]
        
        skycoord_red = SkyCoord(tile_ra_c*u.deg, tile_dec_c*u.deg)

        if restrict:
            tile_pix = {ID:np.intersect1d(astro_hp.cone_search_skycoord(coord, tile_rad), pixmap) 
                        for ID, coord in zip(tileid_red, skycoord_red)}
        else:
            tile_pix = {ID:astro_hp.cone_search_skycoord(coord, tile_rad) for ID, coord in zip(tileid_red, skycoord_red)}

        tile_coord = {ID:coord for ID, coord in zip(tileid_red, skycoord_red)}

        # initializing with the best tile
        # return array will have size of pixmap

        # Dictionaries are explicitly ordered in >Python3.6 so this algorithm relies on that
        pix_in_best_tile = pixmap[np.isin(pixmap, tile_pix[tileid_red[0]], assume_unique = True)]

        best_ids = np.array([tileid_red[0]], dtype = int) #, tileid_red[new_max]], dtype = int)


        best_ids = recursive_pix_filter(map_properties, 
                             pix_in_best_tile, 
                             tile_pix,
                             1, 
                             best_ids,
                             num_pointings)

        if restrict and len(best_ids) != num_pointings:
            tileid_redd = np.setdiff1d(tileid_red, best_ids)
            best_ids = np.append(best_ids, tileid_redd)[:num_pointings]
        
        best_tiles[p] = best_tiles[p][np.isin(best_tiles[p]['TILEID'], best_ids, assume_unique = True)]
    
    return best_tiles

def plot_pointings(ra_map, dec_map, pointing_dict:dict, savename:str):
    
    plt.scatter(ra_map, dec_map, alpha=0.2)
    for p in ('DARK', 'BRIGHT', 'BACKUP'):
        plt.scatter(pointing_dict[p]['RA'], pointing_dict[p]['DEC'], label=p, alpha = 0.8, marker="+")
        
    plt.legend()
    plt.savefig(savename)
    return None

if __name__ == "__main__":
    pass
    # Unit tests????