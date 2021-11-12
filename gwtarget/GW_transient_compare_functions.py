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
# For importing useful functions
sys.path.append('/global/homes/p/portmanm/timedomain/')
_ = [sys.path.append(user_home + '/desi/' + x + '/py/') for x in os.listdir('/global/homes/p/portmanm/desi/')]

from astropy.io import fits
from astropy.table import Table, Column, join, hstack, vstack, unique, setdiff
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, Angle
from astropy.time import Time

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

import sqlite3


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
    mode = 'w' if overwrite else 'a'
    if verbose:
        mode = mode + '+'
    
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
    # meta: {DEPNAM00: desitarget, DEPNAM01: desitarget-git, DEPVER00: 0.53.0.dev4635, DEPVER01: 0.53.0-24-g58c9a719, EXTNAME: TOO, RELEASE: 9999}
    # schema: astropy-2.0
    RA DEC PMRA PMDEC REF_EPOCH CHECKER TOO_TYPE TOO_PRIO OCLAYER MJD_BEGIN MJD_END TOOID\n""")
            
        today = Time.now()
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
def access_alerts(radius = 3600*4, order_by = 'oid', order_mode = 'DESC', classifier = 'stamp_classifier', class_names=['SN', 'AGN'], **kwargs):
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

# This function takes information in from the GW file and returns them in a dictionary
# I promise I'll make my functions properly documented... eventually
def read_gwfile(filepath: str, hdu_num = 1):
    
    properties = {}
    
    try:
        with fits.open(filepath) as hdu1:
    
            hdr = hdu1[hdu_num].header
        
            properties["mjd"] = hdr["MJD-OBS"]
            properties["nside"] = hdr["NSIDE"]
            properties["nest"] = True if hdr["ORDERING"] == "NESTED" else False #save myself some time here
            #print(hdr['ORDERING'])
            
            data_table = Table(hdu1[hdu_num].data) #columns
            properties["prob"] = data_table["PROB"].data
        
            #targ_id = data_table['TARGETID']
            #targ_ra = data_table['TARGET_RA'].data # Now it's a numpy array
            #targ_dec = data_table['TARGET_DEC'].data
            #targ_mjd = data_table['MJD'][0] some have different versions of this so this is a *bad* idea... at least now I know the try except works!
            
    except:
        filename = filepath.split("/")[-1]
        print("Could not open or use:", filename)
        print("In path:", filepath)
        print("Trying the next file if it exists...")
        return properties #np.array([]), np.array([]), np.array([])
    
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
    
    for i in range(len(sort_percentile)):
        area = hp.nside2pixarea(gw_in['nside'], degrees=True)*len(all_idx[i])
        print("The ", percentile[i]*100. ,"% map is ",area," deg^2", sep = '')

    return all_idx


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

def read_fits_info(filepath: str, tile = True):
    
    # Disabling INFO logging temporarily to suppress INFO level output/print from read_spectra
    logging.disable(logging.INFO)
    
    try:
        spec_info = read_spectra(filepath).fibermap
        
        targ_id = spec_info['TARGETID'].data
        targ_ra = spec_info['TARGET_RA'].data # Now it's a numpy array
        targ_dec = spec_info['TARGET_DEC'].data
        targ_mjd = spec_info['LAST_MJD'] #.data
        
        targ_mjd = Time(targ_mjd[0], format = 'mjd')
        #targ_mjd = spec_info['FIRST_MJD'].data[0]
            
    except:
        filename = filepath.split("/")[-1]
        print("Could not open or use:", filename)
        #print("In path:", filepath)
        #print("Trying the next file...")
        return np.array([]), np.array([]), 0, 0
    
    if not tile and not np.all(targ_mjd):
        print("Unable to grab mjd from spectra, taking it from the filename...")
        targ_mjd = filepath.split("/")[-1].split("_")[-2] #to grab the date
        #targ_mjd = targ_mjd[:4]+"-"+targ_mjd[4:6]+"-"+targ_mjd[6:] # Adding dashes for Time
        targ_mjd = Time(targ_mjd, format = 'mjd') #.mjd
        
    # Re-enabling logging for future calls if necessary
    logging.disable(logging.NOTSET)    
    
    return targ_ra, targ_dec, targ_mjd, targ_id

# Grabbing the frame fits files
def glob_frames(exp_d: str):   
    
    # This function grabs the names of all input files in the transient directory and does some python string manipulation
    # to grab the names of the input files with full path and the filenames themselves.

    try:
        filenames_read = glob(exp_d + "/cframe-" + color_band + "*.fits") # Only need one of b, r, z
        # sframes not flux calibrated
        # May want to use tiles... coadd (will need later, but not now)
    
    except:
        print("Could not grab/find any fits in the exposure directory:")
        print(exp_d)
        filenames_read = [] # Just in case
        #filenames_out = [] # Just in case
        raise SystemExit("Exitting.")
        
    #else:
        #filenames_out = [s.split(".")[0] for s in filenames_read]
        #filenames_out = [s.split("/")[-1] for s in filenames_read]
        #filenames_out = [s.replace("in", "out") for s in filenames_out]
        
    return filenames_read #, filenames_out

#path_to_transient = "/global/cfs/cdirs/desi/science/td/daily-search/desitrip/out"
#print(all_candidate_filenames(path_to_transient)[1])


# ## Closer checking (handler) function
# In the algorithm, **initial_check** has already been run and has produced a dictionary of DESI matches and 1-to-1 targetlist matches. Using said dictionary, we only check the dates/tiles in which we found a match and then compare all 5000 fibers to the whole targetlist again to save on operation time and memory.
# 
# This function handles via **glob_frames** and **read_fits_info** FITS I/O and feeds that directly into **inner_matching** to perform the actual match. So technically it doesn't actually perform the matching but it's the handler. 
# 
# It ultimately outputs two lists, one with the DESI exposure matches and the other with matches to the second catalog. In this case, the first (DESI) isn't needed so that will be all but empty while the second has the necessary info.

def closer_check(matches_dict = {}, catalog2_ras = [], catalog2_decs = [], exclusion_list = []):
    
    # Initalizing output structures
    all_exp_matches = {}
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
            continue

        all_exp_matches[date] = []
        file_indices = {}

        all_targ_ras = np.array([])
        all_targ_decs = np.array([])

        # Iterating tile by tile and grabbing all of the necessary info
        for i in row:
            exp_paths = '/'.join((exposure_path, "daily/exposures", str(i['obsdate']), "000"+str(i['expid'])))
            #print(exp_paths)
            
            # Now finally going file by file and putting all of those ras/decs into one giant numpy array
            # It should start making sense why we go date by date now... :)
            for path in glob_frames(exp_paths):
                #print(path)
                targ_ras, targ_decs, _, _ = read_fits_info(path, True)

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
        desi_dict_matches, targlist_fiber_matches = inner_matching(exposure_ras_in = all_targ_ras, exposure_decs_in = all_targ_decs,
                ra_in = catalog2_ras, dec_in = catalog2_decs, 
                max_sep = 1, sep_units = 'arcsec', around = False, query_dates = "")

        if targlist_fiber_matches:
            print(len(targlist_fiber_matches), "matches")
            all_targlist_matches.extend(targlist_fiber_matches)
        else:
            print("No matches found. Continuing...")
            continue
        
    return all_exp_matches, all_targlist_matches


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
            targlist_threshold = vstack(targlist_threshold, Table(io.read_targets_in_hp(hpdirname, nside = nside, pixlist=pix_map, columns=readcols)))

    # targlist90.rename_column('BRICK_OBJID', 'OBJID')
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
        print("No targetlist table fed in. Exitting and returning empty table.")
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
            exp_paths = '/'.join((exposure_path, "daily/exposures", str(i['obsdate']), "000"+str(i['expid'])))
            #print(exp_paths)
            #all_exp_fits[date].extend()
            for path in glob_frames(exp_paths):
                #print(path)
                _, _, _, target_ids = read_fits_info(path, True)

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

# Unit tests????