#!/usr/bin/env python
# coding: utf-8

# Author: Matthew Portman
# Date: 11/11/21

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
#sys.path.append('/global/homes/p/portmanm/desi/desiutil/py/')

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

# Importing all functions here
from GW_transient_compare_functions import *

if __name__ == "__main__":
# *************************************************************************
#
# SETTING SOME DEFAULT VARIABLES AND PREPARING
#
# *************************************************************************
    
    global db_filename
    db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'

    # Hardcoding this for now until I can get the environment set up right
    global exposure_path
    exposure_path = "/global/cfs/cdirs/desi/spectro/redux" #os.environ["DESI_SPECTRO_REDUX"]

    global color_band
    color_band = "r"
    
    global today
    today = Time.now()
    
    targetlists_path = os.path.join(user_home, 'targetlists')
    # Check if targetlists folder exists, if not make one
    if not os.path.isdir(targetlists_path):
        os.makedir(targetlists_path)

# *************************************************************************
#
# SETTING UP AND GRABBING COMMAND LINE ARGS
#
# *************************************************************************
    
    help_text = sys.argv[0] + " handles and compares Gravitational Wave (GW) " \
    "localization maps to DESI observations. Most generally, the code finds those " \
    "observations in the area of a user-specified confidence interval (CI) of the GW map."
    
    parser = argparse.ArgumentParser(description = help_text)
    c_intervals = parser.add_mutually_exclusive_group()
    modes = parser.add_mutually_exclusive_group()

    
    parser.add_argument(dest     = 'gwfile', 
                        type     = str, 
                        help     = 'The input GW file')
    
    c_intervals.add_argument('--confidence_interval=90', '-CI90', 
                        dest      = 'CI_val',
                        action    = 'store_const', 
                        const     = [0.9],
                        help      = 'Specifying 90% CI')
    
    c_intervals.add_argument('--confidence_interval=95', '-CI95', 
                        dest      = 'CI_val',
                        action    = 'store_const', 
                        const     = [0.95],
                        help      = 'Specifying 95% CI')
    
    modes.add_argument('--disruptive', '-D', 
                        dest   = 'mode', 
                        action = 'store_true',
                        help   = 'Enable disruptive mode.')
    
    modes.add_argument('--non-disruptive', '-ND', 
                        dest   = 'mode', 
                        action = 'store_false',
                        help   = 'Enable non-disruptive mode.')

    args = parser.parse_args() # Using vars(args) will call produce the args as a dict
    gwfile = args.gwfile
    
    try:
        CI_val = args.CI_val
    except:
        CI_val = [0.9]
        
    try:
        mode = args.mode
    except:
        mode = False

# *************************************************************************
#
# READING IN SKYMAP, GRABBING INFORMATION, AND CONVERTING TO PIXEL MAP
#
# *************************************************************************
    h = fits.open(gwfile)
    #h=fits.open('skymaps/GW190412_combined_skymap.fits.gz')
    head = h[1].header
    gw_mjd = head['MJD-OBS']

    # Read in GW file, grab its properties, and determine pixels in CI area
    gw_name = gwfile.split('_')[0].split('/')[1]
    gw_properties = read_gwfile(gwfile)
    gw_map = hp.read_map(gwfile, nest = gw_properties["nest"])
    
    # Grab pixel locations for probabilities in x% CI
    pixmap = prob_pixel_locs(gw_properties, percentile = CI_val)
    
    # Converting skymap to pixelmap
    ra_map, dec_map = hp.pix2ang(gw_properties["nside"], pixmap, nest = gw_properties["nest"], lonlat = True)

# *************************************************************************
#
# DEGRADING PIXEL MAP TO CONSTRAIN PIXEL SEARCH SPACE
#
# *************************************************************************
    # Finds Supernova (SN) and Active Galactic Nuclei (AGN) matches in ALERCE alert broker data and matches to x% probability angles from above
    gw_degraded_properties = deepcopy(gw_properties)
    gw_degraded_properties["nside"] = 32
    print("Degrading...")

    # By using power = -2, we keep the sum of the map invariant (per documentation)
    degrade_map = hp.ud_grade(gw_map, nside_out = gw_degraded_properties["nside"], order_in = gw_properties["nest"], order_out = gw_properties["nest"], power = -2)

    gw_degraded_properties["prob"] = degrade_map
    pix_degraded = prob_pixel_locs(gw_degraded_properties, percentile = CI_val)

    ra_degraded, dec_degraded = hp.pix2ang(gw_degraded_properties["nside"], pix_degraded, nest = gw_degraded_properties["nest"], lonlat = True)
    
# *************************************************************************
#
# GRABBING TRANSIENT ALERT MATCHES 'WITHIN' DEGRADED PIXEL MAP PIXELS
#
# *************************************************************************
# mean pixel spacing 1.8323 deg, approximating to 2 deg or 7200 arcseconds (radius)
# We expand this to 4 degrees to account for edge cases/diameter
# Also doing it this way allows for non-continuous CI contours

# As a reminder, we rely on their classifier to determine SN or AGN

    rad = 3600*4
    days_forward = 30

    try:
        agn_old_alerts = access_alerts(order_by = 'firstmjd',
                                       order_mode = 'ASC',
                                       class_names = ['AGN'],
                                       firstmjd = [gw_mjd, gw_mjd + days_forward],
                                       ra = ra_degraded[0],
                                       dec = dec_degraded[0],
                                       radius = rad
                                       )
        #mask=(agn_old_alerts['firstmjd']<gw_mjd+30)
        #alerts_agn_2019=agn_old_alerts[mask]
    except:
        agn_old_alerts = pd.DataFrame()

    # I use a try because if there are no matches, it fails then everything falls apart
    try: 
        alerts_sn = access_alerts(class_names = ['SN'],
                                  order_by = 'firstmjd',
                                  order_mode = 'ASC',
                                  firstmjd = [gw_mjd, gw_mjd + days_forward],
                                  ra = ra_degraded[0],
                                  dec = dec_degraded[0],
                                  radius = rad
                                  )

        #mask=(alerts_sn['firstmjd']<gw_mjd+30)
        #alerts_sn=alerts_sn[mask]'
    except:
        alerts_sn = pd.DataFrame()

    for ra, dec in zip(ra_degraded[1:], dec_degraded[1:]):
        try:
            temp_sn = access_alerts(class_names = ['SN'],
                                    order_by = 'firstmjd',
                                    order_mode = 'ASC',
                                    firstmjd = [gw_mjd, gw_mjd + days_forward],
                                    ra = ra_degraded[0],
                                    dec = dec_degraded[0],
                                    radius = rad
                                    )
            #mask=(temp_sn['firstmjd']<gw_mjd+30)
            #temp_sn=temp_sn[mask]
            alerts_sn = alerts_sn.append(temp_sn, ignore_index=True)
        except:
            pass
        
        try:
            temp_agn = access_alerts(class_names = ['SN'],
                                    order_by = 'firstmjd',
                                    order_mode = 'ASC',
                                    firstmjd = [gw_mjd, gw_mjd + days_forward],
                                    ra = ra_degraded[0],
                                    dec = dec_degraded[0],
                                    radius = rad
                                    )
            #mask=(temp_sn['firstmjd']<gw_mjd+30)
            #temp_sn=temp_sn[mask]
            alerts_agn = alerts_agn.append(temp_agn, ignore_index=True)
        except:
            pass
        

    if alerts_sn.empty and alerts_agn.empty:
        print("No transient alert matches (SN or AGN) found within 4 degrees of CI contour.")
        print("... this is unlikely. Something probably went wrong with call to broker!")
        print("Exitting.")
        sys.exit()
        
    alerts_sn = alerts_sn.drop_duplicates(subset='oid')
    alerts_agn = agn_old_alerts.drop_duplicates(subset='oid')

    alerts_sn_ra = alerts_sn['meanra'].to_numpy()
    alerts_sn_dec = alerts_sn['meandec'].to_numpy()

    alerts_agn_ra = alerts_agn['meanra'].to_numpy()
    alerts_agn_dec = alerts_agn['meandec'].to_numpy()

    # print("Number of SN alerts:", alerts_sn_ra.size)
    # print("Number of AGN alerts:", alerts_agn_ra.size)

# *************************************************************************
#
# FIND TRANSIENT ALERTS IN PROBABILITY REGION
#
# *************************************************************************

    # Convert ALERCE data to pixel locations on GW map for ease of comparison (since CI pixel range is continuous, no fancy matching necessary)
    transient_pix_locs_sn = hp.ang2pix(gw_properties["nside"], alerts_sn_ra, alerts_sn_dec, lonlat = True, nest = gw_properties["nest"])
    transient_pix_locs_agn = hp.ang2pix(gw_properties["nside"], alerts_agn_ra, alerts_agn_dec, lonlat = True, nest = gw_properties["nest"])

    # Checks for matches in region 
    trans_in_prob_region_sn = np.isin(transient_pix_locs_sn, pixmap)
    trans_in_prob_region_agn = np.isin(transient_pix_locs_agn, pixmap)

    if np.any(trans_in_prob_region_sn) or np.any(trans_in_prob_region_agn):
        print(trans_in_prob_region_sn.sum(), "SN matches")
        print(trans_in_prob_region_agn.sum(), "AGN matches")

    alerce_sn_rows = alerts_sn.loc[trans_in_prob_region_sn, :].sort_values(by='firstmjd')
    alerce_agn_rows = alerts_agn.loc[trans_in_prob_region_agn, :].sort_values(by='firstmjd')
    
# *************************************************************************
#
# PLOT TRANSIENT ALERTS IN CI CONTOUR
#
# *************************************************************************

    plt.scatter(ra_map, dec_map, alpha=0.6)
    plt.scatter(alerce_agn_rows['meanra'], alerce_agn_rows['meandec'], label='AGN', s=4, alpha=0.7)
    plt.scatter(alerce_sn_rows['meanra'], alerce_sn_rows['meandec'], label='SN', alpha=0.8)
    plt.legend()
    
    plt.savefig(user_home + "/GW_DESI_Alerce.png")
    plt.clf()


# *************************************************************************
#
# WRITE TRANSIENT ALERTS ToO LEDGER TO ECSV
#
# *************************************************************************
    # ## Write ToO ledger to ecsv!
    
    '''Checker is Matthew Portman/Antonella Palmese 
    So if anything seems off... blame Antonella
    .
    .
    .
    if you dare.
    '''

    ## Only add SNe for now. AGNs should pass redshift cut. Also, Try rejecting stars
    #write_too_ledger(user_home + '/testing_ToO-Alerce.ecsv', alerce_sn_rows, checker='MP/AP', overwrite=True, verbose=True, tabformat='ALERCE')
    
# *************************************************************************
# *************************************************************************
# *************************************************************************
#
#
# PART 2 - 
# MATCHING DESI OBSERVATIONS TO X% CI CONTOUR
# BUILDING DR9 TARGETLIST
#
#
# *************************************************************************
# *************************************************************************
# *************************************************************************

# ## Matching DESI observations to x% CI contour
# This merely matches skymap angles to tile pointings and indicates which tiles matched and their program (bright/dark).
# 
# *The code does not currently write anything to file*

# ## The two *major* matching functions, **initial_check** and **inner_matching**
# ## Initial check
# Performs the initial match with tile pointings. Grabs info from the exposures sql table, filters it, condenses by date, then sends all the good stuff to be matched via **inner_matching**. Finally outputs the results from that into a dictionary with the dates as keys and the sql table information as the values (sqlite3 row type, similar to a namedtuple/well-keyed dictionary). So there can be multiple elements as the values but only one key hence it's all organized by date. A convenient thing but also necessary for ALERCE efficiency.
# 
# e.g. {date:(exposure info)}

# ## Matching CI interval pixel locations to DESI tile pointings

# TODO: Does this need to be done? *********************************************************************

#    desi_matches, _ = initial_check(ra_map, dec_map)

    # ## Print matched tile info (and Bright/Dark)
    # To be used later for plotting

    #print('tilera, tiledec, obsdate, obsmjd, expid, tileid')

    #tile_ras = list(set(tile_ras))
    #tile_decs = list(set(tile_decs))
    #print(len(bd_tile_ras))


    # Plotting the CI pixel map against the matched tiles from DESI
    # These are only matched within 2 deg so edge matches are to be expected
    #plt.scatter(ra_map, dec_map)
    #plt.scatter(bd_tile_ras, bd_tile_decs)

# *************************************************************************
#
# DR9 TARGETLIST BUILDING/REBUILDING/READING-IN
#
# *************************************************************************
# This uses the same matching philosophy as in *light_transient_matching* in that it finds matches to tile pointings and then performs a second stage check to find individual targets. Two methods are offered, matching by RA/DEC and matching by TargetID. As mentioned earlier, RA/DEC seems more robust.
# 
# The target lists are built in the same way as [Segev's code](https://github.com/desihub/timedomain/blob/master/gwtarget/gw_dr9.ipynb)

    trg_file = targetlists_path + '/' + gw_name + str(int((100*CI_val[0]))) + '_dr9.ecsv'
    
    try:
        targlist = Table.read(trg_file)
    except FileNotFoundError:
        targlist = build_targlist_table(gw_properties["nside"], pixmap)
        targlist_write(targlist, trg_file, overwrite)
        
        
# *************************************************************************
#
# INITIAL (VIA TILE) AND CLOSE MATCH (VIA FIBERS) TO TARGETLIST
#
# *************************************************************************

    m_dict, _ = initial_check(np.array(targlist['RA']), np.array(targlist['DEC']))
    
    if not m_dict:
        print("No matches found in initial matching (4 degree).")
        print("Quitting.")
        sys.exit()
    
    # As a reminder, uses original targlist data to find 1 arcsecond matches to individual targets via fibers
    desi_target_matches, targlist_target_matches = closer_check(matches_dict = m_dict, catalog2_ras = np.array(targlist['RA']),     catalog2_decs =  np.array(targlist['DEC']))
    if not targlist_target_matches:
        print("No matches found in closer matching (1 arcsecond).")
        print("Quitting.")
        sys.exit()

    # Some data reduction to avoid repeats (there's quite a lot!)
    unique_targlist_target_matches = SkyCoord(list(set([(val.ra.deg, val.dec.deg) for val in targlist_target_matches])), unit = 'deg')
    
    #print(len(targlist))
    #print(len(targlist_target_matches))
    #print(len(unique_targlist_target_matches))
    
    # For when 'around = True'... but as commented earlier, that's a *bad* idea in this case
    
    #m_dict_around, _ = \
    #    initial_check(np.array(targlist['RA']), np.array(targlist['DEC']), around = True)
    
    #closer_check(matches_dict = m_dict, catalog2_ras = targlist_matches.ra.deg, \
    #catalog2_decs = targlist_matches.ra.deg)


    # Generating lists of ras, decs, and program info for later graphing
    # ************** Need to ask about the below ***************
    """
    bd_tile_ras = []
    bd_tile_decs = []
    bd_program_info = []
    bd_tileids = []
    
    targetlist_tile_ras = []
    targetlist_tile_decs = []
    targetlist_program_info = []
    targetlist_tileids = []

    program_list = ['BRIGHT', 'DARK', 'BGS', 'ELG', 'QSO', 'LRG']

    # Iterate through intial_check dictionary outputs
    # I do both at the same time to be a little more efficient
    # .values()[i] is the sqlite3 row info for each 2 degree matched tile
    # the keys are the dates which we don't need here
    for i in range(np.max((len(desi_matches), len(m_dict)))):
        try:
            for j in list(desi_matches.values())[i]:
                p_name = j['program']
                # We check against tileids now instead of taking them out later so that
                # everything is retained in order
                if any(substring in p_name.upper() for substring in program_list) and j['tileid'] not in bd_tileids:
                    bd_tile_ras.append(j['tilera'])
                    bd_tile_decs.append(j['tiledec']) 
                    bd_program_info.append(p_name)
                    bd_tileids.append(j['tileid'])
        except:
            pass
        try:
            for j in list(m_dict.values())[i]:
                p_name = j['program']
                if any(substring in p_name.upper() for substring in program_list) and j['tileid'] not in targetlist_tileids:
                    targetlist_tile_ras.append(j['tilera'])
                    targetlist_tile_decs.append(j['tiledec']) 
                    targetlist_program_info.append(p_name)
                    targetlist_tileids.append(j['tileid'])
        except:
            pass
    """
    # ************** Need to ask about the above ***************

# *************************************************************************
#
# REDUCING TARGETLIST BY PREVIOUS OBSERVATIONS
#
# *************************************************************************

    # Convert targetlist matches to a Table to take advantage of astropy's setdiff method
    tlist_matches_table = Table([SkyCoord(unique_targlist_target_matches).ra.deg, 
                                 SkyCoord(unique_targlist_target_matches).dec.deg], names = ('RA', 'DEC'))

    #print(Table(targlist_matches))
    targlist_radec_reduced = setdiff(targlist['RA', 'DEC'], tlist_matches_table)

    assert len(targlist_radec_reduced) == len(targlist) - len(unique_targlist_target_matches), "Something went wrong masking the dr9 target list! Stopping."
    print('{:.2f}% of targets have already been observed within 2 degrees of DESI tile pointing in {} CI.'.format(100*len(unique_targlist_target_matches)/len(targlist), CI_val))
    
# *************************************************************************
#
# WRITING TO ToO ECSV
#
# *************************************************************************

    # This may take awhile, it's likely a lot of targets that we have yet to observe.

    #targlist_write(targlist_radec_reduced, "dr9_targlist" + CI_val + "_reduced_radec.ecsv", overwrite = True)

    # Write dr9 targets to ecsv
    write_too_ledger(filename = user_home + '/testing_ToO-bgs.ecsv', too_table = targlist_radec_reduced.to_pandas(), checker='MP/AP', overwrite=True, verbose=False, tabformat='LEGACY')

# *************************************************************************
#
# REDUCING TARGETLIST USING TARGETID
#
# *************************************************************************
    """
    # Seemingly less robust than RA/DEC but still worth holding onto. Somewhat faster but not as fast as one might hope (there is probably a quicker method to doing this however so don't take my word on speed until that alternative has been fully investigated (if it is ever investigated)) 

     targlist_targetid_reduced = grab_desi_targetid(matches_dict = m_dict, targlist_table = targlist, exclusion_list = [])

     #print(len(targlist))

     print('{:.2f}% of targets have already been observed within 2 degrees of DESI tile pointing in {} CI.'.format(100*(len(targlist) - len(targlist_targetid_reduced))/len(targlist), CI_val))

     # Print reduced targlist to file
     #targlist_write(targlist_targetid_reduced, "dr9_targlist" + str(100*CI_val[0]) +"_reduced_targetid.ecsv", overwrite = True)
     write_too_ledger(filename = user_home + '/testing_by_targetid_ToO-bgs.ecsv', too_table = targlist_targetid_reduced.to_pandas(), checker='MP/AP', overwrite=True, verbose=False, tabformat='LEGACY')
    """
# *************************************************************************
#
# PLOTTING RESULTS
#
# *************************************************************************

    # fig = plot_cartmap_tiles(gwfile, tile_ra = [218.9,217.5], tile_dec = [36.6,35.8], targ_ra = bd_tile_ras, targ_dec = bd_tile_decs, angsize = 5, program_names = program_info)
    # ax = fig.gca()
    # ax.set(xlim=(222, 213), ylim=(33,40))
    # #ax.scatter(exp_ras, exp_decs)
    # #ax.scatter(targlist_radec_reduced['RA'], targlist_radec_reduced['DEC'], alpha = 0.2, s = 0.3)
    # plt.savefig(user_home + '/' + gw_name + '_desi_tile-matches.png', dpi=120)
    
    # Plot targetlist tile matches and DESI FOV, just in case 
    fig = plot_cartmap_tiles(gwfile, tile_ra = [218.9,217.5], tile_dec = [36.6,35.8], \
                             targ_ra = targetlist_tile_ras, targ_dec = targetlist_tile_decs, angsize = 5, program_names = [])
    ax = fig.gca()
    ax.set(xlim=(222, 213), ylim=(33,40))
    ax.scatter(targlist_radec_reduced['RA'], targlist_radec_reduced['DEC'], alpha = 0.2, s = 0.3)
    #ax.scatter(targlist['RA'], targlist['DEC'])
    plt.savefig(user_home + '/' + gw_name + '_desi_tile-matches.png', dpi=120)
