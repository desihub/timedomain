#!/usr/bin/env python
"""Convert a native GW followup json obs queue to a DES wide compatible one
"""

import sys
import json
import logging

from argparse import ArgumentParser
from collections import defaultdict
from functools import partial
from math import sin, cos, radians, degrees, acos
from copy import deepcopy
from dateutil.parser import parse as parse_date
import datetime
import numpy as np
from astropy.time import Time
from astropy.coordinates import Longitude
from astropy.coordinates import Latitude
from astropy.coordinates import Angle
from astropy.units import deg
from astropy.units import radian

from spherical_geometry.polygon import SphericalPolygon
from spherical_geometry.vector import radec_to_vector
from .telescope import blancoHorizonLimits as load_blanco_limits

# constants

max_match_angle = 3.0/60.0
camera_radius = 1.3 # This should depend on the camera no?
ctio_latitude = Latitude(-30.16527778, unit=deg)
ctio_longitude = Longitude(-70.815, unit=deg)
ctio_elevation = 2215.0
kpno_latitude    = Latitude(31.9600784, unit=deg)
kpno_longitude   = Longitude(-111.598169, unit=deg)
kpno_elevation   =  2067.0

# exception classes
# interface functions

def gwwide(wide_queue, gw_queue, start_time, dt, sort, camera="decam"):
    """Replace entries in a GW quere with nearby entries in the DES wide survey

    :Parameters:
        - `wide_queue`: the DES wide survey queue (list of dicts)
        - `gw_queue`: the naieve GW queue (list of dicts)
        - `start_time`: the start time in mjd
        - `dt`: time buffer (seconds)
        - `sort`: sort by set time

    @returns: a queue (list of dicts) with entries replaced
    """
    if camera == "decam" :
        obs_longitude = ctio_longitude
        obs_latitude = ctio_latitude
        obs_elevation = ctio_elevation
    elif camera == "desi" :
        obs_longitude = kpno_longitude
        obs_latitude = kpno_latitude
        obs_elevation = kpno_elevation
    start_time = Time(start_time, format='mjd', scale='utc', location=(obs_longitude, obs_latitude))
    start_time.format="iso"
    start_time = start_time.value

    fixed_queue = [fix_obs(gw_obs, wide_queue, start_time, dt, camera) for gw_obs in gw_queue]

    if sort:
        fixed_queue = sorted(fixed_queue, key = lambda e: e['lateTime'])

    return fixed_queue

# classes
# internal functions & classes

def angle(obs1, obs2):
    """Angular distance between two queue entries (deg)

    :Parameters:
        - `obs1`: queue entry 1
        - `obs2`: queue entry 2

    @returns: angular distance in degress
    """

    ra1, decl1 = obs1['RA'], obs1['dec']
    ra2, decl2 = obs2['RA'], obs2['dec'] 
    if abs(ra1-ra2)<(1.0/(60*60*10000)) and abs(decl1-decl2)<(1.0/(60*60*10000)):
        return 0.0

    sep = degrees( acos( sin(radians(decl1))*sin(radians(decl2))
                         + cos(radians(decl1))*cos(radians(decl2))*cos(radians(ra1-ra2)) ) )
    return sep

def wide_covered(gw_obs, wide_queue):
    """Return whether a queue entry is covered by wide survey hexes

    Assume it is covered if the pointing is covered by hexes
    in all 10 DES wide tilings.

    :Parameters:
        - `gw_obs`: the queue entry to check
        - `wide_queue`: the list of all DES wide queue entrues

    @returns: true iff an obs is covered by wide survey hexes
    """
    # If we are within a camera radius of a pointing in every tiling, we are in the footprint
    for tiling_id in range(1,11):
        wide_in_tiling = [w for w in wide_queue 
                          if w['tiling_id']==tiling_id 
                          and w['filter']==gw_obs['filter']]
        nearest_in_tiling = min(wide_in_tiling, 
                                key=partial(angle, gw_obs))
        angle_to_nearest = angle(gw_obs, nearest_in_tiling)
        if angle_to_nearest > camera_radius:
            return False
    return True
    
def fix_obs(gw_obs, all_filter_wide_queue, start_time, dt, camera):
    """Adjust an exposure queue entry to match a wide field entry, if close

    :Parameters:
        - `gw_obs`: the GW obs queue entry
        - `all_filter_wide`: a list of all DES wide survey queue entrues
        - `start_time`: the start time
        - `dt`: time buffer (seconds)

    @returns: a queue entry, adjusted if necessary
    """
    if gw_obs['RA'] < 0:
        gw_obs['RA'] = gw_obs['RA']+ 360.0

    if gw_obs['RA'] > 360:
        gw_obs['RA'] = gw_obs['RA'] - 360.0

    gw_obs['earlyTime'] = rising_iso8601(gw_obs['RA'], gw_obs['dec'], start_time, dt, camera)
    gw_obs['lateTime'] = setting_iso8601(gw_obs['RA'], gw_obs['dec'], start_time, dt, camera)
    
    wide_queue = [w for w in all_filter_wide_queue 
                  if w['filter']==gw_obs['filter']]

    logging.debug("GW exposure:   %3.4f, %3.4f", gw_obs['RA'], gw_obs['dec'])

    no_match = False

    try:
        nearest_wide_obs = min(wide_queue, key=partial(angle, gw_obs))

        nearest_angle = angle(gw_obs, nearest_wide_obs)
        logging.debug("wide exposure: %3.4f, %3.4f (tiling %d) is %3.4f degrees away",
                      nearest_wide_obs['RA'], nearest_wide_obs['dec'], 
                      nearest_wide_obs['tiling_id'],
                      nearest_angle)

        if not 'exptime' in gw_obs:
            gw_obs['exptime'] = gw_obs['expTime']
    except:
        no_match = True
        nearest_angle = 180
        
    if nearest_angle > max_match_angle or gw_obs['exptime']!=90 or no_match:
        #logging.warning("Exposure at %(RA)3.4f %(dec)3.4f not matched", gw_obs)
        #for kw in ['count','seqtot','seqnum']:
        #    gw_obs[kw] = int(gw_obs[kw])

        return gw_obs

    fixed_obs = deepcopy(nearest_wide_obs)
    #for kw in ['count','seqid','seqnum','seqtot','note','comment','propid','earlyTime','lateTime']:
    for kw in ['note','comment','propid','earlyTime','lateTime']:
        fixed_obs[kw] = gw_obs[kw]

    #for kw in ['count','seqtot','seqnum']:
    #    fixed_obs[kw] = int(fixed_obs[kw])

    return fixed_obs

def file_gwwide(gw_fnames, wide_fname, start_time, dt, fixed_fname, sort):
    """Adjust exposure queue file to match a wide field entries, if close

    :Parameters:
        - `gw_fnames`: the list of GW obs queue json files
        - `wide_fname`: the json file with all DES wide field queue entries
        - `start_time`: the start date/time in ISO 8601 format
        - `dt`: time buffer (seconds)
        - `fixed_fname`: the json obs queue file into which to write entries
        - `sort`: sort exposures by setting time 

    """

    with open(wide_fname, 'r') as f:
        wide_queue = json.load(f)
    logging.info("Loaded %d wide survey entries from %s", len(wide_queue), wide_fname)

    if isinstance(gw_fnames, str):
        gw_fnames = [gw_fnames]

    gw_queue = []
    for gw_fname in gw_fnames:
        with open(gw_fname, 'r') as f:
            new_gw_queue = json.load(f)
        logging.info("Loaded %d GW followup entries from %s", len(new_gw_queue), gw_fname)
        gw_queue += new_gw_queue
        
    fixed_queue = gwwide(wide_queue, gw_queue, start_time, dt, sort)
    logging.info("Finished correcting queue entries")

    with open(fixed_fname, 'w') as f:
        json.dump(fixed_queue, f, indent=4)
    logging.info("Saved %d GW followup entries in %s", len(fixed_queue), fixed_fname)

    return 0

def lst(mjd, camera):
    """Return the local sidereal time at CTIO

    :Parameters:
        - `mjd`: the floating point MJD

    @Return: LST in degrees
    """
    if camera == "decam" :
        obs_longitude = ctio_longitude
        obs_latitude = ctio_latitude
        obs_elevation = ctio_elevation
    elif camera == "desi" :
        obs_longitude = kpno_longitude
        obs_latitude = kpno_latitude
        obs_elevation = kpno_elevation
    t = Time(mjd, scale='utc', format='mjd', location=(obs_longitude, obs_latitude))
    t.delta_ut1_utc = 0.
    local_sidereal_time = t.sidereal_time('mean', longitude=obs_longitude)
    lst_deg = float(local_sidereal_time / deg)
    return lst_deg

def in_blanco_limits(ra, decl, mjd):
    """Return whether a given pointing is within Blanco limits.

    :Parameters:
        - `ra`: the right ascension (decimal degrees)
        - `dec`: the declination (decimal degrees)
        - `mjd`: the modified julian date

    @Return: boolean true of false
    """
    ha = lst(mjd,"decam") - ra
    ha_limits, decl_limits = load_blanco_limits()
    ha_limits = np.degrees(ha_limits)
    decl_limits = np.degrees(decl_limits)
    
    ha_limits = np.append(ha_limits, ha_limits[0])
    decl_limits = np.append(decl_limits, decl_limits[0])
    blanco_limits = SphericalPolygon.from_radec(
        ha_limits, decl_limits,
        center=(180*deg, 0*deg), degrees=True)
    pointing = radec_to_vector(ha, decl, degrees=True)
    in_limits = blanco_limits.contains_point(pointing)
    return in_limits

def polygon_setting_mjd(ra, decl, mjd):
    dd = 1.0/(24*60)
    while not in_blanco_limits(ra, decl, mjd):
        mjd += 0.011

    while in_blanco_limits(ra, decl, mjd):
        before_set = mjd
        mjd += dd

    return before_set

def polygon_rising_mjd(ra, decl, mjd):
    dd = 1.0/(24*60)
    setting_mjd = polygon_setting_mjd(ra, decl, mjd)
    test_mjd = setting_mjd-0.012

    while in_blanco_limits(ra, decl, test_mjd):
        test_mjd -= 0.011

    while not in_blanco_limits(ra, decl, test_mjd):
        test_mjd += dd

    return test_mjd

def setting_ha(decl):
    """Return the HA at which a given declination sets

    :Parameters:
        - `decl`: the declination (in decimal degrees)

    @Returns: the HA (in decimal degrees) at which dec leaves Blanco limits
    """
    rad_ha_limits, rad_decl_limits = load_blanco_limits()
    ha_limits = np.degrees(rad_ha_limits)
    decl_limits = np.degrees(rad_decl_limits)
    setting_ixs = ( ha_limits < 0) & (decl_limits >= -90)  # Eric had -88; I do not know why.
    setting_has = ha_limits[setting_ixs]
    setting_decls = decl_limits[setting_ixs]
    setting_ixs = np.argsort(setting_decls)
    setting_decls = setting_decls[setting_ixs]
    setting_has = setting_has[setting_ixs]
    #assert np.all(np.diff(setting_decls) > 0)
    # some declinations never set
    ha_deg = np.interp(decl, setting_decls, setting_has)

    ha = float( Angle(ha_deg*deg).wrap_at(360*deg)/deg )
        
    return ha

def rising_ha(decl):
    """Return the HA at which a given declination sets

    :Parameters:
        - `decl`: the declination (in decimal degrees)

    @Returns: the HA (in decimal degrees) at which dec leaves Blanco limits
    """
    rad_ha_limits, rad_decl_limits = load_blanco_limits()
    ha_limits = np.degrees(rad_ha_limits)
    decl_limits = np.degrees(rad_decl_limits)
    rising_ixs = ( ha_limits < 0) & (decl_limits >= -90)  # Eric had -88; I do not know why.
    rising_has = ha_limits[rising_ixs]
    rising_decls = decl_limits[rising_ixs]
    if not rising_decls.min()<decl<rising_decls.max():
        raise ValueError

    rising_ixs = np.argsort(rising_decls)
    rising_decls = rising_decls[rising_ixs]
    rising_has = rising_has[rising_ixs]
    
    assert np.all(np.diff(rising_decls) >= 0)
    ha_deg = np.interp(decl, rising_decls, rising_has)

    ha = float( Angle(ha_deg*deg).wrap_at(360*deg)/deg )
        
    return ha

def setting_lst(ra, decl):
    """Return the LST (in degrees) at which a given set of coordiantes sets.

    :Parameters:
        - `ra`: the right ascension (in decimal degrees)
        - `decl`: the declination (in decimal degrees)

    @Returns: the LST (in decimal degrees) at which dec leaves Blanco limits
    """
    lst_deg = setting_ha(decl) + ra
    lst = float( Angle(lst_deg*deg).wrap_at(360*deg)/deg )
    return lst

def rising_lst(ra, decl):
    """Return the LST (in degrees) at which a given set of coordiantes sets.

    :Parameters:
        - `ra`: the right ascension (in decimal degrees)
        - `decl`: the declination (in decimal degrees)

    @Returns: the LST (in decimal degrees) at which dec leaves Blanco limits
    """
    lst_deg = rising_ha(decl) + ra
    lst = float( Angle(lst_deg*deg).wrap_at(360*deg)/deg )
    return lst

def setting_mjd(ra, decl, mjd, dt=0.0, camera="decam"):
    """Return the next setting time after MJD

    :Parameters:
        - `ra`: the right ascension (in decimal degrees)
        - `decl`: the declination (in decimal degrees)
        - `mjd`: the mjd after which to find the next setting
        - `dt`: seconds to leave before hitting the Blanco limits

    @Returns: the setting MJD
    """
    set_lst = Angle(setting_lst(ra, decl), unit=deg)
    current_lst = Angle(lst(mjd,camera), unit=deg)
    deg_in_future = float( (set_lst-current_lst).wrap_at(360*deg)/deg )

    # factor of 1/360 converts from degrees to sidereal days,
    # 365.2422/366.2422 from sidereal days to solar days
    set_mjd = mjd+(deg_in_future*(365.2422/366.2422)/360.0)


    # these don't handle near the south pole where things never set
    # JTA Feb 2020
    #assert in_blanco_limits(ra, decl, set_mjd-1.0/(24*60*60))
    #assert (not in_blanco_limits(ra, decl, set_mjd + 1.0/(24*60*60)))

    set_mjd = set_mjd - dt/(24*60*60.0)
    return set_mjd

def rising_mjd(ra, decl, mjd, dt=0.0, camera="decam"):
    """Return the next rising time after MJD

    :Parameters:
        - `ra`: the right ascension (in decimal degrees)
        - `decl`: the declination (in decimal degrees)
        - `mjd`: the mjd after which to find the next rising
        - `dt`: seconds to leave before hitting the Blanco limits

    @Returns: the rising MJD
    """
    rise_lst = Angle(rising_lst(ra, decl), unit=deg)
    current_lst = Angle(lst(mjd,camera), unit=deg)
    deg_in_past = float( (current_lst-rise_lst).wrap_at(360*deg)/deg )

    # factor of 1/360 converts from degrees to sidereal days,
    # 365.2422/366.2422 from sidereal days to solar days
    rise_mjd = mjd-(deg_in_past*(365.2422/366.2422)/360.0)

    assert in_blanco_limits(ra, decl, rise_mjd+1.0/(24*60*60))
    assert (not in_blanco_limits(ra, decl, rise_mjd - 1.0/(24*60*60)))

    rise_mjd = rise_mjd + dt/(24*60*60.0)
    return rise_mjd

def setting_iso8601(ra, decl, tstring, dt=0.0, camera="decam"):
    """Return the ISO 8601 string representation of the setting tem

    :Parameters:
        - `ra`: the right ascension (in decimal degrees)
        - `decl`: the declination (in decimal degrees)
        - `tstring`: the string representation of the time
        - `dt`: seconds to leave before hitting the Blanco limits

    @Returns: the setting date string
    """
    t = Time(parse_date(tstring))
    set_mjd = setting_mjd(ra, decl, t.mjd, dt, camera)
    set_t = Time(set_mjd, scale='tt', format='mjd', out_subfmt='date_hms')
    set_date = set_t.iso[:10]
    set_time = set_t.iso[11:19]
    set_iso = set_date + 'T' + set_time + 'Z'
    return set_iso
    

def rising_iso8601(ra, decl, tstring, dt=0.0, camera="decam"):
    """Return the ISO 8601 string representation of the rising time

    :Parameters:
        - `ra`: the right ascension (in decimal degrees)
        - `decl`: the declination (in decimal degrees)
        - `tstring`: the string representation of the time
        - `dt`: seconds to leave before hitting the Blanco limits

    @Returns: the rising date string
    """
    t = Time(parse_date(tstring))
    rise_mjd = rising_mjd(ra, decl, t.mjd, dt, camera)
    rise_t = Time(rise_mjd, scale='tt', format='mjd', out_subfmt='date_hms')
    rise_date = rise_t.iso[:10]
    rise_time = rise_t.iso[11:19]
    rise_iso = rise_date + 'T' + rise_time + 'Z'
    return rise_iso
    

def main():
    current_timestring = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    
    parser = ArgumentParser(
        description='Fix a GW followup script to match wide-survey exposures')    
    parser.add_argument('wide_script', help='Filename of wide obs script')
    parser.add_argument('fixed_script', help='Filename we write fixed script to')
    parser.add_argument('gw_scripts', nargs='+', help='Filename of GW obs script')
    parser.add_argument('-t', '--start_time', help='Start time (ISO 8601)', default=current_timestring)
    parser.add_argument('-d', '--dt', help='time buffer before setting (seconds)', type=float, default=300)
    parser.add_argument('-s', '--sort', action="store_true")
    parser.add_argument('-v', '--verbose', action="count", 
                        help="be verbose")

    args = parser.parse_args()

    log_levels = defaultdict(lambda:logging.WARNING)
    log_levels[1] = logging.INFO
    log_levels[2] = logging.DEBUG
    logging.basicConfig(format='%(asctime)s %(message)s',
                        level=log_levels[args.verbose])

    # Fix format of start time
    st = Time(parse_date(args.start_time), out_subfmt='date_hms').iso
    st_date = st[:10]
    st_time = st[11:19]
    st_iso = st_date + 'T' + st_time + 'Z'
    
    logging.info("Begin")

    status = file_gwwide(args.gw_scripts, args.wide_script, args.start_time, args.dt,
                         args.fixed_script, args.sort)
    return status


if __name__ == '__main__':
    status = main()
    sys.exit(status)    

