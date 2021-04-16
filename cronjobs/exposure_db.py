import os
import numpy as np

from astropy.io import fits
from astropy.time import Time

import ephem

from glob import glob

import fitsio

import re
import sqlite3

from argparse import ArgumentParser


class ReduxDB:

    def __init__(self, db_file):
        """Initialize DB with tile and exposure tables.
    
        Parameters
        ----------
        db_file: str
            Database file name.
        """
        self.conn = None
        try:
            self.conn = sqlite3.connect(db_file)
        except sqlite3.Error as e:
            print(e)

        # Create the tiles table.
        sql = """ CREATE TABLE IF NOT EXISTS tiles (
                        tileid integer NOT NULL PRIMARY KEY,
                        tilera float NOT NULL,
                        tiledec float NOT NULL,
                        faflavor text NOT NULL
                    ); """
        self.create_table(sql)

        # Create the exposures table.
        sql = """ CREATE TABLE IF NOT EXISTS exposures (
                        expid integer NOT NULL PRIMARY KEY,
                        tileid integer NOT NULL,
                        program text NOT NULL,
                        obsdate integer NOT NULL,
                        obstime text NOT NULL,
                        obsmjd float NOT NULL,
                        exptime float NOT NULL,
                        airmass float NOT NULL,
                        moonalt float NOT NULL,
                        moonsep float NOT NULL,
                        moonfrac float NOT NULL,
                        FOREIGN KEY (tileid) REFERENCES tiles (tileid)
                    ); """
        self.create_table(sql)

    def add_tile(self, tileid, tilera, tiledec, faflavor):
        """Add a tile to the tiles table.

        Parameters
        ----------
        tileid : int
            Integer ID of the tile. Should always be unique.
        tilera : float
            RA of the tile, as defined in fiberassign file.
        tiledec : float
            Dec of the tile, as defined in fiberassign file.
        faflavor : str
            Fiberassign flavor: cmx, sv1lrgqso, etc.
        """
        sql = """ INSERT OR IGNORE INTO tiles(tileid,tilera,tiledec,faflavor)
                      VALUES(?,?,?,?) """
        self.insert_row(sql, (tileid, tilera, tiledec, faflavor))

    def add_exposure(self, expid, tileid, program, obsdate, obstime, obsmjd, exptime, airmass, moonalt, moonsep, moonfrac):
        """Add an exposure to the exposures table.

        Parameters
        ----------
        expid : int
            Integer ID of the exposure. Should always be unique.
        tileid : int
            Integer ID of the tile. Should always be unique.
        program : str
            Program defined when exposure is executed.
        obsdate : int
            Observation date at start of night, in YYYYMMDD.
        obstime : str
            Observation time: YYYY-MM-DDTHH:MM:SS.SSSS [UT]
        obsmjd : float
            MJD of observation.
        exptime : float
            Exposure time in seconds.
        airmass : float
            Airmass at the start of the exposure.
        moonalt : float
            Moon altitude (degrees) at the start of the exposure.
        moonsep : float
            Moon separation (degrees) at the start of the exposure.
        moonfrac : float
            Moon illumination fraction (phase) in [0,1].
        """
        sql = """ INSERT OR IGNORE INTO exposures(expid,tileid,program,obsdate,obstime,obsmjd,exptime,airmass,moonalt,moonsep,moonfrac)
                      VALUES(?,?,?,?,?,?,?,?,?,?,?) """
        self.insert_row(sql, (expid, tileid, program, obsdate, obstime, obsmjd, exptime, airmass, moonalt, moonsep, moonfrac))

        
    
    def create_table(self, create_table_sql):
        """ Create a table from the create_table_sql statement.
        
        Parameters
        ----------
        create_table_sql: str
            A CREATE TABLE statement.
        """
        try:
            c = self.conn.cursor()
            c.execute(create_table_sql)
        except sqlite3.Error as e:
            print(e)
        

    def insert_row(self, sql, rowinfo):
        """Insert a row into an SQL table.
        
        Parameters
        ----------
        sql : str
            SQL insertion string.
        rowinfo : tuple
            Tuple of data to be inserted into the row.
        
        Returns
        -------
        rowid : int
            Last inserted row ID.
        """
        cur = self.conn.cursor()
        cur.execute(sql, rowinfo)
        self.conn.commit()
        return cur.lastrowid


if __name__ == '__main__':
    p = ArgumentParser(description='DB of spectroscopic reductions.')
    p.add_argument('redux', nargs=1,
                   help='reduction version: daily, andes, blanc, ...')
    args = p.parse_args()

    redux = args.redux[0]

    # Create the spectroscopic reductions DB. It will contain two tables:
    # 1) tiles, with relatively static tile data.
    # 2) exposures, with more dynamic information.
    #Move '{}.db'.format(redux) to tablename
    db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'

    db = ReduxDB(db_filename)
    
    #Save the expids that are already in the db
    conn = sqlite3.connect(db_filename)
    c = conn.cursor()    
    check_expids_query="select expid from exposures;"
    cursor = c.execute(check_expids_query)
    expids_list = cursor.fetchall()
    expids_done  = np.array(expids_list).flatten()
    conn.close()

    # Observatory object for some astronomical calculations.
    mayall = None

    # Loop through date folders and identify subfolders with exposures.
    
    dateexpfolder = os.environ['DESI_SPECTRO_REDUX'] + '/{}/exposures'.format(redux)
    
    # Loop through observation date subfolders.
    dates = sorted(glob('{}/*'.format(dateexpfolder)))
    
    for date in dates:
        # Extract obsdate from the folder name.
        obsdate = date.split('/')[-1]

        # Check that obsdate is of format 20YYMMDD
        if not re.match(r'20\d{6}', obsdate):
            continue
        obsdate = int(obsdate)
        
        #Only look from March 1st 2021
        if obsdate<20210301:
            continue
        print('  + {}'.format(date))
    
        expfolder = os.environ['DESI_SPECTRO_REDUX'] + '/{}/exposures/{:08d}/'.format(redux,obsdate)
        exps = sorted(glob('{}/*'.format(expfolder)))

        expnums=[]
        for exp in exps:
            # Extract tileid from the folder name.
            expnum = exp.split('/')[-1]
            try:
                expnum = int(expnum)
                expnums.append(expnum)
            except:
                print("Expnum probably not a number")
            print(expnum)
            
            
        expids_all = np.array(expnums)

        # Now check if this expid is already in db
        # If it is don't bother with opening files
        expids_arr = np.setdiff1d(expids_all,expids_done)
        expids = set(expids_arr)

        if (expids_arr.shape[0]>0):
            print('Found new unique expids',expids_arr)

            # Access fiberassign file and DESI exposure file for each exposure.
            for expid in expids:
                cframe = os.environ['DESI_SPECTRO_REDUX'] + '/{}/exposures/{:08d}/{:08d}/cframe-b0-{:08d}.fits'.format(redux,obsdate,expid, expid)
                #desiexp = os.environ['DESI_SPECTRO_DATA'] + '/{}/{:08d}/desi-{:08d}.fits.fz'.format(expid, expid)
                print('trying to open ',cframe)

                if os.path.exists(cframe):
                    # Add tile info.
                    h=fits.open(cframe)
                    fhdr = h[1].header
                    fa_details = h[5].header

                    _id, _ra, _dec = [fa_details[_] for _ in ['TILEID', 'TILERA', 'TILEDEC']]
                    if 'FAFLAVOR' in fa_details:
                        _flav = fa_details['FAFLAVOR']
                    elif 'FA_SURV' in fa_details:
                        _flav = fa_details['FA_SURV']
                    else:
                        _flav = 'None'

                    db.add_tile(_id, _ra, _dec, _flav)

                    mra, mdec = fa_details['MOONRA'], fa_details['MOONDEC']
                    if 'MOONSEP' in fa_details:
                        msep = fa_details['MOONSEP']
                    else:
                        cosA = np.sin(_dec)*np.sin(mdec) + np.cos(_dec)*np.cos(mdec)*np.cos(_ra - mra) 
                        msep = np.degrees(np.arccos(cosA))

                    mjd, date, time = fa_details['MJD-OBS'], fa_details['NIGHT'], fa_details['DATE-OBS']
                    lat, lon, elev = fa_details['OBS-LAT'], fa_details['OBS-LONG'], fa_details['OBS-ELEV']

                    if mayall is None:
                        mayall = ephem.Observer()
                        mayall.lat = lat
                        mayall.lon = lon
                        mayall.elev = elev
                    mayall.date = Time(mjd, format='mjd', scale='utc').iso
                    moon = ephem.Moon(mayall)
                    mfrac = moon.moon_phase
                    malt = moon.alt * 180/np.pi
                    
                    
                    prog = fa_details['PROGRAM']
                    et = fa_details['EXPTIME']
                    am = fa_details['EXPTIME']

                    print('Adding   - {}'.format((expid, _id, prog, date, time, mjd, et, am, malt, msep, mfrac)))
                    db.add_exposure(expid.item(), _id, prog, date, time, mjd, et, am, malt, msep, mfrac)

