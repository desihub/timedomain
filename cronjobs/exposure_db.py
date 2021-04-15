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

    # Loop through tile folders and identify subfolders with exposures.
    tilesfolder = os.environ['DESI_SPECTRO_REDUX'] + '/{}/tiles'.format(redux)
    tiles = sorted(glob('{}/*'.format(tilesfolder)))

    for tile in tiles:
        # Extract tileid from the folder name.
        tileid = tile.split('/')[-1]

        #Switching from 5 digit tiled id in SV3
        # Now we check it can be translated into int
        
        # Check that tile ID is a 5-digit integer.
        #if not re.match(r'\d{5}', tileid):
        #    continue
        try:
            tileid = int(tileid)
        except:
            print("Tileid probably not a number")
        print(tileid)

        # Loop through observation date subfolders.
        dates = sorted(glob('{}/*'.format(tile)))
        for date in dates:
            # Extract obsdate from the folder name.
            obsdate = date.split('/')[-1]

            # Check that obsdate is of format 20YYMMDD
            if not re.match(r'20\d{6}', obsdate):
                continue
            obsdate = int(obsdate)
            print('  + {}'.format(date))

            cframefiles = sorted(glob('{}/cframe*.fits'.format(date)))

            # Extract unique list of expids from the cframe filenames.
            expids_all = np.array([int(cframefile.split('-')[-1][:-5]) for cframefile in cframefiles])
            
            # Now check if this expid is already in db
            # If it is don't bother with opening files
            expids_arr = np.setdiff1d(expids_all,expids_done)
            print(expids_all)
            
            expids = set(expids_arr)
            
            if (expids_arr.shape[0]>0):
                print('Found new unique expids',expids_arr)
                spectfiles = sorted(glob('{}/spectra*.fits'.format(date)))
                coaddfiles = sorted(glob('{}/coadd*.fits'.format(date)))

                # Access fiberassign file and DESI exposure file for each exposure.
                for expid in expids:
                    fassign = os.environ['DESI_SPECTRO_DATA'] + '/{}/{:08d}/fiberassign-{:06d}.fits'.format(obsdate, expid, tileid)
                    desiexp = os.environ['DESI_SPECTRO_DATA'] + '/{}/{:08d}/desi-{:08d}.fits.fz'.format(obsdate, expid, expid)

                    if os.path.exists(fassign) and os.path.exists(desiexp):
                        # Add tile info.
                        fhdr = fitsio.read_header(fassign, 'PRIMARY')

                        _id, _ra, _dec = [fhdr[_] for _ in ['TILEID', 'TILERA', 'TILEDEC']]
                        if 'FAFLAVOR' in fhdr:
                            _flav = fhdr['FAFLAVOR']
                        elif 'FA_SURV' in fhdr:
                            _flav = fhdr['FA_SURV']
                        else:
                            _flav = 'None'

                        db.add_tile(_id, _ra, _dec, _flav)

                        # Add exposure info.
                        dhdr = fitsio.read_header(desiexp, 1)

                        mra, mdec = dhdr['MOONRA'], dhdr['MOONDEC']
                        if 'MOONSEP' in dhdr:
                            msep = dhdr['MOONSEP']
                        else:
                            cosA = np.sin(_dec)*np.sin(mdec) + np.cos(_dec)*np.cos(mdec)*np.cos(_ra - mra) 
                            msep = np.degrees(np.arccos(cosA))

                        mjd, date, time = dhdr['MJD-OBS'], obsdate, dhdr['DATE-OBS']
                        lat, lon, elev = dhdr['OBS-LAT'], dhdr['OBS-LONG'], dhdr['OBS-ELEV']

                        if mayall is None:
                            mayall = ephem.Observer()
                            mayall.lat = lat
                            mayall.lon = lon
                            mayall.elev = elev
                        mayall.date = Time(mjd, format='mjd', scale='utc').iso
                        moon = ephem.Moon(mayall)
                        mfrac = moon.moon_phase
                        malt = moon.alt * 180/np.pi

                        prog = dhdr['PROGRAM']
                        et = dhdr['EXPTIME']
                        am = dhdr['EXPTIME']

                        print('    - {}'.format((expid, tileid, prog, date, time, mjd, et, am, malt, msep, mfrac)))
                        db.add_exposure(expid.item(), tileid, prog, date, time, mjd, et, am, malt, msep, mfrac)

