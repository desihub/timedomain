"""Access classes for a database of external ToO alerts."""

import sqlite3
import numpy as np


class TooDB:
    """Class-level access to database with processed Target-of-Opportunity
    (ToO) alerts.
    """

    def __init__(self, filename='too_list.db'):
        try:
            self.filename = filename
            self.conn = sqlite3.connect(self.filename)
            self.curs = self.conn.cursor()
            self.conn.execute("""CREATE TABLE IF NOT EXISTS toolist
              ([too_id] INTEGER PRIMARY KEY,
               [obj_id] TEXT,
               [instrument] TEXT,
               [discovery_date] TEXT,
               [discovery_mjd] INTEGER,
               [event_number] INTEGER,
               [ra] REAL,
               [dec] REAL)
              """)
            self.conn.commit()
        except sqlite3.Error as err:
            print(f'Failed to open DB {filename}', err)
            raise SystemExit

    def add_alert(self, objid, instrument, date, mjd, ra, dec):
        """Add a ToO alert to the ToO sqlite database.

        Parameters
        ----------
        objid : str
            Alert name or object ID.
        instrument : str
            Name of observing instrument.
        date : str
            Date string in ISO format YYYY-MM-DDTHH:MM.SS.S.
        mjd : int
            Modified Julian Date.
        ra : float
            Right ascension, in degrees.
        dec : float
            Declination, in degrees.

        Returns
        -------
        tooid : int
            ToO ID if successful addition; 0 if not.
        """
        try:
            # Check if this specific alert is already loaded.
            search_query = f"""
                SELECT too_id, obj_id, discovery_mjd, event_number FROM toolist
                WHERE obj_id='{objid}';
            """
            self.curs.execute(search_query)
            results = self.curs.fetchall()
            if results:
                print(f'Alert {objid} already in ToO DB.')
                return 0
            
            # Check for other alerts from this MJD.
            search_query = f"""
                SELECT too_id, obj_id, discovery_mjd, event_number FROM toolist
                WHERE discovery_mjd={mjd};
            """
            self.curs.execute(search_query)
            results = self.curs.fetchall()
            number = len(results) + 1

            # Encode a ToO ID and write the entry to the database.
            tooid = self.encode_tooid(mjd, number)

            insert_query = """
                INSERT INTO toolist
                (too_id, obj_id, instrument, discovery_date, discovery_mjd, event_number, ra, dec)
                VALUES
                (?, ?, ?, ?, ?, ?, ?, ?);
            """
            data = (tooid, objid, instrument, date, mjd, number, ra, dec)
            self.curs.execute(insert_query, data)
            self.conn.commit()

            return tooid

        except sqlite3.Error as err:
            print(err)
            return 0

    def get_data(self):
        """Access all data in the ToO database.

        Returns
        -------
        data : np.recarray
            Array of data from the database.
        """
        try:
            search_query = 'SELECT * FROM toolist;'
            self.curs.execute(search_query)
            results = self.curs.fetchall()
            if results:
                return np.rec.array(results, names=[
                    'TOO_ID',
                    'OBJ_ID',
                    'INSTRUMENT',
                    'DISCOVERY_DATE',
                    'DISCOVERY_MJD',
                    'EVENT_NUMBER',
                    'RA',
                    'DEC'])
            else:
                return np.asarray(results)
        except sqlite3.Error as err:
            print(f'Could not access data in {self.filename}', err)
            return np.asarray([])

    def encode_tooid(self, mjd, number):
        """Helper to encode MJD and alert number in MJD into a 32b ID.

        Parameters
        ----------
        mjd : int
            Modified Julian Date of the alert (valid thru 9 Nov 2054).
        number : int
            Alert number within given MJD (up to 262144).

        Returns
        -------
        tooid : int
            32b ID number for TOO.
        """
        # MJD = 14 MSBs, number = 18 LSBs. Total TOOID = 32 bits.
        tooid = ((((mjd - 55197) & 0x3fff) << 18) + number) & 0xffffffff
        return tooid

    def decode_tooid(self, tooid):
        """Helper to decode MJD and alert number in MJD into a 32b ID.

        Parameters
        ----------
        tooid : int
            32b ID number for TOO.

        Returns
        -------
        mjd : int
            Modified Julian Date of the alert.
        number : int
            Alert number within given MJD.
        """
        mjd = 55197 + (tooid >> 18)
        number = tooid & ((1 << 18) - 1)
        return mjd, number

