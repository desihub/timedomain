import sqlite3

from astropy.table import Table
from astropy.time import Time

import numpy as np

from argparse import ArgumentParser

class TooDB:

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

    def add_alert(objid, instrument, date, mjd, ra, dec):
        try:
            # Check if this specific alert is already loaded.
            search_query = f"""
                SELECT too_id, obj_id, discovery_mjd, event_number FROM toolist
                WHERE obj_id='{objid}';
            """
            self.curs.execute(search_query)
            results = self.curs.fetchall()
            if results:
                print(f'Alert {obj_id} already in ToO DB.')
                return 0
            
            # Check for other alerts from this MJD.
            search_query = f"""
                SELECT too_id, obj_id, discovery_mjd, event_number FROM toolist
                WHERE discovery_mjd='{mjd};'
            """
            self.curs.execute(search_query)
            results = self.curs.fetchall()
            number = len(results) + 1

            # Encode a ToO ID and write the entry to the database.
            tooid = encode_tooid(mjd, number)

            insert_query = """
                INSERT INTO toolist
                (too_id, obj_id, instrument, discovery_date, discovery_mjd, event_number, ra, dec)
                VALUES
                (?, ?, ?, ?, ?, ?);
            """
            data = (tooid, objid, instrument, date, mjd, number, ra, dec))
            toodb.curs.execute(insert_query, data)
            toodb.conn.commit()

            return tooid

        except sqlite.Error as err:
            print(err)
            return 0

#try:
#    c = sqlite3.connect('too_list.db')
#    c.execute("""
#              CREATE TABLE IF NOT EXISTS toolist
#              ([too_id] INTEGER PRIMARY KEY,
#               [obj_id] TEXT,
#               [instrument] TEXT,
#               [discovery_date] TEXT,
#               [discovery_mjd] INTEGER,
#               [event_number] INTEGER,
#               [ra] REAL,
#               [dec] REAL)
#              """)
#    c.commit()
#except sqlite3.Error as err:
#    print('Failed to create table "toolist"', err)
#finally:
#    if c:
#        c.close()

def encodetooid(mjd, number):
    """Encode MJD and alert number in MJD into a 32b ID.

    Parameters
    ----------
    mjd : int
        Modified Julian Date of the alert (valid thru 10 Nov 2054).
    number : int
        Alert number within given MJD (up to 262144).

    Returns
    -------
    tooid : int
        32b ID number for TOO.
    """
    tooid = ((mjd - 55197) << 18) + number
    return tooid

def decodetooid(tooid):
    """Encode MJD and alert number in MJD into a 32b ID.

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

if __name__ == '__main__':
    p = ArgumentParser(description='ToO Ledger maker')
    p.add_argument('inputfile', nargs=1, help='Input TOO file.')
    args = p.parse_args()

    alert_count = {}

    data = Table.read(args.inputfile[0])
    toodb = TooDB()

#    if 'XWIN_WORLD' in data.columns:

    try:
#        c = sqlite3.connect('too_list.db')

        for entry in data:
            objid, ra, dec = entry['ObjectID'], entry['XWIN_WORLD'], entry['YWIN_WORLD']
            yr, mo, dy = objid[1:5], objid[5:7], objid[7:9]
            date = f'{yr}-{mo}-{dy}T00:00:00.1'
            t = Time(f'{date}', format='isot', scale='utc')
            mjd = int(np.floor(t.mjd))

            if mjd in alert_count:
                alert_count[mjd] += 1
            else:
                alert_count[mjd] = 1
            number = alert_count[mjd]

            tooid = encodetooid(mjd, number)
            print(objid, ra, dec, t, mjd, number, tooid)
            print(decodetooid(tooid))

            try:
                insert_query = f"""
                    INSERT INTO toolist
                    (too_id, obj_id, instrument, discovery_date, discovery_mjd, event_number, ra, dec)
                    VALUES
                    ({tooid}, '{objid}', 'DECam', '{date}', {mjd}, {number}, {ra}, {dec})
                """
                toodb.curs.execute(insert_query)
                toodb.conn.commit()
            except sqlite3.Error as err:
                print('Failed to insert data into sqlite table "toolist"', err)

            print(insert_query)
    except sqlite3.Error as err:
        print(err)
        raise SystemExit
#    finally:
#        if c:
#            c.close()
