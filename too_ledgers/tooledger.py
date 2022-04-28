"""This module implements TooAlertList handlers to manage a large variety of
possible ToO file formats from different sources, such as TNS, DECam, etc. To
add a new format, create a subclass of TooAlertList and add a format check to
the AlertListFactory to use the new format.
"""
import os

from toodb import TooDB

from abc import ABC, abstractmethod

import astropy.units as u
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord

import numpy as np


class AlertListFactory:
    """Factory class used to choose between ToO file formats."""

    def get_alert_handler(self, too_input_file):
        """Open a ToO input file and decide on the format.

        Parameters
        ----------
        too_input_file : str
            Path to ToO input file, in a format readable by astropy.table.

        Returns
        -------
        alerthandler : TooAlertList
            Handler object to parse the ToO list.
        """
        data = Table.read(too_input_file)

        if 'XWIN_WORLD' in data.columns:
            # DECam ToO format.
            return DECamAlertList(data)
        elif 'Reporting Group/s' in data.columns and 'Discovery Data Source/s' in data.columns:
            # TNS ToO format.
            return TNSAlertList(data)
        else:
            raise RuntimeError('Unrecognized ToO alert list format.')

factory = AlertListFactory()


class TooLedgerMaker:
    """Generic alert handler to manage ToO lists in different formats."""

    def build_too_ledger(self, too_input_file, too_ledger_file, verbose=False):
        """Open a ToO input file and build an alert ledger.

        Parameters
        ----------
        too_input_file : str
            Path to ToO input file, in a format readable by astropy.table.
        too_ledger_file : str
            ToO ledger, in a format appropriate for fiberassign.
        """
        # Alert factory decides on the file format and acts appropriately.
        alerthandler = factory.get_alert_handler(too_input_file)

        # Appropriate alert handler writes the output.
        alerthandler.write_too_list(too_ledger_file, verbose)


class TooAlertList(ABC):
    """Abstract ToO list handler that implements ToO ledger output.
    """
    def __init__(self, too_table):
        # Set up the ToO data table and database connection.
        self.too_table = too_table
        self.toodb = TooDB()

    @abstractmethod
    def generate_too_list(self):
        """Find unique alerts not in the ToO database.
        """
        pass

    def write_too_list(self, outputfile, verbose=False):
        """Append ToO data to an ecsv file for processing by fiberassign.

        Parameters
        ----------
        outputfile : str
            Path to output ecsv file.
        verbose : bool
            Print additional info if true.
        """
        too_list = self.generate_too_list()

        with open(outputfile, 'a') as outf:
            header = """# %ECSV 0.9
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
RA DEC PMRA PMDEC REF_EPOCH CHECKER TOO_TYPE TOO_PRIO OCLAYER MJD_BEGIN MJD_END TOOID"""
            if os.path.exists(outputfile):
                if os.path.getsize(outputfile) == 0:
                    outf.write(f'{header}\n')
            else:
                outf.write(f'{header}\n')

            for entry in too_list:
                ra, dec, pmra, pmdec, epoch, checker, too_type, prio, prog, mjd0, mjd1, too_id = entry
                output = '{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {}'.format(ra, dec, pmra, pmdec, epoch, checker, too_type, prio, prog, mjd0, mjd1, too_id)
                if verbose:
                    print(output)
                outf.write(f'{output}\n')


class DECamAlertList(TooAlertList):
    """ToO input handling from DECam/DESIRT.
    """
    def __init__(self, too_table):
        super().__init__(too_table)

    def generate_too_list(self):
        """Find unique alerts not in the ToO database.

        Returns
        -------
        too_list : tuple or list
            Per-alert info needed for the ToO ledger.
        """

        too_list = []

        for entry in self.too_table:
            # Extract the relevant data for the ToO alert.
            objid, ra, dec = entry['ObjectID'], entry['XWIN_WORLD'], entry['YWIN_WORLD']
            yr, mo, dy = objid[1:5], objid[5:7], objid[7:9]
            date = f'{yr}-{mo}-{dy}T00:00:00.1'
            t = Time(f'{date}', format='isot', scale='utc')
            mjd = int(np.floor(t.mjd))

            # Enter data into the DB. 
            tooid = self.toodb.add_alert(objid, 'DECam', date, mjd, ra, dec)
            if tooid == 0:
                continue

            # Accumulate data for output:
            # RA, DEC, PMRA, PMDEC, EPOCH, CHECKER, TYPE, PRIO, PROG, MJD_START, MJD_STOP, TOO_ID
            too_list.append(
                [ra, dec, 0., 0., 2000.0, 'SB/AP', 'FIBER', 'LO', 'BRIGHT', mjd, mjd + 14, tooid]
            )

        return too_list


class TNSAlertList(TooAlertList):
    """ToO input handling from the Transient Name Server (TNS), https://www.wis-tns.org/.
    """
    def __init__(self, too_table):
        super().__init__(too_table)

    def generate_too_list(self):
        """Find unique alerts not in the ToO database.

        Returns
        -------
        too_list : tuple or list
            Per-alert info needed for the ToO ledger.
        """

        too_list = []

        for entry in self.too_table:
            # Extract the relevant data for the ToO alert.
            objid = entry['Name']
            coord = SkyCoord(entry['RA'], entry['DEC'], unit=(u.hourangle, u.degree), frame='icrs')
            ra, dec = coord.ra.to_value('deg'), coord.dec.to_value('deg')
            date = entry['Discovery Date (UT)']
            t = Time(date, format='iso', scale='utc')
            mjd = int(np.floor(t.mjd))
            instrument = entry['Disc. Instrument/s']

            # Enter data into the DB. 
            tooid = self.toodb.add_alert(objid, instrument, date, mjd, ra, dec)
            if tooid == 0:
                continue

            # Accumulate data for output:
            # RA, DEC, PMRA, PMDEC, EPOCH, CHECKER, TYPE, PRIO, PROG, MJD_START, MJD_STOP, TOO_ID
            too_list.append(
                [ra, dec, 0., 0., 2000.0, 'SB/AP', 'FIBER', 'LO', 'BRIGHT', t.mjd, t.mjd + 14, tooid]
            )

        return too_list

