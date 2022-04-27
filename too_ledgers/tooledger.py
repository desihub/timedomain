"""This module implements TooAlertList handlers to manage a large variety of
possible ToO file formats from different sources, such as TNS, DECam, etc. To
add a new format, create a subclass of TooAlertList and add a format check to
the AlertListFactory to use the new format.
"""

from toodb import TooDB

from abc import ABC, abstractmethod

from astropy.table import Table
from astropy.time import Time

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
        else:
            raise RuntimeError('Unrecognized ToO alert list format.')

factory = AlertListFactory()


class TooLedgerMaker:
    """Generic alert handler to manage ToO lists in different formats."""

    def build_too_ledger(self, too_input_file, too_ledger_file):
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
        alerthandler.write_too_list(too_ledger_file)


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

    def write_too_list(self, outputfile):
        """Append ToO data to an ecsv file for processing by fiberassign.

        Parameters
        ----------
        outputfile : str
            Path to output ecsv file.
        """
        too_list = self.generate_too_list()

        with open(outputfile, 'a') as outf:
            for entry in too_list:
                print(entry)
                ra, dec, pmra, pmdec, epoch, checker, too_type, prio, prog, mjd0, mjd1, too_id = entry
                outf.write('{:<10.6f} {:>10.6f} {:>8.6f} {:>8.6f} {:>6.1f} {} {} {} {} {:>13.8f} {:>13.8f} {}\n'.format(ra, dec, pmra, pmdec, epoch, checker, too_type, prio, prog, mjd0, mjd1, too_id))


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
#    if 'XWIN_WORLD' in data.columns:

