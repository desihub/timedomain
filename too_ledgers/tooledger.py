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
from astropy.coordinates import SkyCoord, angular_separation

from datetime import datetime

import numpy as np


class DESICalibField :
    """A simple class that stores the name, location, and observing dates of
    DESI calibration fields. See https://desi.lbl.gov/trac/wiki/SurveyOps/CalibrationFields.
    """

    def __init__(self, name, ra, dec, dates):
        """Initialize a calibration field.

        Parameters
        ----------
        name : str
            Name of the field.
        ra : float
            Field central RA.
        dec : float
            Field central declination.
        dates : str
            Months when the field is observed.
        """
        self.name = name
        self.ra = ra
        self.dec = dec
        self.months = dates.split()

    def in_fov(self, ra, dec, sepcut=1.6):
        """Return true if a sky position is within the calibration field.

        Parameters
        ----------
        ra : float
            RA of input coordinates.
        dec : float
            Declination of input coordinates.
        sepcut : float
            Separation angle cut for FOV selection.

        Returns
        -------
        infield : bool
            True if target is inside a calibration field.
        """
        sep = angular_separation(self.ra*u.deg, self.dec*u.deg, ra*u.deg, dec*u.deg)
        infield = sep < sepcut*u.deg
        return infield

    def in_obstime(self, obstime):
        """Return true if an observation is requested during the coverage of
        this calibration field.

        Parameters
        ----------
        obstime : int, float, str, or datetime
            If float, expect date in MJD. If str, expect month in standard
            short notation (Jan, Feb, Mar, etc.).

        Returns
        -------
        inmonth : bool
            True if date is one of the months observed for this field.
        """
        month = None

        if isinstance(obstime, (int, float, np.floating)):
            month = Time(obstime, format='mjd').to_datetime().strftime('%b')
        elif isinstance(obstime, str):
            month = obstime
        elif isinstance(obstime, datetime):
            month = obstime.strftime('%b')
        else:
            raise TypeError('obstime not type float, str, or datetime')

        return month in self.months


def in_calibration_field(ra, dec, mjd0, mjd1):
    """Return true if a sky position is within a calibration field.

    Parameters
    ----------
    ra : float
        RA of input coordinates.
    dec : float
        Declination of input coordinates.
    mjd0 : float
        Start of observation, in MJD.
    mjd1 : float
        End of observation, in MJD.

    Returns
    -------
    infield : bool
        True if target is inside a calibration field.
    field : str or None
        Field name if target is inside a calibration field, None otherwise.
    """
    desi_calib_fields = [
        DESICalibField('COSMOS',  150.1, 2.182, 'Dec Jan Feb'),
        DESICalibField('M-BHB 1', 203.5, 17.5, 'Mar Apr'),
        DESICalibField('GAMA 15', 215.7, -0.7, 'May Jun Jul'),
        DESICalibField('XMM LSS',  35.7, -4.75, 'Aug Sep Oct Nov')
    ]

    for calfield in desi_calib_fields:
        # DESI field angular radius is 1.6 deg but the calibration tiles are
        # also dithered by up to 1 degree from the center of the field, so
        # perform a search radius of 2.6 degrees.
        fov_cut = 2.6 # degrees

        if calfield.in_fov(ra, dec, fov_cut):
            
            # Check if observation overlaps with this tile at a given time.
            in_obstime = calfield.in_obstime(mjd0) or calfield.in_obstime(mjd1)
            if in_obstime:
                return in_obstime, calfield.name

    return False, None


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

        if 'Candidate' in data.columns and 'Pipeline' in data.columns:
            # DECam ToO format from RK pipeline.
            return DECamAlertListRK(data)
        elif 'XWIN_WORLD' in data.columns:
            # DECam ToO format from TAMU pipeline.
            return DECamAlertListTAMU(data)
        elif 'objid' in data.columns and 'skycoord_obj' in data.columns and 'ra_obj' in data.columns and 'dec_obj' in data.columns:
            # DECam ToO alert list from Lei Hu's pipeline.
            return DECamAlertListLH(data)
        elif 'Field' in data.columns and 'PROGRAM' in data.columns:
            return SMBBHAlertList(data)
        elif 'Reporting Group/s' in data.columns and 'Discovery Data Source/s' in data.columns:
            # TNS ToO format.
            return TNSAlertList(data)
        else:
            raise RuntimeError('Unrecognized ToO alert list format.')

factory = AlertListFactory()


class TooLedgerMaker:
    """Generic alert handler to manage ToO lists in different formats."""

    def build_too_ledger(self, too_input_file, too_ledger_file, add_header=False, verbose=False):
        """Open a ToO input file and build an alert ledger.

        Parameters
        ----------
        too_input_file : str
            Path to ToO input file, in a format readable by astropy.table.
        too_ledger_file : str
            ToO ledger, in a format appropriate for fiberassign.
        add_header : bool
            Add header to output file.
        verbose : bool
            Enable verbose output.
        """
        # Alert factory decides on the file format and acts appropriately.
        alerthandler = factory.get_alert_handler(too_input_file)

        # Appropriate alert handler writes the output.
        alerthandler.write_too_list(too_ledger_file, add_header, verbose)


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

    def write_too_list(self, outputfile, add_header=False, verbose=False):
        """Append ToO data to an ecsv file for processing by fiberassign.

        Parameters
        ----------
        outputfile : str
            Path to output ecsv file.
        add_header : bool
            Add header to output file.
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
            if add_header:
                # Only write the header to an empty file.
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


class DECamAlertListRK(TooAlertList):
    """ToO input handling from DECam using Rob Knop's pipeline.
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
            objid, ra, dec = entry['Candidate'], entry['ra'], entry['dec']
            date = entry['Min_Date']
            t = Time(date, format='isot', scale='utc')
            mjd = int(np.floor(t.mjd))

            # Enter data into the DB. 
            tooid = self.toodb.add_alert(objid, 'DECam', date, mjd, ra, dec)
            if tooid == 0:
                continue

            # Compute the observation window. If the discovery date is old
            # enough that the window will be <9 days, starting today, move the
            # window up.
            now = Time.now().mjd
            dt = now - mjd
            if dt >= 5:
                print(f'Shifting time window for alert {tooid} on {mjd}')
                mjd0, mjd1 = now, now+14
            else:
                mjd0, mjd1 = mjd, mjd+14

            # Accumulate data for output:
            # RA, DEC, PMRA, PMDEC, EPOCH, CHECKER, TYPE, PRIO, PROG, MJD_START, MJD_STOP, TOO_ID
            too_list.append(
                [ra, dec, 0., 0., 2000.0, 'SB/AP', 'TILE', 'HI', 'BRIGHT', mjd0, mjd1, tooid]
            )

        return too_list


class DECamAlertListTAMU(TooAlertList):
    """ToO input handling from DECam/DESIRT using the TAMU pipeline.
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

            # Compute the observation window. If the discovery date is old
            # enough that the window will be <9 days, starting today, move the
            # window up.
            now = Time.now().mjd
            dt = now - mjd
            if dt >= 5:
                print(f'Shifting time window for alert {tooid} on {mjd}')
                mjd0, mjd1 = now, now+14
            else:
                mjd0, mjd1 = mjd, mjd+14

            # Accumulate data for output:
            # RA, DEC, PMRA, PMDEC, EPOCH, CHECKER, TYPE, PRIO, PROG, MJD_START, MJD_STOP, TOO_ID
            in_cal, calname = in_calibration_field(ra, dec, mjd0, mjd1)
            if in_cal:
                # RA, Dec, time puts this observation in a DESI calibration
                # field; set it up for TILE fiberassignment.
                too_list.append(
                    [ra, dec, 0., 0., 2000.0, 'SB/AP', 'TILE', 'HI', 'BRIGHT', mjd0, mjd1, tooid]
                )
            else:
                # Normal observation: FIBER mode, LO priority.
                too_list.append(
                    [ra, dec, 0., 0., 2000.0, 'SB/AP', 'FIBER', 'LO', 'BRIGHT', mjd0, mjd1, tooid]
                )

        return too_list


class DECamAlertListLH(TooAlertList):
    """ToO input handling from DECam/DESIRT using Lei Hu's pipeline.
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
            objid, ra, dec = entry['objid'], entry['ra_obj'], entry['dec_obj']
            date = entry['date_first_alert']
            t = Time(f'{date}', format='isot', scale='utc')
            mjd = int(np.floor(t.mjd))

            # Enter data into the DB. 
            tooid = self.toodb.add_alert(objid, 'DECam', date, mjd, ra, dec)
            if tooid == 0:
                continue

            # Compute the observation window. If the discovery date is old
            # enough that the window will be <9 days, starting today, move the
            # window up.
            now = Time.now().mjd
            dt = now - mjd
            if dt >= 5:
                print(f'Shifting time window for alert {tooid} on {mjd}')
                mjd0, mjd1 = now, now+14
            else:
                mjd0, mjd1 = mjd, mjd+14

            # Accumulate data for output:
            # RA, DEC, PMRA, PMDEC, EPOCH, CHECKER, TYPE, PRIO, PROG, MJD_START, MJD_STOP, TOO_ID
            in_cal, calname = in_calibration_field(ra, dec, mjd0, mjd1)
            if in_cal:
                # RA, Dec, time puts this observation in a DESI calibration
                # field; set it up for TILE fiberassignment.
                too_list.append(
                    [ra, dec, 0., 0., 2000.0, 'SB/AP', 'TILE', 'HI', 'BRIGHT', mjd0, mjd1, tooid]
                )
            else:
                # Normal observation: FIBER mode, LO priority.
                too_list.append(
                    [ra, dec, 0., 0., 2000.0, 'SB/AP', 'FIBER', 'LO', 'BRIGHT', mjd0, mjd1, tooid]
                )

        return too_list


class SMBBHAlertList(TooAlertList):
    """ToO input handling from DECam using Rob Knop's pipeline.
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
            objid, ra, dec = entry['Id'], entry['RA'], entry['Dec']
            instrum = entry['Survey']
            t = Time.now()
            date = t.to_datetime().isoformat(sep='T', timespec='seconds')
            mjd = int(np.floor(t.mjd))

            print(date)

            # Enter data into the DB. 
            tooid = self.toodb.add_alert(objid, instrum, date, mjd, ra, dec)
            if tooid == 0:
                continue

            # Compute the observation window: now to now + 4 yr.
            mjd0, mjd1 = mjd, mjd + 4*365

            # Accumulate data for output:
            # RA, DEC, PMRA, PMDEC, EPOCH, CHECKER, TYPE, PRIO, PROG, MJD_START, MJD_STOP, TOO_ID
            too_list.append(
                [ra, dec, 0., 0., 2000.0, 'SMBBH', 'TILE', 'HI', 'BRIGHT', mjd0, mjd1, tooid]
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

