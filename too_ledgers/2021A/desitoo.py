"""ToO ledger code for DESI.
"""

import os

from bs4 import BeautifulSoup
import requests
import json

import pandas as pd

from astropy.time import Time

from abc import ABC, abstractmethod


class ToOLedgerMaker(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def write_too_ledger(self, filename, checker, overwrite=False, verbose=False):
        pass


class DECamLedgerMaker(ToOLedgerMaker):

    def __init__(self, url, usecached=True):
        """Download reduced DECam transient data from Texas A&M.
        Cache the data to avoid lengthy and expensive downloads.
        
        Parameters
        ----------
        url : str
            URL for accessing the data.
        usecached : bool
            If False, download new data and overwrite cached data.
            
        Returns
        -------
        decam_transients : pandas.DataFrame
            Table of transient data.
        """
        folders = url.split('/')
        thedate = folders[-1] if len(folders[-1]) > 0 else folders[-2]
        outfile = '{}.csv'.format(thedate)
        
        if os.path.exists(outfile) and usecached:
            # Access cached data.
            decam_transients = pd.read_csv(outfile)
        else:
            # Download the DECam data index.
            # A try/except is needed because the datahub SSL certificate isn't playing well with URL requests.
            try:
                decam_dets = requests.get(url).text
            except:
                requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)
                decam_dets = requests.get(url, verify=False).text
                
            # Convert transient index page into scrapable data using BeautifulSoup.
            soup = BeautifulSoup(decam_dets)
            
            # Loop through transient object summary JSON files indexed in the main transient page.
            # Download the JSONs and dump the info into a Pandas table.
            decam_transients = None
            j = 0

            for a in soup.find_all('a', href=True):
                if 'object-summary.json' in a:
                    link = a['href'].replace('./', '')
                    summary_url  = url + link        
                    summary_text = requests.get(summary_url, verify=False).text
                    summary_data = json.loads(summary_text)

                    j += 1
                    print('Accessing {:3d}  {}'.format(j, summary_url))

                    if decam_transients is None:
                        decam_transients = pd.DataFrame(summary_data, index=[0])
                    else:
                        decam_transients = pd.concat([decam_transients, pd.DataFrame(summary_data, index=[0])])
                        
            # Cache the data for future access.
            print('Saving output to {}'.format(outfile))
            decam_transients.to_csv(outfile, index=False)
            
        self.too_table = decam_transients

    def write_too_ledger(self, filename, checker, overwrite=False, verbose=False):
        """Write ToO ledger in the ECSV format specified by Adam Meyers.
        These can be passed to fiberassign for secondary targeting.

        Parameters
        ----------
        filename : str
            Output filename of the ledger (can be an absolute path).
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
            for idx, row in self.too_table.iterrows():
                # Ledger format:
                # datatype:
                # - {name: RA, unit: deg, datatype: float64}
                # - {name: DEC, unit: deg, datatype: float64}
                # - {name: PMRA, unit: mas / yr, datatype: float32}
                # - {name: PMDEC, unit: mas / yr, datatype: float32}
                # - {name: REF_EPOCH, unit: yr, datatype: float32}
                # - {name: CHECKER, datatype: string}
                # - {name: TOO_TYPE, datatype: string}
                # - {name: OCLAYER, datatype: string}
                # - {name: MJD_BEGIN, unit: d, datatype: float64}
                # - {name: MJD_END, unit: d, datatype: float64}
                ra = row['RA-OBJECT']
                dec = row['DEC-OBJECT']
                t_disc = Time(row['Discovery-Time'], format='isot', scale='utc')
                t_last = Time(row['Latest-Time'], format='isot', scale='utc')

                mag_found  = row['Discovery-Magnitude']
                mag_latest = row['Latest-Magnitude']
                prog = 'BRIGHT' if mag_latest < 21 else 'DARK'

                epoch = 2000.0

                outf.write('{:<10.6f} {:>10.6f} {:>9.6f} {:>9.6f} {:>7.1f}  {}  FIBER  {:7s} {:>15.8f} {:>15.8f}\n'.format(
                        ra, dec, 0, 0, epoch, checker, prog, t_disc.mjd, t_last.mjd+14))

            if verbose:
                outf.seek(0)
                for line in outf:
                    print(line.strip())
