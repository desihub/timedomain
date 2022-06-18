#!/usr/bin/env python

from astropy import units as u
from astropy.table import Table, vstack, hstack, join
import psycopg2
from tqdm import tqdm
import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def get_pv_targets(pvtargfile):
    """Get PV targets file. Join with SGA data if needed.
    
    Parameters
    ----------
    pvtargfile : str
        Absolute or relative path to PV targets FITS file.
    
    Returns
    -------
    ttarg : astropy.Table
        Targets table, including underlying galaxy data.
    """
    ttarg = Table.read(pvtargfile)
    if 'SGA_ID' not in ttarg.columns:
        if 'sv' in pvtargfile:
            pvfullfile = pvtargfile.replace('sv', 'full')
        else:
            pvfullfile = pvtargfile.replace('.fits', '_full.fits')
        tfull = Table.read(pvfullfile)
        ttarg = join(ttarg, tfull['OBJID','BRICKID','BRICKNAME','SGA_ID'],
                     keys=['OBJID','BRICKID','BRICKNAME'])
    return ttarg


def match_targets(pvtargtab, redux='daily', search='healpix'):
    """Match PV targets against the redshift DB for a particular spectroscopic
    reduction.

    Parameters
    ----------
    pvtargtab : astropy.Table
        Table of PV target info; need the RA, DEC, PVTYPE, and SGA_ID fields.
    redux : str
        Spectroscopic reduction: 'daily', 'everest', 'fuji', 'guadalupe', ...
    search : str
        'healpix' to search HEALPix tables, 'tiles' to search tiles tables.

    Returns
    -------
    desi_targets : astropy.Table
        Joined table of DESI redshifts and PV targets for all matches.
    """
    # Accumulate data in this table.
    desi_targets = None

    try:
        db = psycopg2.connect(host='decatdb.lbl.gov', database='desidb', user='desi')
        cursor = db.cursor()
        # cursor.execute('SET search_path TO da, public;')

        # Loop over all TNS alerts and perform a coordinate match with DESI
        # observations.
        N = len(pvtargtab)
        n = 0
        with tqdm(total=N) as progress_bar:

            for i, obj in enumerate(pvtargtab):
                ra, dec = obj['RA'], obj['DEC']

                # Enable search in HEALPix tables.
                if search == 'healpix':
                    query = 'SELECT f.targetid,f.target_ra,f.target_dec,h.healpix,h.survey,r.z,r.zerr,r.zwarn,r.deltachi2,h.filename\n' \
                            f'FROM {redux}.healpix_fibermap f\n' \
                            f'INNER JOIN {redux}.healpix h ON f.healpix_id=h.id\n' \
                            f'INNER JOIN {redux}.healpix_redshifts r ON r.healpix_id=h.id AND r.targetid=f.targetid\n' \
                            f'WHERE q3c_radial_query( f.target_ra, f.target_dec, {ra}, {dec}, 1./3600. );'

                    colnames = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'HEALPIX', 'SURVEY', 'Z', 'ZERR', 'ZWARN', 'DELTACHI2', 'FILENAME']
                # Enable search in tiles tables.
                elif search == 'tiles':
                    query = 'SELECT f.targetid,f.target_ra,f.target_dec,c.tileid,c.night,r.z,r.zerr,r.zwarn,r.deltachi2,c.filename\n' \
                            f'FROM {redux}.tiles_fibermap f\n' \
                            f'INNER JOIN {redux}.cumulative_tiles c ON f.cumultile_id=c.id\n' \
                            f'INNER JOIN {redux}.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid\n' \
                            f'WHERE q3c_radial_query( f.target_ra, f.target_dec, {ra}, {dec}, 1./3600. );'
                    colnames = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'TILEID', 'NIGHT', 'Z', 'ZERR', 'ZWARN', 'DELTACHI2', 'FILENAME']
                else:
                    raise ValueError(f'Search {search} not recognized; use "healpix" or "tiles."')

                cursor.execute(query)
                rows = cursor.fetchall()

                if rows:
                    # Convert postgresql row output to an astropy Table.
                    data = Table(list(map(list, zip(*rows))),
                                 names=colnames)

                    # hstack the postgresql rows with the PV target info.
                    # The following vstack loop ensures every row gets a match.
                    pv_data = obj
                    if len(data) > 1:
                        for j in range(1, len(data)):
                            pv_data = vstack([pv_data, obj])
                    data = hstack([data, pv_data['PVTYPE', 'SGA_ID', 'RA', 'DEC']])

                    # Accumulate matched targets.
                    if desi_targets is None:
                        desi_targets = data
                    else:
                        desi_targets = vstack([desi_targets, data], join_type='outer')

                if (i+1) % 50 == 0:
                    progress_bar.update(50)
                    n += 50

            if n < N:
                progress_bar.update(N - n)

    except (Exception, psycopg2.Error) as error:
        print(error)
    finally:
        if db is not None:
            db.close()

    return desi_targets


if __name__ == '__main__':
    p = ArgumentParser(description='Redshift DB target search',
                       formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('infiles', nargs='+',
                   help='')
    p.add_argument('-r', '--redux', default='daily',
                   choices=['daily', 'everest', 'fuji', 'guadalupe'],
                   help='Spectroscopic reduction')
    p.add_argument('-s', '--search', default='healpix',
                   choices=['healpix', 'tiles'],
                   help='Search healpix or tiles tables')
    p.add_argument('-o', '--output', default='desi_targets.fits',
                   help='Output file name')

    args = p.parse_args()

    output = None
    for infile in args.infiles:
        data = get_pv_targets(infile)
        print(f'Matching targets from {infile}...')
        matched_data = match_targets(data, redux=args.redux, search=args.search)

        if output is None:
            output = matched_data
        else:
            output = vstack([output, matched_data])

    output.write(args.output, overwrite=True)

