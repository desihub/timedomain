import numpy as np
from datetime import datetime
from astropy.io import ascii
from astropy.time import Time
from argparse import ArgumentParser

from antares_client.search import search, download

def build_query(ra0, dec0, fov, date):
    """Generate a query (a Python dictionary) to submit to the ANTARES client.

    Parameters
    ----------
    ra0 : float or None
        Central RA for object search, in deg.
    dec0 : float or None
        Central declination for object search, in deg.
    fov : float
        Side length of box for search, in deg.
    date : str
        Start date for search; format is YYYY-MM-DD.

    Returns
    -------
    query : dict
        An ElasticSearch dictionary.
    """

    # Build up the query.
    query = { 'query': { 'bool': { 'must': [] } } }

    # desi_candidate_test data stream:
    # snfilter_last_proc_status should have a string like "Locus has two or
    # more detections and is in DESI brightness range. Triggering."
    query['query']['bool']['must'].append(
        { 'match':{ 'properties.snfilter_last_proc_status': '*DESI*' } })

    # Set up the declination search.
    if dec0 is not None:
        ddec = 0.5 * fov
        # dra / cos(dec) ensures an equal-area search rectangle.
        dra = 0.5*fov / np.cos(np.radians(dec0))

        query['query']['bool']['must'].append(
            {'range': {'dec':{ 'gte':dec0-ddec, 'lte':dec0+ddec, } } })
    else:
        dra = 0.5*fov

    # Set up the RA search.
    if ra0 is not None:
        query['query']['bool']['must'].append(
            {'range': {'ra':{ 'gte':(ra0-dra)%360., 'lte':(ra0+dra)%360., } } })

    # Set up the cumulative date search.
    if date is not None:
        tobs = Time(date).mjd
        query['query']['bool']['must'].append(
            {'range': {'mjd':{ 'gte':tobs, } } })

    return query

if __name__ == '__main__':

    today = datetime.today()

    parser = ArgumentParser(description='Client API to query ANTARES alert DB')
    parser.add_argument('--ra', default=None, type=float,
                        help='RA (J2000), in deg')
    parser.add_argument('--dec', default=None, type=float,
                        help='Dec (J2000), in deg')
    parser.add_argument('--tobs', default=datetime.today().strftime('%Y-%m-%d'),
                        help='Obs date [YYYY-MM-DD]')
    args = parser.parse_args()

    # Create query dict for ANTARES stream search.
    query = build_query(ra0=args.ra, dec0=args.dec, fov=3.2, date=args.tobs)
    print(query)

    #result_set = search(query)
    #print(result_set)

    outfile = 'results_antares'
    if args.ra is not None:
        outfile = '{}_ra{:03.1f}'.format(outfile, args.ra)
    if args.dec is not None:
        outfile = '{}_dec{:03.1f}'.format(outfile, args.dec)
    if args.tobs is not None:
        outfile = '{}_{}'.format(outfile, args.tobs)
    outfile += '.csv'

    result_set = download(query, outfile, output_format='csv', decompress=True)

