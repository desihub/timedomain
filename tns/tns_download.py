#!/usr/bin/env python
"""Access TNS using a URL GET. Save the data to CSV.

Adapted from the tns_watcher script by Dima Duev,
https://github.com/dmitryduev/kowalski/blob/master/kowalski/tns_watcher.py
which was ported here by Anna Ho.

Example usage:
> python tns_download.py --classified -n 3

This will download 3 months of alerts classified as SNe and save them to CSV.
"""

import argparse
import requests
import pandas as pd
import io


def download_data(args):
    """Download data from the Transient Name Server.

    Parameters
    ----------
    args : dict
        Dictionary of key-value pairs to pass to the TNS URL GET query.

    Returns
    -------
    tns_data : pandas.DataFrame
        Accumulated data tables from TNS.
    """

    # Build up the URL GET.
    base_url = 'https://www.wis-tns.org/search?format=csv'
    for k, v in args.items():
        if v is None:
            continue
        base_url = f'{base_url}&{k}={v}'

    # User-agent needed for request to go through.
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.55 Safari/537.36'
    }

    # Start the request. Set up a loop through pages returned by TNS.
    # There is a maximum of 500 results per page.
    loop = True
    page = 0
    tns_data = None
    npage = 500
    
    while loop:
        url = f'{base_url}&num_page={npage}&page={page}'
        print(url)
        f = requests.get(url, headers=headers)
        data = pd.read_csv(io.StringIO(f.content.decode('utf-8')), error_bad_lines=False)
    
        # We're done if we have no data...
        if len(data) == 0:
            break
    
        # ...or if the data on this page is less than npage.
        loop = len(data) >= npage
    
        # Accumulate data from each page.
        if tns_data is None:
            tns_data = data
        else:
            tns_data = pd.concat([tns_data, data])
    
        page += 1

    return tns_data


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='TNS data access script',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--classified', dest='classified_sne',
                   default=None, const=1, nargs='?', type=int,
                   help='Require classification of alerts.')
    p.add_argument('-n', '--nperiod', dest='discovered_period_value', default=2,
                   help='Number of periods (days, weeks, months) to check for alerts.')
    p.add_argument('-u', '--units', dest='discovered_period_units', default='months',
                   choices=['days', 'weeks', 'months'],
                   help='Number of periods (days, weeks, months) to check for alerts.')
    p.add_argument('-o', '--output', dest='output', default='tns_search.csv',
                   help='Output CSV file from TNS query.')

    args = p.parse_args()

    # Clean up output arg so it's not passed to download_data.
    output = args.output
    del args.output

    # Access data and save to output.
    tns_data = download_data(vars(args))
    tns_data.to_csv(output)

    print(f'Saved {len(tns_data)} alerts to {output}')
