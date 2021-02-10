""" Download entries from TNS into a sqlite3 database

Eventually this could be run on a regular basis

Anna Ho, 10 Feb 2021
"""

# Import packages
import numpy as np
import pandas as pd
import io
import sqlite3
import requests
from tns_watcher_utils import *


# Create a table to ingest tns sources
# In the future this could instead add rows to the table, instead of creating a new one from scratch
conn = sqlite3.connect('tns.db')
c = conn.cursor()
c.execute('''CREATE TABLE tns_sources 
             (ID text, RA real, Dec real)''')


# Download all sources on TNS from the last three months
base_url = "https://www.wis-tns.org/search?&discovered_period_value=3&format=csv"
url = f"{base_url}&page=0"
csv_data = requests.get(url, verify=False).content
data = pd.read_csv(io.StringIO(csv_data.decode("utf-8")), error_bad_lines=False)
npages = int(data['ID'].max()/len(data))+1

# Have to loop through TNS page by page
for page in np.arange(npages):
    url = f"{base_url}&page={page}"
    csv_data = requests.get(url, verify=False).content
    data = pd.read_csv(io.StringIO(csv_data.decode("utf-8")), error_bad_lines=False)

    # Now loop through all the objects on a single page
    # (there might be a way to ingest all at once...)
    for ii in np.arange(len(data)):
        insert_id = data['Name'][ii]
        insert_ra = data['RA'][ii]
        insert_dec = data['DEC'][ii]

        # Convert 
        # store RA and Dec as decimal degrees
        ra_rad,dec_rad = radec_str2rad(insert_ra, insert_dec)
        ra_deg = ra_rad*180/np.pi
        dec_deg = dec_rad*180/np.pi
        
        # Add row to the table
        c.execute("INSERT INTO tns_sources (ID,RA,Dec) VALUES ('%s',%s,%s)" %(insert_id,ra_deg,dec_deg))
        conn.commit()

conn.close()
