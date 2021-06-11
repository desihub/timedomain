"""

Given a date and a targetid, return all pairs of spectra from that date and preceeding dates

Default is to compare relative to the most recent release, not necessarily the "subdir"

Internally iterates over:

- input (#tile we don't care which tileid the spectrum came from) targetid/date pairs
  -  preceding dates that contain the targetid
     - panels that exist in both dates

"""
"""
To make a new (reference) coadd, you collect the relevent spectra into a single Spectra object and run the function coadd.
To make the difference spectrum, put the new and -reference spectrum into a Spectra object and run the function coadd.  This should produce the subtraction.

return a pair of spectra (new/old) as a Spectra object

create new_spectra_obj from new_spectra with targetid/(new)date and ref_spectra with targetid/(ref > 30 days)date
coadd(new_spectra_obj)


what are we iterating over? targetid/date pairs found in processed_daily
why are we iterating? so that we can get pairs of spectra of tid/new_date and tid/ref_date (how far back to reference?) for all ref_dates
return this pair of spectra as a Spectra object (which has a fits file? of its own?)
move on to next tid.

with this return value where do we store or coadd or ??

what is the relationship between what we return in tid_coadd.py and read_spectra(filename)? How do I query a spectra obj?
"""

new_spectra_coadd[dict]
for targetid in new_spectra_coadd:
    new_spectra_obj = Spectra()
    new_spectra_obj.update(new_spectra_coadd) 
    new_spectra_obj.update(ref_spectra.select(target=[targetid],date=[i]))
return(new_spectra_obj)
#move on to another ref_date
increment date i++

import sqlite3
import os
import time
from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra
from desispec.coaddition import coadd
from desiutil.log import get_logger, DEBUG
log = get_logger()
log.setLevel(DEBUG)

start_time=time.time()

filename = "/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db"

con = sqlite3.connect(filename)
cur = con.cursor()

for i in coadd_arr[]:
    #create a Spectra object to contain a pair
    for nos. of 
    spectra_pair=Spectra()
    
    
datapath = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative/"

paired_spectra=[] #do we need an array of paired spectra?
for a in iter_arr[]:
    i=100
    spectra=Spectra()
    for targetid, yyyymmdd, tileid, petal in a:
        filename = os.path.join(datapath,str(tileid),str(yyyymmdd),f"spectra-{petal}-{tileid}-thru{yyyymmdd}.fits")
        spec = read_spectra(filename)
        spectra.update(spec.select(targets=[targetid]), date=[yyyymmdd-i]) #What to do about Feb., and months with 31 where preceeding month does not have 31 days? 
    if spectra.num_spectra()>1:
        coadd(spectra)
    paired_spectra.append(spectra)
    
#iterator
    
    return((read_spectra(filename) , read_spectra(filename2)))

