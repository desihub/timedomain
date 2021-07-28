from astropy.table import Table, vstack, unique, SortedArray
import glob
import time
from datetime import date, timedelta, datetime
import os
import numpy





def readInFile(yyyymmdd):
    #for a in query_arr:
    start_time = time.time()
    #get time duration stamp on process that reads in all spectra files for one night
   
    
    #read in and store in one place all the fibermap information in the spectra files
    dats=[]
    for filename in glob.glob(f"/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative/*/{yyyymmdd}/spectra-*.fits"):
        t = Table.read(filename, format='fits',hdu=1, memmap=True)[['TARGETID','TARGET_RA','TARGET_DEC','TILEID','OBJTYPE','PETAL_LOC','NIGHT','MJD']]
        t=t[t['OBJTYPE']=='TGT']
        dats.append(t)
    dats=vstack(dats, join_type='inner',metadata_conflicts='silent')
    
    
    # group all the observations by TARGET_RA and TARGET_DEC
    # note that this is more reliable than grouping by TARGETID as TARGETID is NOT a unique identifier of RA and DEC
    dats_group = dats.group_by(['TARGET_RA','TARGET_DEC'])
    print("--- read in file:  %s seconds ---" % (time.time() - start_time))
    #get RA_DEC, NIGHT for multiple pairs of tile/petal
    return(dats_group)