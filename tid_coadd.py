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

filename_trans = "/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db"

con = sqlite3.connect(filename_trans)
cur = con.cursor()

query = """
select distinct targetid, yyyymmdd, tileid, petal
from processed_daily;
"""

cur.execute(query)
query_arr = cur.fetchall()
cur.close()

coadd_arr=[]
#new array to hold targetids
tid_night=[query_arr[0]]
for i in range(1,len(query_arr)):
    if query_arr[i][0]== tid_night[0][0] and query_arr[i][1]== tid_night[0][1]:
        tid_night.append(query_arr[i])
    else:
        coadd_arr.append(tid_night)
        tid_night=[query_arr[i]]
        if i == len(query_arr)-1:
             coadd_arr.append(tid_night)
                
# for i in range(len(coadd_arr)):
#         print(coadd_arr[i])
                
#datapath = os.path.join(os.environ.get('DESI_REDUX'),'daily','tiles','cumulative')

datapath = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative/"

new_spectra_coadd=[] #do we need an array of new spectra?


for a in coadd_arr:
#loop through targetids (exposures) from processed_daily and create a SpecObj from each
    new_spec=Spectra()
    #new SpecObj per targetid
    ref_spec=Spectra()
    for targetid, yyyymmdd, tileid, petal in a:
    #loop through targetids  
        new_filename = os.path.join(datapath,str(tileid),str(yyyymmdd),f"spectra-{petal}-{tileid}-thru{yyyymmdd}.fits")
        #find file already made? or make file? for that exposure
        hold_spec = read_spectra(new_filename)
        #filename becomes a SpecObj
        new_spec.update(hold_spec.select(targets=[targetid]))
        #update new_spec (which is initially empty) with info from hold_spec based on targetid -> created a SpecObj (new_spec) which holds targetid, tile, yyyymmdd, petal
    new_spectra_coadd.append(new_spec)
        

       

