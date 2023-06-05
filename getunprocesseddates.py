import time
from datetime import date, timedelta, datetime
import psycopg2
import sqlite3
import os
import numpy

#get array of yyyymmdds to loop through
#getUnprocessedDates.py
def getUnprocessedDates():
    start_time = time.time()
    # Connect to desi.db with POSTGRES to get latest observed and unprocessed yyyymmdd
    f = open('/global/cscratch1/sd/clepart/desi_pg.txt') #what is more valid?
    file = f.read()
    db_name, db_user, db_pwd, db_host = file.split()
    conn = psycopg2.connect(dbname=db_name, user=db_user, password=db_pwd, host=db_host)
    cur = conn.cursor()
    cur.execute("""SELECT DISTINCT yyyymmdd from fibermap_daily WHERE yyyymmdd > 20210604""") #most recent, remove in future
    desi_arr = cur.fetchall()
    cur.close()
    conn.close()
    
    # Open to transients_search.db for latest processed yyyymmdd to do comparison with unprocessed
    filename_conn = "/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db"
    conn = sqlite3.connect(filename_conn)
    trans_arr = conn.execute("""SELECT DISTINCT yyyymmdd from unprocessed_exposures""").fetchall()
    conn.close()

    # Compare yyyymmdd from fibermap_daily with unprocessed_exposures, retaining those in fibermap_daily and not in unprocessed_daily
    night_arr = numpy.setdiff1d(desi_arr, trans_arr)

    print('len(night_arr): ' + str(len(night_arr)))
    print("--- get unprocessed dates took:  %s seconds ---" % (time.time() - start_time))
    
    return(night_arr) 