#import libraries
import time
import psycopg2
import sqlite3
import numpy
import sys


#function to get data from desi databases
def getz(ra, dec):
    start_time = time.time()
    f = open('/global/cfs/cdirs/desi/science/td/secrets/desi_pg.txt') 
    file = f.read()
    db_name, db_user, db_pwd, db_host = file.split()
    db_name = 'desidb'
    conn = psycopg2.connect(dbname=db_name, user=db_user, password=db_pwd, host=db_host)
    cur = conn.cursor()
    curr = """SELECT f.targetid,c.tileid,c.night,r.z, f.target_ra, f.target_dec, r.zerr,r.zwarn,c.filename
FROM daily.tiles_fibermap f
INNER JOIN daily.cumulative_tiles c ON f.cumultile_id=c.id
INNER JOIN daily.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid
WHERE q3c_radial_query( f.target_ra, f.target_dec, {0}, {1}, 10./3600. )""".format(ra, dec)
    cur.execute(curr) 
    desi_arr = cur.fetchall()
    
    if (len(desi_arr) == 0):
        curr = """SELECT f.targetid,c.tileid,c.night,r.z, f.target_ra, f.target_dec, r.zerr,r.zwarn,c.filename
FROM fuji.tiles_fibermap f
INNER JOIN fuji.cumulative_tiles c ON f.cumultile_id=c.id
INNER JOIN fuji.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid
WHERE q3c_radial_query( f.target_ra, f.target_dec, {0}, {1}, 10./3600. )""".format(ra, dec)
        cur.execute(curr) 
        desi_arr = cur.fetchall()

    if (len(desi_arr) == 0):
        curr = """SELECT f.targetid,c.tileid,c.night,r.z, f.target_ra, f.target_dec, r.zerr,r.zwarn,c.filename
FROM guadalupe.tiles_fibermap f
INNER JOIN guadalupe.cumulative_tiles c ON f.cumultile_id=c.id
INNER JOIN guadalupe.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid
WHERE q3c_radial_query( f.target_ra, f.target_dec, {0}, {1}, 10./3600. )""".format(ra, dec)
        cur.execute(curr) 
        desi_arr = cur.fetchall()
        
    cur.close()
    conn.close()
    return desi_arr
    

#helper function to create output and save in a text file
def helper(input_info):
    #helper_info is a dic of ZTF id, input_ra, input_dec
    
    with open('output_test.txt', 'w') as f:
        original_stdout = sys.stdout
        sys.stdout = f # Change the standard output to the file we created.
        ZTFids = input_info['id']
        ras = input_info['RA']
        decs = input_info['DEC']
        print("ZTF_ID" + "\t \t Input_RA" + "\t Input_DEC" + "\t DESI_TID" + "\t     DESI_Night" + "\t \t DESI_Z" + "\t \t \t DESI_RA" + "\t \t DESI_DEC")
        for i in range(len(ras)):
            ra = ras[i]
            dec = decs[i]
            ztfid = ZTFids[i]

            zinfo = getz(ra, dec)

            if (len(zinfo) == 0):
                print(str(ztfid) + "\t" +  str(ra) + "\t" + str(dec) + "\t No_Info." )
            else:
                for i in range(len(zinfo)):
                    tid = zinfo[i][0]
                    night = zinfo[i][2]
                    z = zinfo[i][3]
                    zerr = zinfo[i][6]
                    desi_ra = zinfo[i][4]
                    desi_dec =zinfo[i][5]
                    if (zerr != 0): 
                        print(str(ztfid) + "\t" +  str(ra) + "\t" + str(dec) + "\t zerr != 0" )
                        break
                    if (z > 0.15):
                        print(str(ztfid) + "\t" +  str(ra) + "\t" + str(dec) + "\t" + str(tid) + "\t"  + str(night) + "\t" + str(z) + "\t" + str(desi_ra)+ "\t" + str(desi_dec))

                    elif (z < 0.15):
                        print(str(ztfid) + "\t" +  str(ra) + "\t" + str(dec) + "\t" + str(tid) + "\t"  + str(night) + "\t N/A" )
    sys.stdout = original_stdout # Reset the standard output to its original value
    

#extract the data from ZTF text file
data = {'id': [], 'RA': [], 'DEC': []}
i = 1
with open('test.txt') as infile:
    for line in infile:
        if (i == 1):
            i += 1
            continue # to not let column names get into the data
        data['id'].append(line.split()[0])
        data['RA'].append(line.split()[1])
        data['DEC'].append(line.split()[2])


#output data into output text file
helper(data)
