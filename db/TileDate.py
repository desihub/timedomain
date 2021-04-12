import os
import fitsio

from desispec.io import read_spectra, write_spectra
from timedomain.sp_utils import *
from timedomain.fs_utils import *
from timedomain.iterators import *
from timedomain.filters import *

import matplotlib.pyplot as plt
import numpy
import sqlite3
import glob
from astropy.io import ascii

class TileDate:

    filename = "/global/cfs/cdirs/desi/science/td/db/secondary.db"
    
    # init method or constructor   
    def __init__(self):
        self.con = sqlite3.connect(TileDate.filename)
        self.cur = self.con.cursor()
        
    @staticmethod
    def type_py2sql(name):
        if name.kind =='b' or name.kind: return 'INTEGER'
        if name.kind =='f': return 'REAL'
        raise Exception("Boy are we screwed")
        
    def load_targets(self, createTable=False):
        dir_root = "/global/cfs/cdirs/desi/target/secondary/sv3/outdata/0.57.0/"
        lunations = ["bright", "dark"]
        table_name = "targets"
        
        if createTable:
            # draw one random file to pull information about the format and create the table
            command = f"CREATE TABLE {table_name} (LUNATAION TEXT, PROGRAM TEXT, TARGETID INTEGER, RA REAL, DEC REAL, UNIQUE(TARGETID))"
            self.cur.execute(command)

        # Fill in the table
        for lunation in lunations:
            for root, dirs, files in os.walk(os.path.join(dir_root,lunation)):
                for file in files:
                    if file.endswith('.fits'):
                        filename = os.path.join(root,file)
                        print(file)
                        program = file[:-5]
                        fits = fitsio.FITS(filename)
                        data = fits[1]['TARGETID','RA','DEC']
                        for t in zip(data['TARGETID'].read(),data['RA'].read(),data['DEC'].read()):
                            converted_t = ["'{}'".format(str(element)) for element in t]
                            command = "INSERT INTO {} VALUES ('{}', '{}', {});".format(table_name,lunation,program, ",".join(converted_t))
                            try:
                                self.cur.execute(command)
                            except sqlite3.IntegrityError:
                                pass
                                
#         if createTable:
#             # draw one random file to pull information about the format and create the table
#             command = f"CREATE TABLE {table_name} (LUNATAION TEXT, PROGRAM TEXT"
#             filename = glob.glob(os.path.join(dir_root,lunations[0]+'/*.fits'))[0]
#             fits = fitsio.FITS(filename)
#             for col in fits['TARGETS'].get_colnames():
#                 command += ", {} {}".format(col,DB.type_py2sql(fits['TARGETS'].get_rec_dtype()[0][col])) 
#             command +=" UNIQUE (TARGETID))"
#             self.cur.execute(command)

#         # Fill in the table
#         for lunation in lunations:
#             for root, dirs, files in os.walk(os.path.join(dir_root,lunation)):
#                 for file in files:
#                     if file.endswith('.fits'):
#                         filename = os.path.join(root,file)
#                         program = file[:-5]
#                         fits = fitsio.FITS(filename)
#                         for t in fits['TARGETS']:
#                             converted_t = ["'{}'".format(str(element)) for element in t]
#                             command = "INSERT INTO {} VALUES ('{}', '{}', {});".format(table_name,lunation,program, ",".join(converted_t))
#                             try:
#                                 self.cur.execute(command)
#                             except sqlite3.IntegrityError:
#                                 print("couldn't add targetid twice")
        self.con.commit()
                            
    def load_tiles(self, createTable=False):

        table_name = "tiles"
        if createTable:
            command = f"CREATE TABLE {table_name} (TARGETID, TILEID, UNIQUE (TARGETID, TILEID))"
            self.cur.execute(command)      
            
#         # SV 1
#         filename = "/global/cfs/cdirs/desi/spectro/redux/daily/tiles.csv"
#         fbase = "/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk"
#         data = ascii.read(filename, include_names=["TILEID","SURVEY"])  
#         for d in data:
#             if d[1]=='sv1':
#                 ptile = d[0].astype('str')
#                 ptile = ptile.zfill(6)
#                 ffile = os.path.join(fbase,ptile[:3],f'fiberassign-{ptile}.fits.gz')
#                 fits = fitsio.FITS(ffile)
#                 col = fits['FIBERASSIGN']['TARGETID']
#                 for c in col.read():
#                     command = "INSERT INTO {} VALUES ('{}', '{}');".format(table_name,c,d[0])
#                     try:
#                         self.cur.execute(command)
#                     except sqlite3.IntegrityError:
#                         pass
                    
        # SV2 onward
        dir_root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
        subdirs = ['sv2','sv3','sv3/secondary']
        lunations = ["bright", "dark"]    

        # Fill in the table
        for sd in subdirs:
            for lunation in lunations:
                for root, dirs, files in os.walk(os.path.join(dir_root,sd,lunation)):
                    for file in files:
                        if file.endswith('.ecsv'):
                            filename = os.path.join(root,file)
                            print(filename)
                            data = ascii.read(filename, include_names=["TARGETID","ZTILEID"])

                            for t,z in data:
                                command = "INSERT INTO {} VALUES ('{}', '{}');".format(table_name,t,z)
                                try:
                                    self.cur.execute(command)
                                except sqlite3.IntegrityError:
                                    pass
        self.con.commit()

    def load_daily(self,createTable=False):
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
        
        table_name = "daily"
        
        if createTable:
            command = f"CREATE TABLE {table_name} (TILEID, YYYYMMDD, UNIQUE (TILEID, YYYYMMDD))"
            self.cur.execute(command)   
        
        #find the last date
        ans=self.cur.execute("SELECT MAX(YYYYMMDD) FROM daily")
        maxdate = ans.fetchone()[0]
        if maxdate is None: maxdate = ""
        print(maxdate)
        
        for path in glob.glob(f'{dir_root}/*/*'):
            split = path.split('/')
            tileid=split[-2]; date=split[-1]
            if date >= maxdate and tileid.isnumeric():
                command = "INSERT INTO {} VALUES ('{}', '{}');".format(table_name,tileid,date)
                try:
                    self.cur.execute(command)
                except sqlite3.IntegrityError:
                    print("couldn't add daily twice")
        self.con.commit()
        
    @staticmethod
    def byTARGETID(targetid):
        db=TileDate()
        command = f'''
        SELECT daily.TILEID,  daily.YYYYMMDD
        FROM daily
        INNER JOIN tiles
            ON daily.TILEID = tiles.TILEID
        INNER JOIN targets 
            ON tiles.TARGETID = targets.TARGETID
        WHERE tiles.TARGETID = "{targetid}";
        '''
        return db.cur.execute(command).fetchall()

    @staticmethod
    def byProgram(program):
        db=TileDate()
        command = f'''SELECT targets.PROGRAM, targets.RA, targets.DEC, tiles.TARGETID, daily.TILEID,  daily.YYYYMMDD
        FROM daily
        INNER JOIN tiles
            ON daily.TILEID = tiles.TILEID
        INNER JOIN targets 
            ON tiles.TARGETID = targets.TARGETID
        WHERE targets.PROGRAM LIKE "{program}%"
        ORDER BY
            targets.PROGRAM,
            tiles.TARGETID,
            daily.TILEID,
            daily.YYYYMMDD;
        '''
        return db.cur.execute(command).fetchall()
    
if __name__ == "__main__":
    # running as a cron job on desi16
    db = TileDate()
    db.load_daily()