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
#         self.cur = self.con.cursor()
        
    @staticmethod
    def type_py2sql(name):
        if name.kind =='b' or name.kind: return 'INTEGER'
        if name.kind =='f': return 'REAL'
        raise Exception("Boy are we screwed")
        
    @staticmethod
    def update():
        db = TileDate()
        db.load_daily(createTable=False)
        db.load_tiles(createTable=False, latestOnly=True)
        
    def load_targets(self, createTable=False):

        table_name = "targets"        
        if createTable:
            # draw one random file to pull information about the format and create the table
            command = f"CREATE TABLE {table_name} (LUNATAION TEXT, PROGRAM TEXT, TARGETID INTEGER, RA REAL, DEC REAL, UNIQUE(TARGETID));"
            self.cur.execute(command)
           
        # ToOs handled differently
        dir_root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
        subdirs = ['sv3/ToO']  
        # Fill in the table
        for sd in subdirs:
            for root, dirs, files in os.walk(os.path.join(dir_root,sd)):
                for file in files:
                    if file.endswith('.ecsv'):
                        filename = os.path.join(root,file)
                        data = ascii.read(filename, include_names=["TARGETID","RA","DEC"])
                        cur = self.con.cursor()
                        for datum in zip(data['TARGETID'],data['RA'],data['DEC']):
                            converted_t = ["'{}'".format(str(element)) for element in datum]
                            command = "INSERT INTO {} VALUES (NULL, 'ToO', {});".format(table_name,",".join(converted_t))
                            try:
                                cur.execute(command)
                            except sqlite3.IntegrityError:
                                pass
                        cur.close()
                        self.con.commit()      
            
        # Fill in the table
        dir_root = "/global/cfs/cdirs/desi/target/secondary/sv3/outdata/0.57.0/"
        lunations = ["bright", "dark"]
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
        self.con.commit()
                            
    def load_tiles(self, createTable=False, latestOnly=False):

        table_name = "tiles"
        if createTable:
            command = f"CREATE TABLE tiles (TARGETID INTEGER, TILEID INTEGER, UNIQUE (TARGETID, TILEID));"
            self.con.execute(command)      
            
        # SV3 onward
        dir_root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
        subdirs = ['sv3','sv3/secondary']
        lunations = ["bright", "dark"]    

#         # Fill in the table
        for sd in subdirs:
            for lunation in lunations:
                for root, dirs, files in os.walk(os.path.join(dir_root,sd,lunation)):
                    for file in files:
                        if file.endswith('.ecsv'):
                            filename = os.path.join(root,file)
                            data = ascii.read(filename, include_names=["TARGETID","ZTILEID"])
                            cur = self.con.cursor()
                            for datum in data:
                                try:
                                    cur.execute(f"insert into tiles values ({datum[0]}, {datum[1]})")
                                except sqlite3.IntegrityError:
                                    pass
                            cur.close()
                            self.con.commit()

        if not latestOnly:
#             # SV 1
            filename = "/global/cfs/cdirs/desi/spectro/redux/daily/tiles.csv"
            fbase = "/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk"
            data = ascii.read(filename, include_names=["TILEID","SURVEY"])  
            for d in data:
                if d[1]=='sv1': # and d[0] == 80617:
                    ptile = d[0].astype('str')
                    ptile = ptile.zfill(6)
                    ffile = os.path.join(fbase,ptile[:3],f'fiberassign-{ptile}.fits.gz')
                    fits = fitsio.FITS(ffile)
                    col = fits['FIBERASSIGN']['TARGETID']
                    cur = self.con.cursor()
                    for datum in col.read():
                        try:
                            cur.execute(f"insert into tiles values ({datum}, {d[0]})")
                        except sqlite3.IntegrityError:
                            pass    
                    cur.close()
                    self.con.commit()
                    
#             # SV2 onward
            dir_root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
            subdirs = ['sv2']
            lunations = ["bright", "dark"]    

            # Fill in the table
            for sd in subdirs:
                for lunation in lunations:
                    for root, dirs, files in os.walk(os.path.join(dir_root,sd,lunation)):
                        for file in files:
                            if file.endswith('.ecsv'):
                                filename = os.path.join(root,file)
                                data = ascii.read(filename, include_names=["TARGETID","ZTILEID"])
                                cur = self.con.cursor()
                                for datum in data:
                                    try:
                                        cur.execute(f"insert into tiles values ({datum[0]}, {datum[1]})")
                                    except sqlite3.IntegrityError:
                                        pass
                                cur.close()
                                self.con.commit()
                            

    def load_daily(self,createTable=False):
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
        
        table_name = "daily"
        
        if createTable:
            command = f"CREATE TABLE daily (TILEID INTEGER, YYYYMMDD INTEGER, PANELID INTEGER, TARGETID INTEGER, UNIQUE (TILEID, YYYYMMDD, PANELID, TARGETID));"
            self.con.execute(command)   
        
        #find the last date
        with self.con:
            ans=self.con.execute("SELECT MAX(YYYYMMDD) FROM daily")
        maxdate = ans.fetchone()[0]
        if maxdate is None: maxdate = 0
        print(maxdate)

        dates = []
        for path in glob.glob(f'{dir_root}/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
        dates = numpy.unique(dates)

        for date in dates:
            if int(date) >= maxdate:
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        for i in range(10):
                            fn = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-{date}.fits'
                            try:
                                tids=fitsio.read(fn, "FIBERMAP", columns="TARGETID")
                            except:
                                print(f"{fn} not found")
                                break
                                
                            cur = self.con.cursor()
                            for datum in tids:
                                try:
                                    cur.execute(f"insert into daily values ({tile}, {date}, {i}, {datum})")
                                except sqlite3.IntegrityError:
                                    pass 
                            cur.close()
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
        WHERE tiles.TARGETID = {targetid};
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
    # running as a cron job on cori10
    TileDate.update()