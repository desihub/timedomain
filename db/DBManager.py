import os
import glob
from astropy.io import ascii
import pandas
import sqlite3
import re
from astropy.table import Table


import fitsio

import numpy

"""
Tables to read:

Tiles

    SV2, SV3:
    /global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
    
    
Secondary targets:

    SV3
    /global/cfs/cdirs/desi/target/secondary/sv3/outdata/0.57.0/


redshifts-prod:

    /global/project/projectdirs/desi/spectro/redux/*/zcatalog-denali-cumulative.fits
    /global/project/projectdirs/desi/spectro/redux/*/zcatalog-denali-pernight.fits
    
"""

class DBManager:
    
    filename = "/global/cfs/cdirs/desi/science/td/db/desi.db"

    #####
    #
    #  The files have column names with SV?_ so purge them
    #
    @staticmethod
    def load_mtl():

        root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
        runs = ["sv3","sv2"]
        programs = ["","ToO","secondary"]
        lunations = ["bright","dark"]

        # The files have column names with SV?_ so purge them

        dfs=[]
        for run in runs:
            dfall_2 = None
            for program in programs:
                for lunation in lunations:
                    path = os.path.join(root,run,program,lunation)
                    if os.path.isdir(path):
                        for file in glob.glob(path+"/*.ecsv"):
                            print(file)
                            data = ascii.read(file)
                            df = data.to_pandas()

                            df['RUN']=numpy.full(df.shape[0],run)
                            df['PROGRAM']=numpy.full(df.shape[0],program)
                            df['LUNATION']=numpy.full(df.shape[0],lunation)
                            
                            dfs.append(df)
                    else:
                        print (path, 'not exists')

        dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                    
        con = sqlite3.connect(DBManager.filename)
        dfs.to_sql('mtl',con,if_exists='replace')  
        con.close()

    @staticmethod
    def load_secondary(): # createTable=False):

        # Fill in the table
        dir_root = "/global/cfs/cdirs/desi/target/secondary/sv3/outdata/0.57.0/"
        lunations = ["bright", "dark"]

        con = sqlite3.connect(DBManager.filename)         
        for lunation in lunations:
            for root, dirs, files in os.walk(os.path.join(dir_root,lunation)):
                for file in files:
                    if file.endswith('.fits'):
                        filename = os.path.join(root,file)
                        print(file)
                        program = file[:-5]
                        
                        dat = Table.read(filename, format='fits')
                        df = dat.to_pandas()
                        df['PROGRAM']=numpy.full(df.shape[0],program)
                        df['LUNATION']=numpy.full(df.shape[0],lunation)
                        df.to_sql('secondary',con,if_exists='fail')

        con.close()
        
    @staticmethod
    def load_redshifts_prod(prod='denali'):
        root = f"/global/project/projectdirs/desi/spectro/redux/{prod}"
        coadds = ["cumulative","pernight"]
        
        con = sqlite3.connect(DBManager.filename) 
        for coadd in coadds:
            filename = os.path.join(root,f"zcatalog-denali-{coadd}.fits")
            print(filename)
            dat = Table.read(filename, format='fits')
            
            # There is a multidimensional column that needs to be broken up
            for icoeff in range(0,10):
                dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
            dat.remove_column('COEFF')
            
            df = dat.to_pandas()
            df['PRODUCTION']=numpy.full(df.shape[0],prod)
            df['COADD']=numpy.full(df.shape[0],coadd)
            df.to_sql('redshifts_prod',con,if_exists='fail')            
        con.close()
 
    @staticmethod
    def load_zbest(prod="daily"):
        if prod == "daily":
            dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
        else:
            dir_root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/tiles/pernight/"
        
        con = sqlite3.connect(DBManager.filename)        
        #find the last date
        ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM zbest WHERE PRODUCTION = '{prod}'")
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
                            
                            # check to see if this file has been done already and if so break
                            if (int(date) == maxdate):
                                ans=con.execute(f"SELECT count(*) FROM zbest WHERE PRODUCTION = '{prod}' AND YYYYMMDD={date} AND TILE={tile} AND PETAL={i}")
                                if ans.fetchone()[0] !=0:
                                    break
                            
                            filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-{date}.fits'
                            print(filename)
                            try:
                                dat = Table.read(filename, format='fits')
                            except:
                                print(f"{filename} not found")
                                break
                                
                            for icoeff in range(0,10):
                                dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
                            dat.remove_column('COEFF')
                            df = dat.to_pandas()
                            df['PRODUCTION']=numpy.full(df.shape[0],prod)
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            df.to_sql('zbest',con,if_exists='append')
        con.close

    @staticmethod
    def load_spectra(prod="daily"):
        if prod == "daily":
            dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
        else:
            dir_root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/tiles/pernight/"
        
        con = sqlite3.connect(DBManager.filename)        
        #find the last date
        ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM spectra WHERE PRODUCTION = '{prod}'")
        maxdate = ans.fetchone()[0]
        if maxdate is None: maxdate = 0
        print(maxdate)

        dates = []
        for path in glob.glob(f'{dir_root}/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
        dates = numpy.unique(dates)

        dfs=[]
        for date in dates:
            if int(date) >= maxdate:
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        for i in range(10):                            
                            
                            # check to see if this file has been done already and if so break
                            if (int(date) == maxdate):
                                ans=con.execute(f"SELECT count(*) FROM spectra WHERE PRODUCTION = '{prod}' AND YYYYMMDD={date} AND TILE={tile} AND PETAL={i}")
                                if ans.fetchone()[0] !=0:
                                    print(f"skipping {date} {tile} {i}")
                                    break
                            
                            filename = f'{dir_root}/{tile}/{date}/spectra-{i}-{tile}-{date}.fits'
                            print(filename)
                            try:
                                dat = Table.read(filename, format='fits')
                            except:
                                print(f"{filename} not found")
                                break
                                
                            df = dat.to_pandas()
                            df['PRODUCTION']=numpy.full(df.shape[0],prod)
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            dfs.append(df)
        if len(dfs != 0):
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)                    
            dfs.to_sql('spectra',con,if_exists='append')  
        con.close    

    @staticmethod
    def load_daily():
        DBManager.load_zbest(prod="daily")        
        DBManager.load_spectra(prod="daily")

    @staticmethod
    def byProgram(program):

        command = f'''SELECT DISTINCT secondary.PROGRAM, secondary.TARGETID, secondary.RA, secondary.DEC, zbest.PRODUCTION, zbest.YYYYMMDD, zbest.Z, zbest.ZERR, zbest.SPECTYPE, redshifts_prod.PRODUCTION, redshifts_prod.Z, redshifts_prod.ZERR, redshifts_prod.SPECTYPE, redshifts_prod.DELTACHI2, redshifts_prod.ZWARN
        FROM zbest
        INNER JOIN secondary
            ON secondary.targetid = zbest.targetid
        INNER JOIN redshifts_prod
            ON secondary.targetid = redshifts_prod.targetid
        WHERE secondary.PROGRAM LIKE "PV%" and redshifts_prod.COADD="cumulative"
        ORDER BY
            secondary.PROGRAM,
            secondary.TARGETID
        LIMIT 10;
        '''
        
        command = f'''SELECT DISTINCT secondary.TARGETID, secondary.RA, secondary.DEC, spectra.EXPID, spectra.FIBER_RA, spectra.FIBER_DEC
        FROM secondary
        INNER JOIN spectra
            ON secondary.targetid = spectra.targetid
        WHERE secondary.PROGRAM LIKE "PV%"
        ORDER BY
            secondary.PROGRAM,
            secondary.TARGETID
        LIMIT 10;
        '''
        command = f'''SELECT DISTINCT secondary.TARGETID, secondary.RA, secondary.DEC, spectra.EXPID, spectra.FIBER_RA, spectra.FIBER_DEC
        FROM secondary
        INNER JOIN spectra
            ON secondary.targetid = spectra.targetid
        WHERE secondary.PROGRAM LIKE "PV%"
        ORDER BY
            secondary.PROGRAM,
            secondary.TARGETID
        LIMIT 10;
        '''        
    
if __name__ == "__main__":
    # running as a cron job on cori10
    DBManager.load_daily()