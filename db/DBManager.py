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
    
    
Secondary targets:

    Table name: secondary

    Table sources:
    SV3
    /global/cfs/cdirs/desi/target/secondary/sv3/outdata/0.57.0/
    
    MAIN
    /global/cfs/cdirs/desi/target/secondary/main/outdata/1.0.0/


redshifts-prod:

    Table name: redshifts_prod
    
    Table sources:
    /global/project/projectdirs/desi/spectro/redux/*/zcatalog-*-cumulative.fits
    /global/project/projectdirs/desi/spectro/redux/*/zcatalog-*-pernight.fits
    
zbest_prod:

    Table name: zbest_<prod>
    
    (hdu=1 is zbest, don't save hdu=2 for fibermap)
    Table sources:
    /global/project/projectdirs/desi/spectro/redux/{prod}/tiles/cumulative/*/zbest*.fits 
    /global/project/projectdirs/desi/spectro/redux/{prod}/tiles/pernight/*/zbest*.fits
    
zbest_daily:

    Table name: zbest_daily
    (hdu=1 is zbest, don't save hdu=2 for fibermap)
    Table sources:
    /global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative
    
spectra_prod:

    Table name: spectra_<prod>

    Table sources:
    /global/project/projectdirs/desi/spectro/redux/{prod}/tiles/cumulative/*/spectra*.fits
    
"""

class DBManager:
    
    filename = "/global/cfs/cdirs/desi/science/td/db/desi.db"

    # ONLY DO ToO since not obvious others are needed.
    @staticmethod
    def load_mtl():

        root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
        runs = ["main","sv3","sv2"]
        programs = ["","ToO","secondary"]
        lunations = ["bright","dark",""]

        # The files have column names with SV?_ so purge them

        dfs=[]
        for run in runs:
            dfall_2 = None
            for program in programs:
                for lunation in lunations:
                    print(run, program, lunation)
                    path = os.path.join(root,run,program,lunation)
                    if os.path.isdir(path):
                        for file in glob.glob(path+"/*.ecsv"):
#                             print(file)
                            if "input" in file:
                                continue
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
        dfs.to_sql('mtl',con,if_exists='append')  
        con.close() 

    @staticmethod
    def load_secondary(runs=['sv3','main']): # createTable=False):

        # SV3
        if 'sv3' in runs:
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
                            df['DESITARGET_v']= numpy.full(df.shape[0],'0.57.0')
                            df['RUN']=numpy.full(df.shape[0],'sv3')
                            df['PROGRAM']=numpy.full(df.shape[0],program)
                            df['LUNATION']=numpy.full(df.shape[0],lunation)
                            df.to_sql('secondary',con,if_exists='append')

        # main
        
        # SHEMA EVOLUTION
        # New keywords DESI_TARGET   SCND_TARGET
        fix = '''
        ALTER TABLE secondary
        ADD COLUMN DESI_TARGET INT;
        ALTER TABLE secondary
        ADD COLUMN SCND_TARGET INT;
        '''
        if 'main' in runs:
            dir_root = "/global/cfs/cdirs/desi/target/secondary/main/outdata/1.0.0/"
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
                            df['DESITARGET_v']= numpy.full(df.shape[0],'1.0.0')
                            df['RUN']=numpy.full(df.shape[0],'main')
                            df['PROGRAM']=numpy.full(df.shape[0],program)
                            df['LUNATION']=numpy.full(df.shape[0],lunation)
                            df.to_sql('secondary',con,if_exists='append')
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
 

    # not debugged somehow got lost so rewrote
    @staticmethod
    def load_proposals_pv():       
        runs = ["main","sv3","sv1"]
        lunations = ["BRIGHT","DARK"]
        priorities = ["LOW", "MEDIUM", "HIGH"]
        for run in runs:
            for lunation in lunations:
                for priority in priorities:
                    filename = f"/global/cfs/cdirs/desi/target/proposals/proposals_{run}_year1_frozen/indata/PV_{lunation}_{priority}.fits"
                    try:
                        dat = Table.read(filename, format='fits')
                    except:
                        print(f"{filename} not found")
                        continue                        
                    
                    df = dat.to_pandas()
                    df['PRIORITY']=numpy.full(df.shape[0],priority)
                    df['LUNATION']=numpy.full(df.shape[0],lunation)
                    df.to_sql('proposals_pv',con,if_exists='fail')            
            con.close()
            
    @staticmethod
    def load_zbest_prod(prod="denali"):

        subdirs = ['cumulative','pernight']
        substrs = ['thru','']
        con = sqlite3.connect(DBManager.filename)        
        for subdir,substr in zip(subdirs,substrs):
            # per night
            dir_root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/tiles/{subdir}/"
            dates = []
            for path in glob.glob(f'{dir_root}/*/202?????'):
                split = path.split('/')
                dates.append(split[-1])
            dates = numpy.unique(dates)

            for date in dates:
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        for i in range(10):                            
                            filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-{substr}{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits',hdu=1)
                            except:
                                print(f"{filename} not found")
                                continue

                            for icoeff in range(0,10):
                                dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
                            dat.remove_column('COEFF')
                            df = dat.to_pandas()
                            df['GROUPING'] = numpy.full(df.shape[0],subdir)
                            df['PRODUCTION']=numpy.full(df.shape[0],prod)
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            df.to_sql(f'zbest_{prod}',con,if_exists='append')
        con.close

    @staticmethod
    def load_zbest_daily():
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
        
        con = sqlite3.connect(DBManager.filename)        
        #find the last date
        ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM zbest_daily")
        maxdate = ans.fetchone()[0]
        if maxdate is None: maxdate = 0
        print('max data ',maxdate)
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
                        print(date,tile)
                        for i in range(10):                            
                            
                            # check to see if this file has been done already and if so break
                            if (int(date) == maxdate):
                                ans=con.execute(f"SELECT count(*) FROM zbest_daily WHERE YYYYMMDD={date} AND TILE={tile} AND PETAL={i}")
                                if ans.fetchone()[0] !=0:
                                    continue
                            
                            filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-thru{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits',hdu=1)
                            except:
                                print(f"{filename} not found")
                                continue
                            for icoeff in range(0,10):
                                dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
                            dat.remove_column('COEFF')
                            df = dat.to_pandas()
                            df['GROUPING'] = numpy.full(df.shape[0],'cumulative')
                            df['PRODUCTION']=numpy.full(df.shape[0],'daily')
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            df.to_sql(f'zbest_daily',con,if_exists='append')
        con.close

    @staticmethod
    def load_fibermap_daily(schema=False):
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
        
        con = sqlite3.connect(DBManager.filename)        
        #find the last date
        ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM fibermap_daily")
        maxdate = ans.fetchone()[0]
        if maxdate is None: maxdate = 0

        print('max data ',maxdate)
        dates = []
        for path in glob.glob(f'{dir_root}/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
        dates = numpy.unique(dates)
        
        dfs=[]   # used only for schema
        for date in dates:
            if int(date) >= maxdate:
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        print(date,tile)
                        if schema:
                            i=0
                            filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-thru{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits',hdu=2)
                            except:
                                print(f"{filename} not found")
                                continue
                                
                            df = dat.to_pandas()
                            df['GROUPING'] = numpy.full(df.shape[0],'cumulative')
                            df['PRODUCTION']=numpy.full(df.shape[0],'daily')
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            dfs.append(df[0:0])
                        else:                                
                            for i in range(10):      
                                # check to see if this file has been done already and if so break
                                if (int(date) == maxdate):
                                    ans=con.execute(f"SELECT count(*) FROM fibermap_daily WHERE YYYYMMDD={date} AND TILE={tile} AND PETAL={i}")
                                    if ans.fetchone()[0] !=0:
                                        continue

                                filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-thru{date}.fits'
                                try:
                                    dat = Table.read(filename, format='fits',hdu=2)
                                except:
                                    print(f"{filename} not found")
                                    continue
                                df = dat.to_pandas()
                                df['GROUPING'] = numpy.full(df.shape[0],'cumulative')
                                df['PRODUCTION']=numpy.full(df.shape[0],'daily')
                                df['TILE']=numpy.full(df.shape[0],int(tile))
                                df['PETAL']=numpy.full(df.shape[0],i)
                                df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                                df.to_sql(f'fibermap_daily',con,if_exists='append')
                            
        if schema:
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)
            dfs=dfs[0:0]
            dfs.to_sql(f'fibermap_daily',con,if_exists='replace')
            
        con.close
    
    @staticmethod
    def load_exposure_tables_daily():
        
        dir_root='/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/'
#         maxdate='20201214'
        con = sqlite3.connect(DBManager.filename)
                #find the last date
        ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM exposure_tables_daily")
        maxdate = ans.fetchone()[0]
        if maxdate is None: maxdate = 0
        print('maxdate ',maxdate)
        
        dates = []
        for path in glob.glob(f'{dir_root}/*/exposure_table_????????.csv'):
            split = path.split('/')
            date = split[-1][-12:-4]
            dates.append(date)

        dates=numpy.sort(dates)

        dfs=[]
        for date in dates:
            if int(date) <= maxdate:
                continue
            file = f'/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/{date[0:6]}/exposure_table_{date}.csv'
            data = ascii.read(file)
            df = data.to_pandas()
            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
            dfs.append(df)
          
        if len(dfs)>0:
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)
            dfs.to_sql(f'exposure_tables_daily',con,if_exists='append')
            
    @staticmethod
    def load_spectra_prod(prod="denali"):
        # from what I can tell all the information in daily is in cumulative
        subdirs = ['cumulative']
        substrs = ['thru']
        con = sqlite3.connect(DBManager.filename)        
        for subdir,substr in zip(subdirs,substrs):
            # per night
            dir_root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/tiles/{subdir}/"
            dates = []
            for path in glob.glob(f'{dir_root}/*/202?????'):
                split = path.split('/')
                dates.append(split[-1])
            dates = numpy.unique(dates)

            dfs=[]
            for date in dates:
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        for i in range(10):                            
                            filename = f'{dir_root}/{tile}/{date}/spectra-{i}-{tile}-{substr}{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits')
                            except:
                                print(f"{filename} not found")
                                continue
                            df = dat.to_pandas()
                            df['GROUPING'] = numpy.full(df.shape[0],subdir)
                            df['PRODUCTION']=numpy.full(df.shape[0],prod)
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            dfs.append(df)
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)                        
            dfs.to_sql(f'spectra_{prod}',con,if_exists='append')
        con.close
        
#     @staticmethod
#     # deprecated
#     def load_spectra(prod="daily"):
#         if prod == "daily":
#             dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
#         else:
#             dir_root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/tiles/pernight/"
        
#         con = sqlite3.connect(DBManager.filename)        
#         #find the last date
#         ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM spectra WHERE PRODUCTION = '{prod}'")
#         maxdate = ans.fetchone()[0]
#         if maxdate is None: maxdate = 0
#         print(maxdate)

#         dates = []
#         for path in glob.glob(f'{dir_root}/*/202?????'):
#             split = path.split('/')
#             dates.append(split[-1])
#         dates = numpy.unique(dates)

#         dfs=[]
#         for date in dates:
#             if int(date) >= maxdate:
#                 for path in glob.glob(f'{dir_root}/*/{date}'):
#                     split = path.split('/')
#                     tile = split[-2]
#                     if tile.isnumeric():
#                         for i in range(10):                            
                            
#                             # check to see if this file has been done already and if so break
#                             if (int(date) == maxdate):
#                                 ans=con.execute(f"SELECT count(*) FROM spectra WHERE PRODUCTION = '{prod}' AND YYYYMMDD={date} AND TILE={tile} AND PETAL={i}")
#                                 if ans.fetchone()[0] !=0:
#                                     print(f"skipping {date} {tile} {i}")
#                                     break
                            
#                             filename = f'{dir_root}/{tile}/{date}/spectra-{i}-{tile}-{date}.fits'
#                             print(filename)
#                             try:
#                                 dat = Table.read(filename, format='fits')
#                             except:
#                                 print(f"{filename} not found")
#                                 break
                                
#                             df = dat.to_pandas()
#                             df['PRODUCTION']=numpy.full(df.shape[0],prod)
#                             df['TILE']=numpy.full(df.shape[0],int(tile))
#                             df['PETAL']=numpy.full(df.shape[0],i)
#                             df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
#                             dfs.append(df)
#         if len(dfs) != 0:
#             dfs = pandas.concat(dfs, ignore_index=True, sort=False)                    
#             dfs.to_sql('spectra',con,if_exists='append')  
#         con.close    

#     deprecated
#     def load_mtl():

#         root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
#         runs = ["sv3","sv2"]
#         programs = ["","ToO","secondary"]
#         lunations = ["bright","dark"]

#         # The files have column names with SV?_ so purge them

#         dfs=[]
#         for run in runs:
#             dfall_2 = None
#             for program in programs:
#                 for lunation in lunations:
#                     path = os.path.join(root,run,program,lunation)
#                     if os.path.isdir(path):
#                         for file in glob.glob(path+"/*.ecsv"):
#                             print(file)
#                             data = ascii.read(file)
#                             df = data.to_pandas()

#                             df['RUN']=numpy.full(df.shape[0],run)
#                             df['PROGRAM']=numpy.full(df.shape[0],program)
#                             df['LUNATION']=numpy.full(df.shape[0],lunation)
                            
#                             dfs.append(df)
#                     else:
#                         print (path, 'not exists')

#         dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                    
#         con = sqlite3.connect(DBManager.filename)
#         dfs.to_sql('mtl',con,if_exists='replace')  
#         con.close()
        
    @staticmethod
    def load_daily():
        DBManager.load_zbest_daily()
        DBManager.load_fibermap_daily()
        DBManager.load_exposure_tables_daily()
#         DBManager.load_spectra(prod="daily")

    @staticmethod
    def byProgram(program):

        command = f'''SELECT DISTINCT secondary.PROGRAM, secondary.TARGETID, secondary.RA, secondary.DEC, zbest.PRODUCTION, zbest.YYYYMMDD, zbest.Z, zbest.ZERR, zbest.SPECTYPE, redshifts_prod.PRODUCTION, redshifts_prod.Z, redshifts_prod.ZERR, redshifts_prod.SPECTYPE, redshifts_prod.DELTACHI2, redshifts_prod.ZWARN
        FROM secondary
            LEFT JOIN zbest using(TARGETID)
            LEFT JOIN redshifts_prod using(TARGETID)
            WHERE  secondary.PROGRAM LIKE "PV%" and redshifts_prod.COADD="cumulative"
        ORDER BY
            secondary.PROGRAM,
            secondary.TARGETID
        LIMIT 10;
        ''' 

        
if __name__ == "__main__":
    # running as a cron job on cori10
    DBManager.load_daily()
