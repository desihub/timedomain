import os
import sys
import glob
from astropy.io import ascii
import pandas
import psycopg2
import re
from astropy.table import Table
import sqlite3

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
def dtypesToSchema(dtypes):
    for index,value in dtypes.items():
        if numpy.issubdtype(value, numpy.integer):
            nvalue ='INTEGER'
        elif numpy.issubdtype(value, numpy.float):
            nvalue ='REAL'
        else:
            nvalue = 'TEXT'
        print(f'"{index}"  {nvalue},')
        

                                    
class mtl():

    schema="""
        CREATE TABLE IF NOT EXISTS "mtl" (
          "RA" REAL,
          "DEC" REAL,
          "PMRA" REAL,
          "PMDEC" REAL,
          "REF_EPOCH" REAL,
          "FLUX_G" REAL,
          "FLUX_R" REAL,
          "FLUX_Z" REAL,
          "PARALLAX" REAL,
          "GAIA_PHOT_G_MEAN_MAG" REAL,
          "GAIA_PHOT_BP_MEAN_MAG" REAL,
          "GAIA_PHOT_RP_MEAN_MAG" REAL,
          "GAIA_ASTROMETRIC_EXCESS_NOISE" REAL,
          "TARGETID" INTEGER,
          "DESI_TARGET" INTEGER,
          "BGS_TARGET" INTEGER,
          "MWS_TARGET" INTEGER,
          "SCND_TARGET" INTEGER,
          "SV3_DESI_TARGET" INTEGER,
          "SV3_SCND_TARGET" INTEGER,
          "SV3_BGS_TARGET"       INTEGER,
          "SV3_MWS_TARGET"       INTEGER,
          "SCND_ORDER" INTEGER,
          "NUMOBS_MORE"  INTEGER,
          "Z"            REAL,
          "ZWARN"        INTEGER,
          "ZTILEID"      INTEGER,
          "Z_QN"          REAL,
          "IS_QSO_QN"       INTEGER,
          "DELTACHI2"        REAL,
          "PRIORITY_INIT" INTEGER,
          "SUBPRIORITY" REAL,
          "NUMOBS_INIT" INTEGER,
          "NUMOBS" INTEGER,
          "OBSCONDITIONS" INTEGER,
          "CHECKER" TEXT,
          "TOO_TYPE" TEXT,
          "TOO_PRIO" TEXT,
          "OCLAYER" TEXT,
          "MJD_BEGIN" REAL,
          "MJD_END" REAL,
          "TOOID" INTEGER,
          "TARGET_STATE"      TEXT,
          "TIMESTAMP"         TEXT,
          "VERSION"           TEXT,
          "PRIORITY"           INTEGER,         
          "RUN" TEXT,
          "PROGRAM" TEXT,
          "LUNATION" TEXT
        );
        """

    root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
    runs = ["main","sv3","sv2"]
    runs=["sv3"]
    programs = ["","ToO","secondary"]
    lunations = ["bright","dark",""]
#     conn = psycopg2.connect("dbname=test user=postgres")
    conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/temp.db")    

    @staticmethod
    def create_table(overwrite = False):
        cur = mtl.conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS mtl;")
        cur.execute(mtl.schema)
        
        cur.close()
        
    @staticmethod    
    def fill_table():

        # The files have column names with SV?_ so purge them
        for run in mtl.runs:
            for program in mtl.programs:
                for lunation in mtl.lunations:
                    path = os.path.join(mtl.root,run,program,lunation)
                    if os.path.isdir(path):
                        print(run, program, lunation)
                        for file in glob.glob(path+"/*.ecsv"):
                            if "input" in file:
                                continue
                            data = ascii.read(file)
                            df = data.to_pandas()
                            df['RUN']=numpy.full(df.shape[0],run)
                            df['PROGRAM']=numpy.full(df.shape[0],program)
                            df['LUNATION']=numpy.full(df.shape[0],lunation)
                            try:
                                df.to_sql('mtl',mtl.conn,index=False,if_exists='append')
                            except sqlite3.OperationalError as err:
                                dtypesToSchema(df.dtypes)

                                print(df.dtypes)
                                sys.exit()
                    else:
                        print (path, 'not exists')

class secondary:
    
    schema="""
    CREATE TABLE IF NOT EXISTS "secondary" (
        "RA"  REAL,
        "DEC"  REAL,
        "PMRA"  TEXT,
        "PMDEC"  TEXT,
        "REF_EPOCH"  TEXT,
        "OVERRIDE"  TEXT,
        "FLUX_G"  TEXT,
        "FLUX_R"  TEXT,
        "FLUX_Z"  TEXT,
        "PARALLAX"  TEXT,
        "GAIA_PHOT_G_MEAN_MAG"  TEXT,
        "GAIA_PHOT_BP_MEAN_MAG"  TEXT,
        "GAIA_PHOT_RP_MEAN_MAG"  TEXT,
        "GAIA_ASTROMETRIC_EXCESS_NOISE"  TEXT,
        "TARGETID"  INTEGER,
        "SV3_DESI_TARGET"  INTEGER,
        "SV3_SCND_TARGET"  INTEGER,
        "DESI_TARGET" INTEGER,
        "SCND_TARGET" INTEGER,
        "SCND_ORDER"  INTEGER,
        "DESITARGET_v"  TEXT,
        "RUN"  TEXT,
        "PROGRAM"  TEXT,
        "LUNATION"  TEXT
    );
    """
    
    runs=['sv3','main']
    runs=['main']
    conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/temp.db")  
    
    @staticmethod
    def create_table(overwrite = False):
        cur = secondary.conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS secondary;")
        cur.execute(secondary.schema)
        cur.close()
        
    
    @staticmethod
    def fill_table():
        
        # SV3
        if 'sv3' in secondary.runs:
            dir_root = "/global/cfs/cdirs/desi/target/secondary/sv3/outdata/0.57.0/"
            lunations = ["bright", "dark"]
       
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
                            try:
                                df.to_sql('secondary',secondary.conn,index=False,if_exists='append')
                            except sqlite3.OperationalError as err:
                                dtypesToSchema(df.dtypes)
                                sys.exit()
                            

        # main
        
        # SCHEMA EVOLUTION
        # New keywords DESI_TARGET   SCND_TARGET
        fix = '''
        ALTER TABLE secondary
        ADD COLUMN DESI_TARGET INT;
        ALTER TABLE secondary
        ADD COLUMN SCND_TARGET INT;
        '''
        if 'main' in secondary.runs:
            dir_root = "/global/cfs/cdirs/desi/target/secondary/main/outdata/1.0.0/"
            lunations = ["bright", "dark"]
       
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
                            try:
                                df.to_sql('secondary',secondary.conn,index=False,if_exists='append')
                            except sqlite3.OperationalError as err:
                                dtypesToSchema(df.dtypes)
                                sys.exit()
                                
class zcatalog_prod:
    
    schema="""
    CREATE TABLE IF NOT EXISTS "zcatalog_prod" (
        "TARGETID"  INTEGER,
        "CHI2"  REAL,
        "Z"  REAL,
        "ZERR"  REAL,
        "ZWARN"  INTEGER,
        "NPIXELS"  INTEGER,
        "SPECTYPE"  TEXT,
        "SUBTYPE"  TEXT,
        "NCOEFF"  INTEGER,
        "DELTACHI2"  REAL,
        "NUMEXP"  INTEGER,
        "NUMTILE"  INTEGER,
        "PETAL_LOC"  INTEGER,
        "DEVICE_LOC"  INTEGER,
        "LOCATION"  INTEGER,
        "FIBER"  INTEGER,
        "FIBERSTATUS"  INTEGER,
        "TARGET_RA"  REAL,
        "TARGET_DEC"  REAL,
        "PMRA"  TEXT,
        "PMDEC"  TEXT,
        "REF_EPOCH"  TEXT,
        "LAMBDA_REF"  TEXT,
        "FA_TARGET"  INTEGER,
        "FA_TYPE"  INTEGER,
        "OBJTYPE"  TEXT,
        "PRIORITY"  INTEGER,
        "SUBPRIORITY"  REAL,
        "OBSCONDITIONS"  INTEGER,
        "RELEASE"  INTEGER,
        "BRICKID"  INTEGER,
        "BRICK_OBJID"  INTEGER,
        "MORPHTYPE"  TEXT,
        "FLUX_G"  TEXT,
        "FLUX_R"  TEXT,
        "FLUX_Z"  TEXT,
        "FLUX_IVAR_G"  TEXT,
        "FLUX_IVAR_R"  TEXT,
        "FLUX_IVAR_Z"  TEXT,
        "MASKBITS"  INTEGER,
        "REF_ID"  INTEGER,
        "REF_CAT"  TEXT,
        "GAIA_PHOT_G_MEAN_MAG"  TEXT,
        "GAIA_PHOT_BP_MEAN_MAG"  TEXT,
        "GAIA_PHOT_RP_MEAN_MAG"  TEXT,
        "PARALLAX"  TEXT,
        "BRICKNAME"  TEXT,
        "EBV"  TEXT,
        "FLUX_W1"  TEXT,
        "FLUX_W2"  TEXT,
        "FIBERFLUX_G"  TEXT,
        "FIBERFLUX_R"  TEXT,
        "FIBERFLUX_Z"  TEXT,
        "FIBERTOTFLUX_G"  TEXT,
        "FIBERTOTFLUX_R"  TEXT,
        "FIBERTOTFLUX_Z"  TEXT,
        "SERSIC"  TEXT,
        "SHAPE_R"  TEXT,
        "SHAPE_E1"  TEXT,
        "SHAPE_E2"  TEXT,
        "PHOTSYS"  TEXT,
        "PRIORITY_INIT"  INTEGER,
        "NUMOBS_INIT"  INTEGER,
        "SV2_DESI_TARGET"  INTEGER,
        "SV2_BGS_TARGET"  INTEGER,
        "SV2_MWS_TARGET"  INTEGER,
        "SV2_SCND_TARGET"  INTEGER,
        "DESI_TARGET"  INTEGER,
        "BGS_TARGET"  INTEGER,
        "MWS_TARGET"  INTEGER,
        "TILEID"  INTEGER,
        "COADD_NUMEXP"  INTEGER,
        "COADD_EXPTIME"  TEXT,
        "MEAN_DELTA_X"  TEXT,
        "RMS_DELTA_X"  TEXT,
        "MEAN_DELTA_Y"  TEXT,
        "RMS_DELTA_Y"  TEXT,
        "MEAN_FIBER_X"  TEXT,
        "MEAN_FIBER_Y"  TEXT,
        "MEAN_FIBER_RA"  REAL,
        "MEAN_FIBER_DEC"  REAL,
        "MEAN_FIBERASSIGN_X"  TEXT,
        "MEAN_FIBERASSIGN_Y"  TEXT,
        "FIRST_NIGHT"  INTEGER,
        "LAST_NIGHT"  INTEGER,
        "NUM_NIGHT"  INTEGER,
        "FIRST_EXPID"  INTEGER,
        "LAST_EXPID"  INTEGER,
        "NUM_EXPID"  INTEGER,
        "FIRST_TILEID"  INTEGER,
        "LAST_TILEID"  INTEGER,
        "NUM_TILEID"  INTEGER,
        "FIRST_FIBER"  INTEGER,
        "LAST_FIBER"  INTEGER,
        "NUM_FIBER"  INTEGER,
        "FIRST_MJD"  TEXT,
        "LAST_MJD"  TEXT,
        "NUM_MJD"  INTEGER,
        "SV1_SCND_TARGET"  INTEGER,
        "SV1_BGS_TARGET"  INTEGER,
        "SV1_DESI_TARGET"  INTEGER,
        "SV1_MWS_TARGET"  INTEGER,
        "CMX_TARGET"  INTEGER,
        "COEFF_0"  REAL,
        "COEFF_1"  REAL,
        "COEFF_2"  REAL,
        "COEFF_3"  REAL,
        "COEFF_4"  REAL,
        "COEFF_5"  REAL,
        "COEFF_6"  REAL,
        "COEFF_7"  REAL,
        "COEFF_8"  REAL,
        "COEFF_9"  REAL,
        "PRODUCTION"  TEXT,
        "COADD"  TEXT
    );
    """

    conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/temp.db")     
    
    @staticmethod
    def create_table(overwrite = False):
        cur = secondary.conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS zcatalog_prod;")
        cur.execute(zcatalog_prod.schema)
        cur.close()

    @staticmethod
    def fill_table(prod='denali'):
        root = f"/global/project/projectdirs/desi/spectro/redux/{prod}"
        coadds = ["cumulative","pernight"]
        coadds = ["cumulative"]
        
        for coadd in coadds:
            filename = os.path.join(root,f"zcatalog-denali-{coadd}.fits")
            print(filename)
            dat = Table.read(filename, format='fits')
            
            # There is a multidimensional column that needs to be broken up
            for icoeff in range(0,10):
                dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
            dat.remove_column('COEFF')
            dat.convert_bytestring_to_unicode()
            
            df = dat.to_pandas()
            df['PRODUCTION']=numpy.full(df.shape[0],prod)
            df['COADD']=numpy.full(df.shape[0],coadd)
            try:
                df.to_sql('zcatalog_prod',zcatalog_prod.conn,index=False,if_exists='append') 
            except sqlite3.OperationalError as err:
                dtypesToSchema(df.dtypes)
                sys.exit()     

                
class dr9_pv:
    schema="""
        CREATE TABLE IF NOT EXISTS "dr9_pv" (
        "OBJID"  INTEGER,
        "BRICKID"  INTEGER,
        "BRICKNAME"  TEXT,
        "RA"  REAL,
        "DEC"  REAL,
        "TYPE"  TEXT,
        "SERSIC"  TEXT,
        "Z_PHOT_MEDIAN"  TEXT,
        "Z_PHOT_L95"  TEXT,
        "mag_g"  TEXT,
        "mag_r"  TEXT,
        "mag_z"  TEXT,
        "mag_B"  TEXT,
        "mag_g_err"  TEXT,
        "mag_r_err"  TEXT,
        "mag_z_err"  TEXT,
        "fibre_mag_g"  TEXT,
        "fibre_mag_r"  TEXT,
        "fibre_mag_z"  TEXT,
        "uncor_radius"  TEXT,
        "BA_ratio"  TEXT,
        "circ_radius"  TEXT,
        "pos_angle"  TEXT,
        "inSGA"  TEXT,
        "inBGS"  TEXT,
        "inlocalbright"  TEXT,
        "inspecfootprint"  TEXT,
        "SGA_pa"  TEXT,
        "SGA_ba"  TEXT,
        "SB_D25_g"  TEXT,
        "SB_D25_r"  TEXT,
        "SB_D25_z"  TEXT,
        "RADIUS_SB25"  TEXT,
        "SGA_MORPHTYPE"  TEXT,
        "SGA_ID"  INTEGER,
        "SGA_redshift"  TEXT,
        "size_SGA"  TEXT,
        "PMRA"  TEXT,
        "PMDEC"  TEXT,
        "REF_EPOCH"  TEXT,
        "OVERRIDE"  TEXT,
        "PVTYPE"  TEXT,
        "PVPRIORITY"  INTEGER,
        "POINTINGID"  INTEGER,
        "TARGET"  TEXT,
        "SURVEY"  TEXT);
    """

    conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/temp.db")     
    
    @staticmethod
    def create_table(overwrite = False):
        cur = secondary.conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS dr9_pv;")
        cur.execute(dr9_pv.schema)
        cur.close()

    @staticmethod
    def fill_table(prod='denali'):      
        dirs=["savepath_dr9","savepath_dr9_corr"]
        surveys=["sv3","main"]
        targets = ["ext","fp","sga","tf"]
        for di, survey in zip(dirs,surveys):
            for target in targets:
                filename = f"/global/homes/k/ksaid/desi_pv/{di}/pv_{target}_full.fits"
                try:
                    dat = Table.read(filename, format='fits')
                except:
                    print(f"{filename} not found")
                    continue
                dat.convert_bytestring_to_unicode()
                df = dat.to_pandas()
                df['TARGET']=numpy.full(df.shape[0],target)
                df['SURVEY']=numpy.full(df.shape[0],survey)
                
                try:
                    df.to_sql('dr9_pv',dr9_pv.conn,index=False,if_exists='append')
                except sqlite3.OperationalError as err:
                    dtypesToSchema(df.dtypes)
                    sys.exit()                 
                

class exposure_tables_daily:
    schema="""
        CREATE TABLE IF NOT EXISTS "exposure_tables_daily" (
            "EXPID"  INTEGER,
            "EXPTIME"  REAL,
            "OBSTYPE"  TEXT,
            "SPECTROGRAPHS"  INTEGER,
            "CAMWORD"  TEXT,
            "TILEID"  INTEGER,
            "NIGHT"  INTEGER,
            "EXPFLAG"  INTEGER,
            "HEADERERR"  TEXT,
            "SURVEY"  INTEGER,
            "SEQNUM"  INTEGER,
            "SEQTOT"  INTEGER,
            "PROGRAM"  TEXT,
            "MJD-OBS"  REAL,
            "REQRA"  REAL,
            "REQDEC"  REAL,
            "TARGTRA"  REAL,
            "TARGTDEC"  REAL,
            "COMMENTS"  TEXT,
            "YYYYMMDD"  INTEGER,
            "PURPOSE"  TEXT,
            "FA_SURV"  TEXT,
            "BADCAMWORD"  TEXT,
            "BADAMPS"  TEXT,
            "LASTSTEP"  TEXT,
            "EFFTIME_ETC"  REAL,
            "FAPRGRM"  TEXT,
            "GOALTIME"  REAL,
            "GOALTYPE"  TEXT,
            "EBVFAC"  REAL,
            "AIRFAC"  REAL,
            "SPEED"  REAL
    );
    """

    conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/temp.db")     
    
    @staticmethod
    def create_table(overwrite = False):
        cur = exposure_tables_daily.conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS exposure_tables_daily;")
        cur.execute(exposure_tables_daily.schema)
        cur.close()

    @staticmethod
    def fill_table():      
        dir_root='/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/'
#         maxdate='20201214'
        con = sqlite3.connect(DBManager.filename)
                #find the last date
        try:
            ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM exposure_tables_daily")
            maxdate = ans.fetchone()[0]
        except:
            maxdate=0
        
        print(maxdate)
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
            print('saving to db')
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)
            try:
                dfs.to_sql(f'exposure_tables_daily',exposure_tables_daily.conn,index=False,if_exists='append')

            except sqlite3.OperationalError as err:
                dtypesToSchema(dfs.dtypes)
                sys.exit()                 
        con.close()
        
    @staticmethod
    def fill_sga():
        cur = secondary.conn.cursor()
        
        if overwrite:
            cur.execute("UPDATE exposure_tables_daily SET SGA_ID=dr9_pv.SGA_ID FROM dr9 WHERE exposure_tables_daily.OBJID=dr9_pv.OBJID AND exposure_tables_daily.BRICKID=dr9_pv.BRICKID;")
        cur.execute(exposure_tables_daily.schema)
        cur.close()    

        
class DBManager:
    
    filename = "/global/cfs/cdirs/desi/science/td/db/desi.db"


 

# 
    # the first two nights of observation do not have anything in cumulative
    @staticmethod
    def load_zbest_daily():
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
        maxdate=None
        con = sqlite3.connect(DBManager.filename)        
        #find the last date
        ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM zbest_daily")
        maxdate = ans.fetchone()[0]
        if maxdate is None: 
            maxdate = 20201214
        print('max date ',maxdate)
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
                            dat.convert_bytestring_to_unicode()
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
        maxdate=None
        con = sqlite3.connect(DBManager.filename)        
        #find the last date
#         ans=con.execute(f"SELECT MAX(YYYYMMDD) FROM fibermap_daily")
#         maxdate = ans.fetchone()[0]
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
                            dat.convert_bytestring_to_unicode()    
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

        
#    @staticmethod
#     def load_zbest_prod(prod="denali"):

#         subdirs = ['cumulative','pernight']
#         substrs = ['thru','']
#         con = sqlite3.connect(DBManager.filename)        
#         for subdir,substr in zip(subdirs,substrs):
#             # per night
#             dir_root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/tiles/{subdir}/"
#             dates = []
#             for path in glob.glob(f'{dir_root}/*/202?????'):
#                 split = path.split('/')
#                 dates.append(split[-1])
#             dates = numpy.unique(dates)

#             for date in dates:
#                 for path in glob.glob(f'{dir_root}/*/{date}'):
#                     split = path.split('/')
#                     tile = split[-2]
#                     if tile.isnumeric():
#                         for i in range(10):                            
#                             filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-{substr}{date}.fits'
#                             try:
#                                 dat = Table.read(filename, format='fits',hdu=1)
#                             except:
#                                 print(f"{filename} not found")
#                                 continue

#                             for icoeff in range(0,10):
#                                 dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
#                             dat.remove_column('COEFF')
#                             df = dat.to_pandas()
#                             df['GROUPING'] = numpy.full(df.shape[0],subdir)
#                             df['PRODUCTION']=numpy.full(df.shape[0],prod)
#                             df['TILE']=numpy.full(df.shape[0],int(tile))
#                             df['PETAL']=numpy.full(df.shape[0],i)
#                             df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
#                             df.to_sql(f'zbest_{prod}',con,if_exists='append')
#         con.close        
        
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
