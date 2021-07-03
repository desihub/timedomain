import os
import sys
import glob
from astropy.io import ascii
import pandas
import psycopg2
import re
from astropy.table import Table
import sqlite3
from sqlalchemy import create_engine
from sqlalchemy import text
import sqlalchemy
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

try:
    with open('/global/homes/a/akim/secrets/desi_pg.txt', 'r') as f:
        data = f.read()
        uname, passwd = data.split()
        engine = create_engine(f'postgresql://{uname}:{passwd}@decatdb.lbl.gov/desi')
        
#     conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/desi.db") 

except:
    print("No connection")

def dtypesToSchema(dtypes):
    for index,value in dtypes.items():
        if value == numpy.float64:
            nvalue = 'DOUBLE PRECISION'
        elif value == numpy.float32:
            nvalue = 'REAL'
        elif value == numpy.int64:
            nvalue = 'BIGINT'
        elif value == numpy.int32:
            nvalue = 'INTEGER'
        elif value == numpy.int16 or value == numpy.uint8 :
            nvalue = 'SMALLINT'
        elif value == numpy.object:
            nvalue = 'TEXT'
        else:
            print ("unknown ",value)
#         print(f'"{index}"  {nvalue}, {value}')
        print(f'"{index}"  {nvalue},')       
                                    
class mtl:

    schema=f"""
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

    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()
        
        if overwrite:
            cur.execute(f"DROP TABLE IF EXISTS mtl;")
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
                                df.to_sql('mtl',conn,index=False,if_exists='append')
                            except sqlite3.OperationalError as err:
                                dtypesToSchema(df.dtypes)

                                print(df.dtypes)
                                sys.exit()
                    else:
                        print (path, 'not exists')

class secondary:

    schema=f"""
    CREATE TABLE IF NOT EXISTS secondary (
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
    
    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()
        
        if overwrite:
            cur.execute(f"DROP TABLE IF EXISTS secondary;")
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
                                df.to_sql('secondary',conn,index=False,if_exists='append')
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
                                df.to_sql('secondary',conn,index=False,if_exists='append')
                            except sqlite3.OperationalError as err:
                                dtypesToSchema(df.dtypes)
                                sys.exit()
                                
class zcatalog_denali_cumulative:
    
    schema="""
    CREATE TABLE IF NOT EXISTS "zcatalog_denali_cumulative" (
        "TARGETID"  BIGINT,
        "CHI2"  DOUBLE PRECISION,
        "Z"  DOUBLE PRECISION,
        "ZERR"  DOUBLE PRECISION,
        "ZWARN"  BIGINT,
        "NPIXELS"  BIGINT,
        "SPECTYPE"  TEXT,
        "SUBTYPE"  TEXT,
        "NCOEFF"  BIGINT,
        "DELTACHI2"  DOUBLE PRECISION,
        "NUMEXP"  INTEGER,
        "NUMTILE"  INTEGER,
        "PETAL_LOC"  SMALLINT,
        "DEVICE_LOC"  INTEGER,
        "LOCATION"  BIGINT,
        "FIBER"  INTEGER,
        "FIBERSTATUS"  INTEGER,
        "TARGET_RA"  DOUBLE PRECISION,
        "TARGET_DEC"  DOUBLE PRECISION,
        "PMRA"  REAL,
        "PMDEC"  REAL,
        "REF_EPOCH"  REAL,
        "LAMBDA_REF"  REAL,
        "FA_TARGET"  BIGINT,
        "FA_TYPE"  SMALLINT,
        "OBJTYPE"  TEXT,
        "PRIORITY"  INTEGER,
        "SUBPRIORITY"  DOUBLE PRECISION,
        "OBSCONDITIONS"  INTEGER,
        "RELEASE"  SMALLINT,
        "BRICKID"  INTEGER,
        "BRICK_OBJID"  INTEGER,
        "MORPHTYPE"  TEXT,
        "FLUX_G"  REAL,
        "FLUX_R"  REAL,
        "FLUX_Z"  REAL,
        "FLUX_IVAR_G"  REAL,
        "FLUX_IVAR_R"  REAL,
        "FLUX_IVAR_Z"  REAL,
        "MASKBITS"  SMALLINT,
        "REF_ID"  BIGINT,
        "REF_CAT"  TEXT,
        "GAIA_PHOT_G_MEAN_MAG"  REAL,
        "GAIA_PHOT_BP_MEAN_MAG"  REAL,
        "GAIA_PHOT_RP_MEAN_MAG"  REAL,
        "PARALLAX"  REAL,
        "BRICKNAME"  TEXT,
        "EBV"  REAL,
        "FLUX_W1"  REAL,
        "FLUX_W2"  REAL,
        "FIBERFLUX_G"  REAL,
        "FIBERFLUX_R"  REAL,
        "FIBERFLUX_Z"  REAL,
        "FIBERTOTFLUX_G"  REAL,
        "FIBERTOTFLUX_R"  REAL,
        "FIBERTOTFLUX_Z"  REAL,
        "SERSIC"  REAL,
        "SHAPE_R"  REAL,
        "SHAPE_E1"  REAL,
        "SHAPE_E2"  REAL,
        "PHOTSYS"  TEXT,
        "PRIORITY_INIT"  BIGINT,
        "NUMOBS_INIT"  BIGINT,
        "SV2_DESI_TARGET"  BIGINT,
        "SV2_BGS_TARGET"  BIGINT,
        "SV2_MWS_TARGET"  BIGINT,
        "SV2_SCND_TARGET"  BIGINT,
        "DESI_TARGET"  BIGINT,
        "BGS_TARGET"  BIGINT,
        "MWS_TARGET"  BIGINT,
        "TILEID"  INTEGER,
        "COADD_NUMEXP"  SMALLINT,
        "COADD_EXPTIME"  REAL,
        "MEAN_DELTA_X"  REAL,
        "RMS_DELTA_X"  REAL,
        "MEAN_DELTA_Y"  REAL,
        "RMS_DELTA_Y"  REAL,
        "MEAN_FIBER_X"  REAL,
        "MEAN_FIBER_Y"  REAL,
        "MEAN_FIBER_RA"  DOUBLE PRECISION,
        "MEAN_FIBER_DEC"  DOUBLE PRECISION,
        "MEAN_FIBERASSIGN_X"  REAL,
        "MEAN_FIBERASSIGN_Y"  REAL,
        "FIRST_NIGHT"  INTEGER,
        "LAST_NIGHT"  INTEGER,
        "NUM_NIGHT"  SMALLINT,
        "FIRST_EXPID"  INTEGER,
        "LAST_EXPID"  INTEGER,
        "NUM_EXPID"  SMALLINT,
        "FIRST_TILEID"  INTEGER,
        "LAST_TILEID"  INTEGER,
        "NUM_TILEID"  SMALLINT,
        "FIRST_FIBER"  INTEGER,
        "LAST_FIBER"  INTEGER,
        "NUM_FIBER"  SMALLINT,
        "FIRST_MJD"  REAL,
        "LAST_MJD"  REAL,
        "NUM_MJD"  SMALLINT,
        "SV1_SCND_TARGET"  BIGINT,
        "SV1_BGS_TARGET"  BIGINT,
        "SV1_DESI_TARGET"  BIGINT,
        "SV1_MWS_TARGET"  BIGINT,
        "CMX_TARGET"  BIGINT,
        "COEFF_0"  DOUBLE PRECISION,
        "COEFF_1"  DOUBLE PRECISION,
        "COEFF_2"  DOUBLE PRECISION,
        "COEFF_3"  DOUBLE PRECISION,
        "COEFF_4"  DOUBLE PRECISION,
        "COEFF_5"  DOUBLE PRECISION,
        "COEFF_6"  DOUBLE PRECISION,
        "COEFF_7"  DOUBLE PRECISION,
        "COEFF_8"  DOUBLE PRECISION,
        "COEFF_9"  DOUBLE PRECISION
    );
    """    
    
    @staticmethod
    def create_table(overwrite = False):
        
        with engine.connect() as conn:
#             cur = conn.cursor()

            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS zcatalog_denali_cumulative;"))
                             
            conn.execute(text(zcatalog_denali_cumulative.schema))
            conn.close()

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
#             df['PRODUCTION']=numpy.full(df.shape[0],prod)
#             df['COADD']=numpy.full(df.shape[0],coadd)

            try:
                df.to_sql('zcatalog_denali_cumulative',engine,index=False,if_exists='append') 
            except:
                dtypesToSchema(df.dtypes)
                print(df[0])
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

    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()
        
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
                    df.to_sql('dr9_pv',conn,index=False,if_exists='append')
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

    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()   
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS exposure_tables_daily;")
        cur.execute(exposure_tables_daily.schema)
        cur.close()

    @staticmethod
    def fill_table():      
        dir_root='/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/'
#         maxdate='20201214'

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
                dfs.to_sql(f'exposure_tables_daily',conn,index=False,if_exists='append')

            except sqlite3.OperationalError as err:
                dtypesToSchema(dfs.dtypes)
                sys.exit()                 

class proposals_pv:
    schema="""
        CREATE TABLE IF NOT EXISTS "proposals_pv" (
            "OBJID"  INTEGER,
            "BRICKID"  INTEGER,
            "BRICKNAME"  TEXT,
            "RA"  REAL,
            "DEC"  REAL,
            "PMRA"  TEXT,
            "PMDEC"  TEXT,
            "REF_EPOCH"  TEXT,
            "OVERRIDE"  TEXT,
            "PVTYPE"  TEXT,
            "PVPRIORITY"  INTEGER,
            "POINTINGID"  INTEGER,
            "SGA_ID"  INTEGER,
            "PRIORITY"  TEXT,
            "LUNATION"  TEXT,
            FOREIGN KEY (OBJID, BRICKID) REFERENCES DR9_PV (OBJID, BRICKID)
    );
    """
  
    
    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS proposals_pv;")
        cur.execute(proposals_pv.schema)
        cur.close()

    @staticmethod
    def fill_table():      
        runs = ["main_year1","sv3","sv1"]
        lunations = ["BRIGHT","DARK"]
        priorities = ["LOW", "MEDIUM", "HIGH"]
        for run in runs:
            for lunation in lunations:
                for priority in priorities:
                    filename = f"/global/cfs/cdirs/desi/target/proposals/proposals_{run}_frozen/indata/PV_{lunation}_{priority}.fits"
                    try:
                        dat = Table.read(filename, format='fits')
                    except:
                        print(f"{filename} not found")
                        continue                        
                    dat.convert_bytestring_to_unicode()
                    
                    df = dat.to_pandas()
                    df['PRIORITY']=numpy.full(df.shape[0],priority)
                    df['LUNATION']=numpy.full(df.shape[0],lunation)
                    try:
                        df.to_sql('proposals_pv',conn,index=False,if_exists='append')
                    except sqlite3.OperationalError as err:
                        dtypesToSchema(df.dtypes)
                        sys.exit()
        
    @staticmethod
    def fill_sga():

        cur = conn.cursor()
        cur.execute("UPDATE proposals_pv SET SGA_ID=(select dr9_pv.SGA_ID FROM dr9_pv WHERE proposals_pv.OBJID=dr9_pv.OBJID AND proposals_pv.BRICKID=dr9_pv.BRICKID);")
           
        cur.close()    


class zbest_daily:
    schema="""
        CREATE TABLE IF NOT EXISTS "zbest_daily" (
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
        "TILE"  INTEGER,
        "PETAL"  INTEGER,
        "YYYYMMDD"  INTEGER
    );
    """
  
    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS zbest_daily;")
        cur.execute(zbest_daily.schema)
        cur.close()

    @staticmethod
    def fill_table():      
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
     
        #find the last date
        dates_db=[]
        try:
            for row in conn.execute(f"SELECT DISTINCT YYYYMMDD FROM zbest_daily"):
                dates_db.append(row[0])
        except:
            pass
            
        dates = []
        for path in glob.glob(f'{dir_root}/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
        dates = numpy.unique(dates)
        for date in dates:
            if int(date) not in dates_db:
                # Do things in terms of dates
                dfs=[]
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        print(date,tile)
                        for i in range(10):                            
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
                            df['TILE']=numpy.full(df.shape[0],int(tile))
                            df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            dfs.append(df)
                if len(dfs)>0:
                    try:
                        dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                        dfs.to_sql(f'zbest_daily',conn,index=False,if_exists='append')
                    except sqlite3.OperationalError as err:
                        dtypesToSchema(dfs.dtypes)
                        sys.exit()
                        
                        
class fibermap_daily:
    schema="""
        CREATE TABLE IF NOT EXISTS "fibermap_daily" (
        "TARGETID"  INTEGER,
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
        "FIBERASSIGN_X"  TEXT,
        "FIBERASSIGN_Y"  TEXT,
        "NUMTARGET"  INTEGER,
        "PRIORITY"  INTEGER,
        "SUBPRIORITY"  REAL,
        "OBSCONDITIONS"  INTEGER,
        "MORPHTYPE"  TEXT,
        "FLUX_G"  TEXT,
        "FLUX_R"  TEXT,
        "FLUX_Z"  TEXT,
        "FLUX_IVAR_G"  TEXT,
        "FLUX_IVAR_R"  TEXT,
        "FLUX_IVAR_Z"  TEXT,
        "REF_ID"  INTEGER,
        "REF_CAT"  TEXT,
        "GAIA_PHOT_G_MEAN_MAG"  TEXT,
        "GAIA_PHOT_BP_MEAN_MAG"  TEXT,
        "GAIA_PHOT_RP_MEAN_MAG"  TEXT,
        "PARALLAX"  TEXT,
        "EBV"  TEXT,
        "FLUX_W1"  TEXT,
        "FLUX_W2"  TEXT,
        "FLUX_IVAR_W1"  TEXT,
        "FLUX_IVAR_W2"  TEXT,
        "FIBERFLUX_G"  TEXT,
        "FIBERFLUX_R"  TEXT,
        "FIBERFLUX_Z"  TEXT,
        "FIBERTOTFLUX_G"  TEXT,
        "FIBERTOTFLUX_R"  TEXT,
        "FIBERTOTFLUX_Z"  TEXT,
        "MASKBITS"  INTEGER,
        "SV1_DESI_TARGET"  INTEGER,
        "SV1_BGS_TARGET"  INTEGER,
        "SV1_MWS_TARGET"  INTEGER,
        "SV1_SCND_TARGET"  INTEGER,
        "SV2_DESI_TARGET"  INTEGER,
        "SV2_BGS_TARGET"  INTEGER,
        "SV2_MWS_TARGET"  INTEGER,
        "SV2_SCND_TARGET"  INTEGER,
        "SV3_DESI_TARGET"  INTEGER,
        "SV3_BGS_TARGET"  INTEGER,
        "SV3_MWS_TARGET"  INTEGER,
        "SV3_SCND_TARGET"  INTEGER,
        "SERSIC"  TEXT,
        "SHAPE_R"  TEXT,
        "SHAPE_E1"  TEXT,
        "SHAPE_E2"  TEXT,
        "PHOTSYS"  TEXT,
        "CMX_TARGET"  INTEGER,
        "PRIORITY_INIT"  INTEGER,
        "NUMOBS_INIT"  INTEGER,
        "RELEASE"  INTEGER,
        "BRICKID"  INTEGER,
        "BRICKNAME"  TEXT,
        "BRICK_OBJID"  INTEGER,
        "BLOBDIST"  TEXT,
        "FIBERFLUX_IVAR_G"  TEXT,
        "FIBERFLUX_IVAR_R"  TEXT,
        "FIBERFLUX_IVAR_Z"  TEXT,
        "DESI_TARGET"  INTEGER,
        "BGS_TARGET"  INTEGER,
        "MWS_TARGET"  INTEGER,
        "HPXPIXEL"  INTEGER,
        "NUM_ITER"  INTEGER,
        "FIBER_X"  REAL,
        "FIBER_Y"  REAL,
        "DELTA_X"  REAL,
        "DELTA_Y"  REAL,
        "FIBER_RA"  REAL,
        "FIBER_DEC"  REAL,
        "NIGHT"  INTEGER,
        "EXPID"  INTEGER,
        "MJD"  REAL,
        "EXPTIME"  REAL,
        "PSF_TO_FIBER_SPECFLUX"  REAL,
        "TILEID"  INTEGER,
        "YYYYMMDD"  INTEGER
    );
    """
  
    @staticmethod
    def create_table(overwrite = False):
        cur = conn.cursor()
        
        if overwrite:
            cur.execute("DROP TABLE IF EXISTS fibermap_daily;")
        cur.execute(fibermap_daily.schema)
        cur.close()

    @staticmethod
    def fill_table():      
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
     
        #find the last date
        dates_db=[]
        try:
            for row in conn.execute(f"SELECT DISTINCT YYYYMMDD FROM fibermap_daily"):
                dates_db.append(row[0])
        except:
            pass
            
        dates = []
        for path in glob.glob(f'{dir_root}/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
        dates = numpy.unique(dates)
        for date in dates:
            if int(date) not in dates_db:
                print(date)
                # Do things in terms of dates
                dfs=[]
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
#                         print(date,tile)
                        for i in range(10):                            
                            filename = f'{dir_root}/{tile}/{date}/zbest-{i}-{tile}-thru{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits',hdu=2)
                            except:
                                print(f"{filename} not found")
                                continue

                            dat.convert_bytestring_to_unicode()
                            df = dat.to_pandas()
#                             df['TILE']=numpy.full(df.shape[0],int(tile))
#                             df['PETAL']=numpy.full(df.shape[0],i)
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            dfs.append(df)
                if len(dfs)>0:
                    try:
                        dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                        dfs.to_sql(f'fibermap_daily',conn,index=False,if_exists='append')
                    except sqlite3.OperationalError as err:
                        dtypesToSchema(dfs.dtypes)
                        sys.exit()
                                                
        
        
class DBManager:
    
    filename = "/global/cfs/cdirs/desi/science/td/db/desi.db"

            
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
