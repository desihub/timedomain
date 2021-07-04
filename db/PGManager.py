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
        elif value == numpy.bool:
            nvalue = 'BOOLEAN'
        else:
            print ("**unknown** ",value)
#         print(f'{index}  {nvalue}, {value}')
        print(f'{index.lower()}  {nvalue},')       

class secondary:

    schema=f"""
    CREATE TABLE IF NOT EXISTS secondary (
        ra  DOUBLE PRECISION,
        dec  DOUBLE PRECISION,
        pmra  REAL,
        pmdec  REAL,
        ref_epoch  REAL,
        override  BOOLEAN,
        flux_g  REAL,
        flux_r  REAL,
        flux_z  REAL,
        parallax  REAL,
        gaia_phot_g_mean_mag  REAL,
        gaia_phot_bp_mean_mag  REAL,
        gaia_phot_rp_mean_mag  REAL,
        gaia_astrometric_excess_noise  REAL,
        targetid  BIGINT,
        sv3_desi_target  BIGINT,
        sv3_scnd_target  BIGINT,
        scnd_order  INTEGER,
        desitarget_v  TEXT,
        run  TEXT,
        program  TEXT,
        lunation  TEXT
    );
    """
    runs=['main','sv3']
    
    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS secondary;"))                            
            conn.execute(secondary.schema)
            conn.close()
        
    
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
                            df['desitarget_v']= numpy.full(df.shape[0],'0.57.0')
                            df['run']=numpy.full(df.shape[0],'sv3')
                            df['program']=numpy.full(df.shape[0],program)
                            df['lunation']=numpy.full(df.shape[0],lunation)
                            df.columns= df.columns.str.lower()
                            try:
                                df.to_sql('secondary',engine,index=False,if_exists='append')
                            except:
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
                            df['desitarget_v']= numpy.full(df.shape[0],'1.0.0')
                            df['run']=numpy.full(df.shape[0],'main')
                            df['program']=numpy.full(df.shape[0],program)
                            df['lunation']=numpy.full(df.shape[0],lunation)
                            df.columns= df.columns.str.lower()
                            try:
                                df.to_sql('secondary',engine,index=False,if_exists='append')
                            except:
                                dtypesToSchema(df.dtypes)
                                sys.exit()    
    
class zcatalog_denali:
    
    schema="""
    CREATE TABLE IF NOT EXISTS "zcatalog_denali" (
        targetid  BIGINT,
        chi2  DOUBLE PRECISION,
        z  DOUBLE PRECISION,
        zerr  DOUBLE PRECISION,
        zwarn  BIGINT,
        npixels  BIGINT,
        spectype  TEXT,
        subtype  TEXT,
        ncoeff  BIGINT,
        deltachi2  DOUBLE PRECISION,
        numexp  INTEGER,
        numtile  INTEGER,
        petal_loc  SMALLINT,
        device_loc  INTEGER,
        location  BIGINT,
        fiber  INTEGER,
        fiberstatus  INTEGER,
        target_ra  DOUBLE PRECISION,
        target_dec  DOUBLE PRECISION,
        pmra  REAL,
        pmdec  REAL,
        ref_epoch  REAL,
        lambda_ref  REAL,
        fa_target  BIGINT,
        fa_type  SMALLINT,
        objtype  TEXT,
        priority  INTEGER,
        subpriority  DOUBLE PRECISION,
        obsconditions  INTEGER,
        release  SMALLINT,
        brickid  INTEGER,
        brick_objid  INTEGER,
        morphtype  TEXT,
        flux_g  REAL,
        flux_r  REAL,
        flux_z  REAL,
        flux_ivar_g  REAL,
        flux_ivar_r  REAL,
        flux_ivar_z  REAL,
        maskbits  SMALLINT,
        ref_id  BIGINT,
        ref_cat  TEXT,
        gaia_phot_g_mean_mag  REAL,
        gaia_phot_bp_mean_mag  REAL,
        gaia_phot_rp_mean_mag  REAL,
        parallax  REAL,
        brickname  TEXT,
        ebv  REAL,
        flux_w1  REAL,
        flux_w2  REAL,
        fiberflux_g  REAL,
        fiberflux_r  REAL,
        fiberflux_z  REAL,
        fibertotflux_g  REAL,
        fibertotflux_r  REAL,
        fibertotflux_z  REAL,
        sersic  REAL,
        shape_r  REAL,
        shape_e1  REAL,
        shape_e2  REAL,
        photsys  TEXT,
        priority_init  BIGINT,
        numobs_init  BIGINT,
        sv2_desi_target  BIGINT,
        sv2_bgs_target  BIGINT,
        sv2_mws_target  BIGINT,
        sv2_scnd_target  BIGINT,
        desi_target  BIGINT,
        bgs_target  BIGINT,
        mws_target  BIGINT,
        tileid  INTEGER,
        coadd_numexp  SMALLINT,
        coadd_exptime  REAL,
        mean_delta_x  REAL,
        rms_delta_x  REAL,
        mean_delta_y  REAL,
        rms_delta_y  REAL,
        mean_fiber_x  REAL,
        mean_fiber_y  REAL,
        mean_fiber_ra  DOUBLE PRECISION,
        mean_fiber_dec  DOUBLE PRECISION,
        mean_fiberassign_x  REAL,
        mean_fiberassign_y  REAL,
        first_night  INTEGER,
        last_night  INTEGER,
        num_night  SMALLINT,
        first_expid  INTEGER,
        last_expid  INTEGER,
        num_expid  SMALLINT,
        first_tileid  INTEGER,
        last_tileid  INTEGER,
        num_tileid  SMALLINT,
        first_fiber  INTEGER,
        last_fiber  INTEGER,
        num_fiber  SMALLINT,
        first_mjd  REAL,
        last_mjd  REAL,
        num_mjd  SMALLINT,
        sv1_scnd_target  BIGINT,
        sv1_bgs_target  BIGINT,
        sv1_desi_target  BIGINT,
        sv1_mws_target  BIGINT,
        cmx_target  BIGINT,
        coeff_0  DOUBLE PRECISION,
        coeff_1  DOUBLE PRECISION,
        coeff_2  DOUBLE PRECISION,
        coeff_3  DOUBLE PRECISION,
        coeff_4  DOUBLE PRECISION,
        coeff_5  DOUBLE PRECISION,
        coeff_6  DOUBLE PRECISION,
        coeff_7  DOUBLE PRECISION,
        coeff_8  DOUBLE PRECISION,
        coeff_9  DOUBLE PRECISION,
        PRIMARY KEY (targetid, tileid)
    );
    """    
    
    # 
    
    @staticmethod
    def create_table(overwrite = False):
        
        with engine.connect() as conn:
            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS zcatalog_denali;"))                            
            conn.execute(text(zcatalog_denali.schema))
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
            df.columns= df.columns.str.lower()
            try:
                df.to_sql('zcatalog_denali',engine,index=False,if_exists='append') 
            except:
                dtypesToSchema(df.dtypes)
                print(df[0])
                sys.exit()
class mtl:

    schema=f"""
        CREATE TABLE IF NOT EXISTS "mtl" (
            ra  DOUBLE PRECISION,
            dec  DOUBLE PRECISION,
            ref_epoch  REAL,
            parallax  REAL,
            pmra  REAL,
            pmdec  REAL,
            targetid  BIGINT,
            sv3_desi_target  BIGINT,
            sv3_bgs_target  BIGINT,
            sv3_mws_target  BIGINT,
            subpriority  DOUBLE PRECISION,
            obsconditions  INTEGER,
            priority_init  BIGINT,
            numobs_init  BIGINT,
            sv3_scnd_target  BIGINT,
            numobs_more  BIGINT,
            numobs  BIGINT,
            z  DOUBLE PRECISION,
            zwarn  BIGINT,
            ztileid  INTEGER,
            target_state  TEXT,
            timestamp  TEXT,
            version  TEXT,
            priority  BIGINT,
            run  TEXT,
            program  TEXT,
            lunation  TEXT,
            flux_g  REAL,
            flux_r  REAL,
            flux_z  REAL,
            gaia_phot_g_mean_mag  REAL,
            gaia_phot_bp_mean_mag  REAL,
            gaia_phot_rp_mean_mag  REAL,
            gaia_astrometric_excess_noise  REAL,
            scnd_order  INTEGER,
            checker  TEXT,
            too_type  TEXT,
            too_prio  TEXT,
            oclayer  TEXT,
            mjd_begin  DOUBLE PRECISION,
            mjd_end  DOUBLE PRECISION,
            tooid  BIGINT
        );
        """

    root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
    runs = ["main","sv3","sv2"]
    runs=["sv3"]
    programs = ["ToO","","secondary"]
    lunations = ["bright","dark",""] 

    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(f"DROP TABLE IF EXISTS mtl;")
            conn.execute(mtl.schema)
            conn.close()

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
                            df['run']=numpy.full(df.shape[0],run)
                            df['program']=numpy.full(df.shape[0],program)
                            df['lunation']=numpy.full(df.shape[0],lunation)
                            df.columns= df.columns.str.lower()
                            try:
                                df.to_sql('mtl',engine,index=False,if_exists='append')
                            except:
                                dtypesToSchema(df.dtypes)
                                print(df.dtypes)
                                sys.exit()
                    else:
                        print (path, 'not exists')


                                


                
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
