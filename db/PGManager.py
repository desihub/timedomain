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
    
fibermap_daily:

    read from cumulative only those taken that night

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
    
class zbest_daily:
    schema="""
        CREATE TABLE IF NOT EXISTS "zbest_daily" (
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
        tile  BIGINT,
        petal  BIGINT,
        yyyymmdd  BIGINT,
        PRIMARY KEY (targetid, tile, yyyymmdd)
    );
    """
  
    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS zbest_daily;"))                            
            conn.execute(zbest_daily.schema)
            conn.close()

    @staticmethod
    def fill_table():      
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
     
        #find the last date
        dates_db=[]
        try:
            with engine.connect() as conn:
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
                print(date)
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
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
                            df.columns= df.columns.str.lower()
                            dfs.append(df)
                if len(dfs)>0:
                    try:
                        dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                        dfs.to_sql(f'zbest_daily',engine,index=False,if_exists='append')
                    except sqlite3.OperationalError as err:
                        dtypesToSchema(dfs.dtypes)
                        sys.exit()

class dr9_pv:
    schema="""
        CREATE TABLE IF NOT EXISTS "dr9_pv" (
        objid  BIGINT,
        brickid  INTEGER,
        brickname  TEXT,
        ra  DOUBLE PRECISION,
        dec  DOUBLE PRECISION,
        type  TEXT,
        sersic  REAL,
        z_phot_median  REAL,
        z_phot_l95  REAL,
        mag_g  REAL,
        mag_r  REAL,
        mag_z  REAL,
        mag_b  REAL,
        mag_g_err  REAL,
        mag_r_err  REAL,
        mag_z_err  REAL,
        fibre_mag_g  REAL,
        fibre_mag_r  REAL,
        fibre_mag_z  REAL,
        uncor_radius  REAL,
        ba_ratio  REAL,
        circ_radius  REAL,
        pos_angle  REAL,
        insga  BOOLEAN,
        inbgs  BOOLEAN,
        inlocalbright  BOOLEAN,
        inspecfootprint  BOOLEAN,
        sga_pa  REAL,
        sga_ba  REAL,
        sb_d25_g  REAL,
        sb_d25_r  REAL,
        sb_d25_z  REAL,
        radius_sb25  REAL,
        sga_morphtype  TEXT,
        sga_id  BIGINT,
        sga_redshift  REAL,
        size_sga  REAL,
        pmra  REAL,
        pmdec  REAL,
        ref_epoch  REAL,
        override  BOOLEAN,
        pvtype  TEXT,
        pvpriority  INTEGER,
        target  TEXT,
        survey  TEXT,
        pointingid  INTEGER
    );
    """

    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS dr9_pv;"))                            
            conn.execute(dr9_pv.schema)
            conn.close()
            
    @staticmethod
    def fill_table():      
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
                df['target']=numpy.full(df.shape[0],target)
                df['survey']=numpy.full(df.shape[0],survey)
                df.columns= df.columns.str.lower()
                try:
                    df.to_sql('dr9_pv',engine,index=False,if_exists='append')
                except:
                    dtypesToSchema(df.dtypes)
                    sys.exit()
    
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
        lunation  TEXT,
        desi_target  BIGINT,
        scnd_target  BIGINT
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


      

class exposure_tables_daily:
    schema="""
        CREATE TABLE IF NOT EXISTS "exposure_tables_daily" (
        expid  BIGINT PRIMARY KEY,
        exptime  DOUBLE PRECISION,
        obstype  TEXT,
        spectrographs  DOUBLE PRECISION,
        camword  TEXT,
        tileid  BIGINT,
        night  BIGINT,
        expflag  TEXT,
        headererr  TEXT,
        survey  TEXT,
        seqnum  BIGINT,
        seqtot  BIGINT,
        program  TEXT,
        mjd_obs  DOUBLE PRECISION,
        reqra  DOUBLE PRECISION,
        reqdec  DOUBLE PRECISION,
        targtra  DOUBLE PRECISION,
        targtdec  DOUBLE PRECISION,
        comments  TEXT,
        yyyymmdd  BIGINT,
        purpose  TEXT,
        fa_surv  TEXT,
        badcamword  TEXT,
        badamps  TEXT,
        laststep  TEXT,
        efftime_etc  DOUBLE PRECISION,
        faprgrm  TEXT,
        goaltime  DOUBLE PRECISION,
        goaltype  TEXT,
        ebvfac  DOUBLE PRECISION,
        airfac  DOUBLE PRECISION,
        speed  DOUBLE PRECISION
    );
    """

    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(f"DROP TABLE IF EXISTS exposure_tables_daily;")
            conn.execute(exposure_tables_daily.schema)
            conn.close()

    @staticmethod
    def fill_table():      
        dir_root='/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/'

        #find the last date
        dates_db=[]
        try:
            with engine.connect() as conn:
                for row in conn.execute(f"SELECT DISTINCT yyyymmdd FROM exposure_tables_daily"):
                    dates_db.append(row[0])
        except:
            pass
        
        dates = []
        for path in glob.glob(f'{dir_root}/*/exposure_table_????????.csv'):
            split = path.split('/')
            date = split[-1][-12:-4]
            dates.append(date)

        dates=numpy.sort(dates)
#         dates=numpy.flip(dates)

        dfs=[]
        for date in dates:
            if int(date) not in dates_db:
                file = f'/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/{date[0:6]}/exposure_table_{date}.csv'
                data = ascii.read(file)
                df = data.to_pandas()
                df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                df.columns= df.columns.str.lower()
                dfs.append(df)
                

          
        if len(dfs)>0:
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)
            dfs.rename(columns={"mjd-obs": "mjd_obs"}, inplace=True)
                
            df.columns= df.columns.str.lower()
            try:
                dfs.to_sql(f'exposure_tables_daily',engine,index=False,if_exists='append')

            except:
                dtypesToSchema(dfs.dtypes)
                sys.exit()                 

class proposals_pv:
    schema="""
        CREATE TABLE IF NOT EXISTS "proposals_pv" (
            objid  BIGINT,
            brickid  INTEGER,
            brickname  TEXT,
            ra  DOUBLE PRECISION,
            dec  DOUBLE PRECISION,
            pmra  REAL,
            pmdec  REAL,
            ref_epoch  REAL,
            override  BOOLEAN,
            pvtype  TEXT,
            pvpriority  INTEGER,
            pointingid  INTEGER,
            sga_id  BIGINT,
            priority  TEXT,
            lunation  TEXT
    );
    """
  
    
    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS proposals_pv;"))                            
            conn.execute(proposals_pv.schema)
            conn.close()
            
            
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
                    df['priority']=numpy.full(df.shape[0],priority)
                    df['lunation']=numpy.full(df.shape[0],lunation)
                    df.columns= df.columns.str.lower()
                    try:
                        df.to_sql('proposals_pv',engine,index=False,if_exists='append')
                    except sqlite3.OperationalError as err:
                        dtypesToSchema(df.dtypes)
                        sys.exit()
        
    @staticmethod
    def fill_sga():

        cur = conn.cursor()
        cur.execute("UPDATE proposals_pv SET SGA_ID=(select dr9_pv.SGA_ID FROM dr9_pv WHERE proposals_pv.OBJID=dr9_pv.OBJID AND proposals_pv.BRICKID=dr9_pv.BRICKID);")
           
        cur.close()    



                        
                        
class fibermap_daily:
    schema="""
        CREATE TABLE IF NOT EXISTS "fibermap_daily" (
            targetid  BIGINT,
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
            fiberassign_x  REAL,
            fiberassign_y  REAL,
            numtarget  SMALLINT,
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
            flux_ivar_w1  REAL,
            flux_ivar_w2  REAL,
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
            cmx_target  BIGINT,
            sv1_desi_target  BIGINT,
            sv1_bgs_target  BIGINT,
            sv1_mws_target  BIGINT,
            priority_init  BIGINT,
            numobs_init  BIGINT,
            desi_target  BIGINT,
            bgs_target  BIGINT,
            mws_target  BIGINT,
            scnd_target  BIGINT,
            plate_ra  DOUBLE PRECISION,
            plate_dec  DOUBLE PRECISION,
            num_iter  BIGINT,
            fiber_x  DOUBLE PRECISION,
            fiber_y  DOUBLE PRECISION,
            delta_x  DOUBLE PRECISION,
            delta_y  DOUBLE PRECISION,
            fiber_ra  DOUBLE PRECISION,
            fiber_dec  DOUBLE PRECISION,
            exptime  DOUBLE PRECISION,
            psf_to_fiber_specflux  DOUBLE PRECISION,
            night  INTEGER,
            expid  INTEGER,
            mjd  DOUBLE PRECISION,
            tileid  INTEGER,
            yyyymmdd  BIGINT,
            blobdist  REAL,
            fiberflux_ivar_g  REAL,
            fiberflux_ivar_r  REAL,
            fiberflux_ivar_z  REAL,
            hpxpixel  BIGINT
    );
    """
  
    @staticmethod
    def create_table(overwrite = False):
        with engine.connect() as conn:
            if overwrite:
                conn.execute(text("DROP TABLE IF EXISTS fibermap_daily;"))                            
            conn.execute(fibermap_daily.schema)
            conn.close()


    @staticmethod
    def fill_table():      
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
     
        #find the last date
        dates_db=[]
        try:
            with engine.connect() as conn:
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
                            df=df[df['NIGHT']== int(date)]
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            
                            df.columns= df.columns.str.lower()
                            dfs.append(df)
                if len(dfs)>0:
                    try:
                        dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                        dfs.to_sql(f'fibermap_daily',engine,index=False,if_exists='append')
                    except:
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
