import os
import sys
import glob
import argparse
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

# try:
#     with open('/global/homes/a/akim/secrets/desi_pg.txt', 'r') as f:
#         data = f.read()
#         uname, passwd = data.split()
#         engine = create_engine(f'postgresql://{uname}:{passwd}@decatdb.lbl.gov/desi')
# except:
#     print("No connection")

# #     conn = sqlite3.connect("/global/cfs/cdirs/desi/science/td/db/desi.db") 

# global variable engine is the engine connection

engine = None
def create_db_engine( secrets_file, dbhost="decatdb.lbl.gov", dbname="desi" ):
    global engine
    with open( secrets_file, 'r' ) as f:
        data = f.read()
        uname, passwd = data.split()
    engine = create_engine( f'postgresql://{uname}:{passwd}@{dbhost}/{dbname}' )

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


class tablebaseclass:
    """A base class for tables that has some methods that are shared by many"""

    tablename = None
    createstatements = None

    @classmethod
    def create_table( cls, overwrite=False ):
        global engine
        with engine.connect() as conn:
            if overwrite:
                sys.stderr.write( f"Dropping {cls.tablename}\n" )
                conn.execute( f"DROP TABLE IF EXISTS {cls.tablename}" )
            sys.stderr.write( f"Creating {cls.tablename}...\n" )
            for statement in cls.createstatements:
                conn.execute( statement )
            sys.stderr.write( "...done.\n" )
            conn.close()
    
        
class zbest_daily(tablebaseclass):
    tablename = "zbest_daily"
    createstatements = [ """
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
                        ]
  
    @staticmethod
    def fill_table():
        global engine
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
                    
                if date < '20210901':
                    base='zbest'
                else:
                    base='redrock'
                
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        for i in range(10):                            
                            filename = f'{dir_root}/{tile}/{date}/{base}-{i}-{tile}-thru{date}.fits'
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

class dr9_pv(tablebaseclass):
    tablename = "pv.dr9"
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS "pv.dr9" (
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
                        ]

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
                    df.to_sql('pv.dr9',engine,index=False,if_exists='append')
                except:
                    dtypesToSchema(df.dtypes)
                    sys.exit()
                    
    @staticmethod
    def null_sga():
        global engine
        with engine.connect() as conn:
#             conn.execute("create table temp as (select * from dr9_pv);")
            query = """
            UPDATE dr9_pv
            SET sga_id = NULL,
                sga_morphtype  = NULL,
                sga_redshift  = NULL,
                sga_pa  = NULL,
                sga_ba = NULL
            WHERE sga_id = 0
            OR sga_id = -1;
            """
            conn.execute(query)
    
class secondary(tablebaseclass):
    tablename = "secondary"
    createstatements = [ """
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
        )
    """,
                         """
       CREATE INDEX IF NOT EXISTS secondary_q3c_ang2ipix_idx
          ON secondary USING btree(public.q3c_ang2ipix(ra, "dec"))
    """,
                        ]
    runs=['main','sv3']
    
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
    
class zcatalog_denali(tablebaseclass):
    tablename = "zcatalog_denali"
    createstatements = [ """
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
        )
    """
                         ]
    
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
                
class zpix_everest(tablebaseclass):
    """Superclass for everest_zpix_redshifts and everest_zpix_fiberman

    This class doesn't represent a table itself.
    """

    @classmethod
    def fill_table(cls, prod='everest', schema=False):

        print(cls.tablename)

        root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/zcatalog/"
        coadds = ["cumulative"]
        coadds=[""]
        runs=['main','sv1','sv2','sv3']
        programs = ['bright','dark','backup','other']
        
        dfs=[]
        for coadd in coadds:
            for run in runs:
                for program in programs:
                    filename = os.path.join(root,f"zpix-{run}-{program}.fits")
                    print(filename)
                    try:
                        dat = Table.read(filename, format='fits',hdu=cls.hdu)
                    except:
                        continue
                        
                    print(filename)
                    
                    # There is a multidimensional column that needs to be broken up
                    self.handle_coef( dat )

                    dat.convert_bytestring_to_unicode()

                    df = dat.to_pandas()
                    df['run']=numpy.full(df.shape[0],run)
                    df['program']=numpy.full(df.shape[0],program)
                    df.columns= df.columns.str.lower()
#                     dtypesToSchema(df.dtypes)

                    if schema:
                        dfs.append(df[0:1])
            
                    if not schema:
                        try:
                            df.to_sql(cls.tablename,engine,index=False,if_exists='append',schema=prod) 
                        except:
                            dtypesToSchema(df.dtypes)
                            print(df[0])
                            sys.exit()
        if schema:
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)
            dtypesToSchema(dfs.dtypes)


    def handle_coef( dat ):
        pass
            
class everest_zpix_redshifts(zpix_everest):
    tablename = "everest.zpix_redshifts"
    hdu = 1
    prod='everest'
    
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS everest.zpix_redshifts (
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
            coadd_fiberstatus  INTEGER,
            target_ra  DOUBLE PRECISION,
            target_dec  DOUBLE PRECISION,
            pmra  REAL,
            pmdec  REAL,
            ref_epoch  REAL,
            fa_target  BIGINT,
            fa_type  SMALLINT,
            objtype  TEXT,
            subpriority  DOUBLE PRECISION,
            obsconditions  INTEGER,
            release  INTEGER,
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
            priority_init  BIGINT,
            numobs_init  BIGINT,
            desi_target  BIGINT,
            bgs_target  BIGINT,
            mws_target  BIGINT,
            scnd_target  DOUBLE PRECISION,
            plate_ra  DOUBLE PRECISION,
            plate_dec  DOUBLE PRECISION,
            coadd_numexp  SMALLINT,
            coadd_exptime  REAL,
            coadd_numnight  SMALLINT,
            coadd_numtile  SMALLINT,
            mean_delta_x  REAL,
            rms_delta_x  REAL,
            mean_delta_y  REAL,
            rms_delta_y  REAL,
            mean_fiber_ra  DOUBLE PRECISION,
            std_fiber_ra  REAL,
            mean_fiber_dec  DOUBLE PRECISION,
            std_fiber_dec  REAL,
            mean_psf_to_fiber_specflux  REAL,
            tsnr2_gpbdark_b  REAL,
            tsnr2_elg_b  REAL,
            tsnr2_gpbbright_b  REAL,
            tsnr2_lya_b  REAL,
            tsnr2_bgs_b  REAL,
            tsnr2_gpbbackup_b  REAL,
            tsnr2_qso_b  REAL,
            tsnr2_lrg_b  REAL,
            tsnr2_gpbdark_r  REAL,
            tsnr2_elg_r  REAL,
            tsnr2_gpbbright_r  REAL,
            tsnr2_lya_r  REAL,
            tsnr2_bgs_r  REAL,
            tsnr2_gpbbackup_r  REAL,
            tsnr2_qso_r  REAL,
            tsnr2_lrg_r  REAL,
            tsnr2_gpbdark_z  REAL,
            tsnr2_elg_z  REAL,
            tsnr2_gpbbright_z  REAL,
            tsnr2_lya_z  REAL,
            tsnr2_bgs_z  REAL,
            tsnr2_gpbbackup_z  REAL,
            tsnr2_qso_z  REAL,
            tsnr2_lrg_z  REAL,
            tsnr2_gpbdark  REAL,
            tsnr2_elg  REAL,
            tsnr2_gpbbright  REAL,
            tsnr2_lya  REAL,
            tsnr2_bgs  REAL,
            tsnr2_gpbbackup  REAL,
            tsnr2_qso  REAL,
            tsnr2_lrg  REAL,
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
            run  TEXT,
            program  TEXT,
            sv1_desi_target  DOUBLE PRECISION,
            sv1_bgs_target  DOUBLE PRECISION,
            sv1_mws_target  DOUBLE PRECISION,
            sv1_scnd_target  DOUBLE PRECISION,
            numtarget  DOUBLE PRECISION,
            blobdist  REAL,
            fiberflux_ivar_g  REAL,
            fiberflux_ivar_r  REAL,
            fiberflux_ivar_z  REAL,
            hpxpixel  DOUBLE PRECISION,
            sv2_desi_target  DOUBLE PRECISION,
            sv2_bgs_target  DOUBLE PRECISION,
            sv2_mws_target  DOUBLE PRECISION,
            sv2_scnd_target  DOUBLE PRECISION,
            sv3_desi_target  DOUBLE PRECISION,
            sv3_bgs_target  DOUBLE PRECISION,
            sv3_mws_target  DOUBLE PRECISION,
            sv3_scnd_target  DOUBLE PRECISION
        )
    """,
                         """
        CREATE INDEX IF NOT EXISTS zpix_redshifts_idx_targetid
            ON everest.zpix_redshifts USING btree (targetid)
    """,
                         """
        CREATE INDEX IF NOT EXISTS zpix_redshifts_q3c_ang2ipix_idx
            ON everest.zpix_redshifts USING btree(public.q3c_ang2ipix(target_ra, target_dec))
    """,
                        ]

    def handle_coef( self, dat ):
        for icoeff in range(0,10):
            dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
        dat.remove_column('COEFF')
        

class everest_zpix_fibermap(zpix_everest):
    tablename = "everest.zpix_fibermap"
    hdu = 2
    prod = 'everest'
    
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS everest.zpix_fibermap (    
            targetid  BIGINT,
            priority  INTEGER,
            subpriority  DOUBLE PRECISION,
            night  INTEGER,
            expid  INTEGER,
            mjd  DOUBLE PRECISION,
            tileid  INTEGER,
            exptime  DOUBLE PRECISION,
            petal_loc  SMALLINT,
            device_loc  INTEGER,
            location  BIGINT,
            fiber  INTEGER,
            fiberstatus  INTEGER,
            fiberassign_x  REAL,
            fiberassign_y  REAL,
            lambda_ref  REAL,
            plate_ra  DOUBLE PRECISION,
            plate_dec  DOUBLE PRECISION,
            num_iter  BIGINT,
            fiber_x  DOUBLE PRECISION,
            fiber_y  DOUBLE PRECISION,
            delta_x  DOUBLE PRECISION,
            delta_y  DOUBLE PRECISION,
            fiber_ra  DOUBLE PRECISION,
            fiber_dec  DOUBLE PRECISION,
            psf_to_fiber_specflux  DOUBLE PRECISION,
            run  TEXT,
            program  TEXT
        )
    """,
                         """
        CREATE INDEX IF NOT EXISTS zpix_fibermap_idx_targetid
            ON everest.zpix_fibermap USING btree (targetid)
    """,
                         """
        CREATE INDEX IF NOT EXISTS zpix_fibermap_q3c_ang2ipix_idx
            ON everest.zpix_fibermap USING btree(public.q3c_ang2ipix(fiber_ra, fiber_dec))
    """,
                         ]
                         

class ztile_everest:
    
    @staticmethod
    def tablename(prod='everest',coadd='pernight',hdu=1,short=False):
        if hdu==1:
            tag='redshifts'
        elif hdu==2:
            tag='fibermap'
        if not short:    
            return f"{prod}.ztile_{coadd}_{tag}"
        else:
            return f"ztile_{coadd}_{tag}"
    
    
    @staticmethod
    def create_table(prod='everest',overwrite = False, hdu=1,coadd='pernight'):
        global engine
        tablename = ztile_everest.tablename(prod=prod,coadd=coadd,hdu=hdu)
        if hdu==1 and coadd=='pernight':
            schema="""
                CREATE TABLE IF NOT EXISTS {} (
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
                    petal_loc  SMALLINT,
                    device_loc  INTEGER,
                    location  BIGINT,
                    fiber  INTEGER,
                    coadd_fiberstatus  INTEGER,
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
                    priority  INTEGER,
                    subpriority  DOUBLE PRECISION,
                    obsconditions  INTEGER,
                    release  INTEGER,
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
                    priority_init  BIGINT,
                    numobs_init  BIGINT,
                    desi_target  BIGINT,
                    bgs_target  BIGINT,
                    mws_target  BIGINT,
                    scnd_target  DOUBLE PRECISION,
                    plate_ra  DOUBLE PRECISION,
                    plate_dec  DOUBLE PRECISION,
                    tileid  INTEGER,
                    coadd_numexp  SMALLINT,
                    coadd_exptime  REAL,
                    coadd_numnight  SMALLINT,
                    coadd_numtile  SMALLINT,
                    mean_delta_x  REAL,
                    rms_delta_x  REAL,
                    mean_delta_y  REAL,
                    rms_delta_y  REAL,
                    mean_fiber_ra  DOUBLE PRECISION,
                    std_fiber_ra  REAL,
                    mean_fiber_dec  DOUBLE PRECISION,
                    std_fiber_dec  REAL,
                    mean_psf_to_fiber_specflux  REAL,
                    mean_fiber_x  REAL,
                    mean_fiber_y  REAL,
                    tsnr2_gpbdark_b  REAL,
                    tsnr2_elg_b  REAL,
                    tsnr2_gpbbright_b  REAL,
                    tsnr2_lya_b  REAL,
                    tsnr2_bgs_b  REAL,
                    tsnr2_gpbbackup_b  REAL,
                    tsnr2_qso_b  REAL,
                    tsnr2_lrg_b  REAL,
                    tsnr2_gpbdark_r  REAL,
                    tsnr2_elg_r  REAL,
                    tsnr2_gpbbright_r  REAL,
                    tsnr2_lya_r  REAL,
                    tsnr2_bgs_r  REAL,
                    tsnr2_gpbbackup_r  REAL,
                    tsnr2_qso_r  REAL,
                    tsnr2_lrg_r  REAL,
                    tsnr2_gpbdark_z  REAL,
                    tsnr2_elg_z  REAL,
                    tsnr2_gpbbright_z  REAL,
                    tsnr2_lya_z  REAL,
                    tsnr2_bgs_z  REAL,
                    tsnr2_gpbbackup_z  REAL,
                    tsnr2_qso_z  REAL,
                    tsnr2_lrg_z  REAL,
                    tsnr2_gpbdark  REAL,
                    tsnr2_elg  REAL,
                    tsnr2_gpbbright  REAL,
                    tsnr2_lya  REAL,
                    tsnr2_bgs  REAL,
                    tsnr2_gpbbackup  REAL,
                    tsnr2_qso  REAL,
                    tsnr2_lrg  REAL,
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
                    run  TEXT,
                    program  TEXT,
                    sv1_desi_target  DOUBLE PRECISION,
                    sv1_bgs_target  DOUBLE PRECISION,
                    sv1_mws_target  DOUBLE PRECISION,
                    sv1_scnd_target  DOUBLE PRECISION,
                    numtarget  DOUBLE PRECISION,
                    blobdist  REAL,
                    fiberflux_ivar_g  REAL,
                    fiberflux_ivar_r  REAL,
                    fiberflux_ivar_z  REAL,
                    hpxpixel  DOUBLE PRECISION,
                    sv2_desi_target  DOUBLE PRECISION,
                    sv2_bgs_target  DOUBLE PRECISION,
                    sv2_mws_target  DOUBLE PRECISION,
                    sv2_scnd_target  DOUBLE PRECISION,
                    sv3_desi_target  DOUBLE PRECISION,
                    sv3_bgs_target  DOUBLE PRECISION,
                    sv3_mws_target  DOUBLE PRECISION,
                    sv3_scnd_target  DOUBLE PRECISION
                );
                """.format(tablename)
        elif hdu==2 and coadd=='pernight':
            schema="""
                CREATE TABLE IF NOT EXISTS {} (    
                    targetid  BIGINT,
                    priority  INTEGER,
                    subpriority  DOUBLE PRECISION,
                    night  INTEGER,
                    expid  INTEGER,
                    mjd  DOUBLE PRECISION,
                    tileid  INTEGER,
                    exptime  DOUBLE PRECISION,
                    petal_loc  SMALLINT,
                    device_loc  INTEGER,
                    location  BIGINT,
                    fiber  INTEGER,
                    fiberstatus  INTEGER,
                    fiberassign_x  REAL,
                    fiberassign_y  REAL,
                    lambda_ref  REAL,
                    plate_ra  DOUBLE PRECISION,
                    plate_dec  DOUBLE PRECISION,
                    num_iter  BIGINT,
                    fiber_x  DOUBLE PRECISION,
                    fiber_y  DOUBLE PRECISION,
                    delta_x  DOUBLE PRECISION,
                    delta_y  DOUBLE PRECISION,
                    fiber_ra  DOUBLE PRECISION,
                    fiber_dec  DOUBLE PRECISION,
                    psf_to_fiber_specflux  DOUBLE PRECISION,
                    run  TEXT,
                    program  TEXT
                );
                """.format(tablename)
        elif hdu==1 and coadd=='cumulative':
            schema="""
                CREATE TABLE IF NOT EXISTS {} (    
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
                petal_loc  SMALLINT,
                device_loc  INTEGER,
                location  BIGINT,
                fiber  INTEGER,
                coadd_fiberstatus  INTEGER,
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
                priority  INTEGER,
                subpriority  DOUBLE PRECISION,
                obsconditions  INTEGER,
                release  INTEGER,
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
                priority_init  BIGINT,
                numobs_init  BIGINT,
                desi_target  BIGINT,
                bgs_target  BIGINT,
                mws_target  BIGINT,
                scnd_target  DOUBLE PRECISION,
                plate_ra  DOUBLE PRECISION,
                plate_dec  DOUBLE PRECISION,
                tileid  INTEGER,
                coadd_numexp  SMALLINT,
                coadd_exptime  REAL,
                coadd_numnight  SMALLINT,
                coadd_numtile  SMALLINT,
                mean_delta_x  REAL,
                rms_delta_x  REAL,
                mean_delta_y  REAL,
                rms_delta_y  REAL,
                mean_fiber_ra  DOUBLE PRECISION,
                std_fiber_ra  REAL,
                mean_fiber_dec  DOUBLE PRECISION,
                std_fiber_dec  REAL,
                mean_psf_to_fiber_specflux  REAL,
                mean_fiber_x  REAL,
                mean_fiber_y  REAL,
                tsnr2_gpbdark_b  REAL,
                tsnr2_elg_b  REAL,
                tsnr2_gpbbright_b  REAL,
                tsnr2_lya_b  REAL,
                tsnr2_bgs_b  REAL,
                tsnr2_gpbbackup_b  REAL,
                tsnr2_qso_b  REAL,
                tsnr2_lrg_b  REAL,
                tsnr2_gpbdark_r  REAL,
                tsnr2_elg_r  REAL,
                tsnr2_gpbbright_r  REAL,
                tsnr2_lya_r  REAL,
                tsnr2_bgs_r  REAL,
                tsnr2_gpbbackup_r  REAL,
                tsnr2_qso_r  REAL,
                tsnr2_lrg_r  REAL,
                tsnr2_gpbdark_z  REAL,
                tsnr2_elg_z  REAL,
                tsnr2_gpbbright_z  REAL,
                tsnr2_lya_z  REAL,
                tsnr2_bgs_z  REAL,
                tsnr2_gpbbackup_z  REAL,
                tsnr2_qso_z  REAL,
                tsnr2_lrg_z  REAL,
                tsnr2_gpbdark  REAL,
                tsnr2_elg  REAL,
                tsnr2_gpbbright  REAL,
                tsnr2_lya  REAL,
                tsnr2_bgs  REAL,
                tsnr2_gpbbackup  REAL,
                tsnr2_qso  REAL,
                tsnr2_lrg  REAL,
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
                run  TEXT,
                program  TEXT,
                sv1_desi_target  DOUBLE PRECISION,
                sv1_bgs_target  DOUBLE PRECISION,
                sv1_mws_target  DOUBLE PRECISION,
                sv1_scnd_target  DOUBLE PRECISION,
                numtarget  DOUBLE PRECISION,
                blobdist  REAL,
                fiberflux_ivar_g  REAL,
                fiberflux_ivar_r  REAL,
                fiberflux_ivar_z  REAL,
                hpxpixel  DOUBLE PRECISION,
                sv2_desi_target  DOUBLE PRECISION,
                sv2_bgs_target  DOUBLE PRECISION,
                sv2_mws_target  DOUBLE PRECISION,
                sv2_scnd_target  DOUBLE PRECISION,
                sv3_desi_target  DOUBLE PRECISION,
                sv3_bgs_target  DOUBLE PRECISION,
                sv3_mws_target  DOUBLE PRECISION,
                sv3_scnd_target  DOUBLE PRECISION
                );
                """.format(tablename)
        elif hdu==2 and coadd=='cumulative':
            schema="""
                CREATE TABLE IF NOT EXISTS {} (    
                targetid  BIGINT,
                priority  INTEGER,
                subpriority  DOUBLE PRECISION,
                night  INTEGER,
                expid  INTEGER,
                mjd  DOUBLE PRECISION,
                tileid  INTEGER,
                exptime  DOUBLE PRECISION,
                petal_loc  SMALLINT,
                device_loc  INTEGER,
                location  BIGINT,
                fiber  INTEGER,
                fiberstatus  INTEGER,
                fiberassign_x  REAL,
                fiberassign_y  REAL,
                lambda_ref  REAL,
                plate_ra  DOUBLE PRECISION,
                plate_dec  DOUBLE PRECISION,
                num_iter  BIGINT,
                fiber_x  DOUBLE PRECISION,
                fiber_y  DOUBLE PRECISION,
                delta_x  DOUBLE PRECISION,
                delta_y  DOUBLE PRECISION,
                fiber_ra  DOUBLE PRECISION,
                fiber_dec  DOUBLE PRECISION,
                psf_to_fiber_specflux  DOUBLE PRECISION,
                run  TEXT,
                program  TEXT
                );
                """.format(tablename)
            
            
        with engine.connect() as conn:
            if overwrite:
                cmd = f'DROP TABLE IF EXISTS {tablename};'
                conn.execute(text(cmd))                            
            conn.execute(text(schema))
            conn.close()

    @staticmethod
    def fill_table(prod='everest',hdu=1,coadd='pernight',schema=False):
        tablename = ztile_everest.tablename(prod=prod,coadd=coadd,hdu=hdu,short=True)
            
        print(tablename)
        root = f"/global/project/projectdirs/desi/spectro/redux/{prod}/zcatalog/"
        runs=['main','sv1','sv2','sv3']
        programs = ['bright','dark','backup','other']
        
        dfs=[]
        for run in runs:
            for program in programs:
                filename = os.path.join(root,f"ztile-{run}-{program}-{coadd}.fits")
                print(filename)
                try:
                    dat = Table.read(filename, format='fits',hdu=hdu)
                except:
                    continue

                print(filename)

                # There is a multidimensional column that needs to be broken up
                if hdu==1:
                    for icoeff in range(0,10):
                        dat[f'COEFF_{icoeff}']= dat['COEFF'][0:len(dat),icoeff]
                    dat.remove_column('COEFF')

                dat.convert_bytestring_to_unicode()

                df = dat.to_pandas()
                df['run']=numpy.full(df.shape[0],run)
                df['program']=numpy.full(df.shape[0],program)
                df.columns= df.columns.str.lower()
#                     dtypesToSchema(df.dtypes)

                if schema:
                    dfs.append(df[0:1])

                if not schema:
                    try:
                        df.to_sql(tablename,engine,index=False,if_exists='append',schema=prod) 
                    except:
                        dtypesToSchema(df.dtypes)
                        print(df[0])
                        sys.exit()
        if schema:
            dfs = pandas.concat(dfs, ignore_index=True, sort=False)
            dtypesToSchema(dfs.dtypes)
            
            
class mtl(tablebaseclass):
    tablename = "mtl"
    createstatements = [ """
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
        )
        """
                        ]
    root = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/"
    runs = ["main","sv3","sv2"]
    runs=["sv3"]
    programs = ["ToO","","secondary"]
    lunations = ["bright","dark",""] 

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


      

class exposure_tables_daily(tablebaseclass):
    tablename = "exposure_tables_daily"
    createstatements = [ """
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
    )
    """
                         ]

    @staticmethod
    def fill_table():
        global engine
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
                data = ascii.read(file, fill_exclude_names=['BADCAMWORD', 'BADAMPS'])
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

class proposals_pv(tablebaseclass):
    tablename = "proposals_pv"
    createstatements = [ """
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
            lunation  TEXT,
            run       TEXT
    )
    """
                         ]
    
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
                    df['run']=numpy.full(df.shape[0],run)
                    df.columns= df.columns.str.lower()
                    try:
                        df.to_sql('proposals_pv',engine,index=False,if_exists='append')
                    except sqlite3.OperationalError as err:
                        dtypesToSchema(df.dtypes)
                        sys.exit()
        
    @staticmethod
    def fill_sga():
        global engine
        with engine.connect() as conn:
#             conn.execute("create table temp as (select * from proposals_pv);")
            query = """
            UPDATE proposals_pv
            SET sga_id = NULL
            WHERE sga_id <1;
            """
            conn.execute(query)
        
            query = """
            UPDATE proposals_pv
            SET sga_id = dr9_pv.sga_id
            FROM dr9_pv
            WHERE proposals_pv.objid = dr9_pv.objid
            AND proposals_pv.brickid = dr9_pv.brickid;
            """
            conn.execute(query)
  



                        
                        
class fibermap_daily(tablebaseclass):
    tablename = "fibermap_daily"
    createstatements = [ """
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
            hpxpixel  BIGINT,
            sv1_scnd_target  DOUBLE PRECISION,
            sv2_desi_target  DOUBLE PRECISION,
            sv2_bgs_target  DOUBLE PRECISION,
            sv2_mws_target  DOUBLE PRECISION,
            sv2_scnd_target  DOUBLE PRECISION,
            sv3_desi_target  DOUBLE PRECISION,
            sv3_bgs_target  DOUBLE PRECISION,
            sv3_mws_target  DOUBLE PRECISION,
            sv3_scnd_target  DOUBLE PRECISION,
            PRIMARY KEY (targetid, expid)
    )
    """
                         ]
  
    @staticmethod
    def fill_table(): 
        global engine
        magicdate='20210503'  # the last date that daily directory was created

#         dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
        
        
        
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
        for path in glob.glob(f'{dir_root}/cumulative/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
            
        dates = numpy.unique(dates)
        dates=dates[dates >= '20201214']
        dates=dates[dates != '20210213']    # there is something wrong with this date
        
        for date in dates:
            if int(date) not in dates_db:
                print(date)
                # Do things in terms of dates
                dfs=[]
                
                if date<=magicdate and date != '20210315':  # include dates that have missing daily reductions
                    use_root = dir_root
                else:
                    use_root = dir_root+'/cumulative/'
                    
                if date < '20210901':
                    base='zbest'
                else:
                    base='redrock'
                    
                for path in glob.glob(f'{use_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        print(date,tile)
                        for i in range(10):
                            if date<=magicdate and date != '20210315':
                                filename = f'{use_root}/{tile}/{date}/zbest-{i}-{tile}-{date}.fits'
                            else:
                                filename = f'{use_root}/{tile}/{date}/{base}-{i}-{tile}-thru{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits',hdu=3)
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

class mosthosts(tablebaseclass):
    """Maayane's collection of information for the mosthosts project.

    Added by rknop 2021-09-16

    This table is not very dynamic, so it's not part of regular loading;
    it should be run as a one-off.  The input is the csv file from
    Maayane Soumagnac.

    NOTE: For purposes of this loading, I edited the first line of the
    csv file.  I removed the leading #, and ai removed the spaces after
    the commas.

    """

    expected_columns = set( ['RA', 'DEC', 'PMRA', 'PMDEC', 'REF_EPOCH',
                             'OVERRIDE', 'snname', 'index', 'hemisphere',
                             'sn_ra', 'sn_dec', 'sn_z', 'program',
                             'priority', 'tns_name', 'iau_name', 'ptfiptf_name' ] )
    
    # hemisphere should be turned into an ENUM type

    tablename = "mosthosts.mosthosts"
    createstatements = ["""
      CREATE TABLE IF NOT EXISTS mosthosts.mosthosts (
          ra            DOUBLE PRECISION,
          dec           DOUBLE PRECISION,
          pmra          DOUBLE PRECISION,
          pmdec         DOUBLE PRECISION,
          ref_epoch     DOUBLE PRECISION,
          override      BOOLEAN,
          snname        TEXT,
          index         DOUBLE PRECISION,
          hemisphere    TEXT,
          sn_ra         DOUBLE PRECISION,
          sn_dec        DOUBLE PRECISION,
          sn_z          DOUBLE PRECISION,
          program       TEXT,
          priority      DOUBLE PRECISION,
          tns_name      TEXT,
          iau_name      TEXT,
          ptfiptf_name  TEXT
      )
    """,
              """
      CREATE INDEX mosthosts_q3c_idx ON mosthosts.mosthosts
         USING btree (public.q3c_ang2ipix(ra, dec));
    """,
              ]
              
              
    @staticmethod
    def fill_table( filename=None ):
        global engine
        if filename is None:
            raise ValueError( "Must pass a filename to mosthosts.fill_table" )
        sys.stderr.write( f'Reading {filename}\n' )
        tab = pandas.read_csv( filename )
        csvcolumns = set( tab.columns )
        if csvcolumns != mosthosts.expected_columns:
            raise ValueError( f'Columns in CSV file {csvcolumns} don\'t match what\'s expected.' )
        renames = {}
        for col in csvcolumns:
            if col != col.lower():
                renames[col] = col.lower()
        tab = tab.rename( renames, axis='columns' )
        sys.stderr.write( 'Loading database...\n' )
        tab.to_sql( 'mosthosts', con=engine, schema='mosthosts', index=False, if_exists="append" )
        sys.stderr.write( '...done.\n' )
        sys.stderr.write( 'Clustering q3c...\n' )
        with engine.connect() as conn:
            conn.execute( text( "CLUSTER mosthosts_q3c_idx ON mosthosts.mosthosts" ) )
            conn.execute( text( "ANALYZE mosthosts.mosthosts" ) )
        sys.stderr.write( '...done\n' )
        

class mosthosts(tablebaseclass):
    """Maayane's collection of information for the mosthosts project.
    Added by rknop 2021-09-16
    This table is not very dynamic, so it's not part of regular loading;
    it should be run as a one-off.  The input is the csv file from
    Maayane Soumagnac.
    NOTE: For purposes of this loading, I edited the first line of the
    csv file.  I removed the leading #, and ai removed the spaces after
    the commas.
    """

    expected_columns = set( ['RA', 'DEC', 'PMRA', 'PMDEC', 'REF_EPOCH',
                             'OVERRIDE', 'snname', 'index', 'hemisphere',
                             'sn_ra', 'sn_dec', 'sn_z', 'program',
                             'priority', 'tns_name', 'iau_name', 'ptfiptf_name' ] )
    
    # hemisphere should be turned into an ENUM type

    tablename = "mosthosts.mosthosts"
    createstatements = ["""
      CREATE TABLE IF NOT EXISTS mosthosts.mosthosts (
          ra            DOUBLE PRECISION,
          dec           DOUBLE PRECISION,
          pmra          DOUBLE PRECISION,
          pmdec         DOUBLE PRECISION,
          ref_epoch     DOUBLE PRECISION,
          override      BOOLEAN,
          snname        TEXT,
          index         DOUBLE PRECISION,
          hemisphere    TEXT,
          sn_ra         DOUBLE PRECISION,
          sn_dec        DOUBLE PRECISION,
          sn_z          DOUBLE PRECISION,
          program       TEXT,
          priority      DOUBLE PRECISION,
          tns_name      TEXT,
          iau_name      TEXT,
          ptfiptf_name  TEXT
      )
    """,
              """
      CREATE INDEX mosthosts_q3c_idx ON mosthosts.mosthosts
         USING btree (public.q3c_ang2ipix(ra, dec));
    """,
              ]
              
              
    @staticmethod
    def fill_table( filename=None ):
        global engine
        if filename is None:
            raise ValueError( "Must pass a filename to mosthosts.fill_table" )
        sys.stderr.write( f'Reading {filename}\n' )
        tab = pandas.read_csv( filename )
        csvcolumns = set( tab.columns )
        if csvcolumns != mosthosts.expected_columns:
            raise ValueError( f'Columns in CSV file {csvcolumns} don\'t match what\'s expected.' )
        renames = {}
        for col in csvcolumns:
            if col != col.lower():
                renames[col] = col.lower()
        tab = tab.rename( renames, axis='columns' )
        sys.stderr.write( 'Loading database...\n' )
        tab.to_sql( 'mosthosts', con=engine, schema='mosthosts', index=False, if_exists="append" )
        sys.stderr.write( '...done.\n' )
        sys.stderr.write( 'Clustering q3c...\n' )
        with engine.connect() as conn:
            conn.execute( text( "CLUSTER mosthosts_q3c_idx ON mosthosts.mosthosts" ) )
            conn.execute( text( "ANALYZE mosthosts.mosthosts" ) )
        sys.stderr.write( '...done\n' )

class fastspecphoto(tablebaseclass):
    tablename = "fastspecphoto"
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS fastspecphoto (  
        TARGETID BIGINT,
        RA DOUBLE PRECISION,
        DEC DOUBLE PRECISION,
        SURVEY TEXT,
        FAPRGRM TEXT,
        HPXPIXEL INTEGER,
        DESI_TARGET BIGINT,
        BGS_TARGET BIGINT,
        MWS_TARGET BIGINT,
        SV1_DESI_TARGET BIGINT,
        SV1_BGS_TARGET BIGINT,  
        SV1_MWS_TARGET BIGINT,
        SV2_DESI_TARGET BIGINT,
        SV2_BGS_TARGET BIGINT,  
        SV2_MWS_TARGET BIGINT,
        SV3_DESI_TARGET BIGINT,
        SV3_BGS_TARGET BIGINT,  
        SV3_MWS_TARGET BIGINT,
        SCND_TARGET BIGINT,
        SV1_SCND_TARGET BIGINT,  
        SV2_SCND_TARGET BIGINT,
        SV3_SCND_TARGET BIGINT,
        Z DOUBLE PRECISION,
        ZWARN BIGINT,
        DELTACHI2 DOUBLE PRECISION,
        SPECTYPE TEXT,
        PHOTSYS TEXT,  
        MW_TRANSMISSION_G REAL,
        MW_TRANSMISSION_R REAL,
        MW_TRANSMISSION_Z REAL,
        MW_TRANSMISSION_W1 REAL,
        MW_TRANSMISSION_W2 REAL,
        FIBERFLUX_G REAL,
        FIBERFLUX_R REAL,
        FIBERFLUX_Z REAL,
        FIBERTOTFLUX_G REAL,
        FIBERTOTFLUX_R REAL,
        FIBERTOTFLUX_Z REAL,    
        FLUX_G REAL,
        FLUX_R REAL,
        FLUX_Z REAL,
        FLUX_W1 REAL,
        FLUX_W2 REAL,    
        FLUX_IVAR_G REAL,
        FLUX_IVAR_R REAL,
        FLUX_IVAR_Z REAL,
        FLUX_IVAR_W1 REAL,
        FLUX_IVAR_W2 REAL,    
        CONTINUUM_Z DOUBLE PRECISION,
        CONTINUUM_CHI2 REAL,
        CONTINUUM_AGE REAL,
        CONTINUUM_AV REAL,
        CONTINUUM_AV_IVAR REAL,
        CONTINUUM_VDISP REAL,
        CONTINUUM_VDISP_IVAR REAL,
        CONTINUUM_SNR_B REAL,
        CONTINUUM_SNR_R REAL,
        CONTINUUM_SNR_Z REAL,
        CONTINUUM_SMOOTHCORR_B REAL,
        CONTINUUM_SMOOTHCORR_R REAL,
        CONTINUUM_SMOOTHCORR_Z REAL,
        DN4000 REAL,
        DN4000_IVAR REAL,
        DN4000_MODEL REAL,
        FLUX_SYNTH_G REAL,
        FLUX_SYNTH_R REAL,
        FLUX_SYNTH_Z REAL,
        FLUX_SYNTH_MODEL_G REAL,
        FLUX_SYNTH_MODEL_R REAL,
        FLUX_SYNTH_MODEL_Z REAL,
        BALMER_Z DOUBLE PRECISION,
        FORBIDDEN_Z DOUBLE PRECISION,
        BROAD_Z DOUBLE PRECISION,
        BALMER_SIGMA REAL,
        FORBIDDEN_SIGMA REAL,
        BROAD_SIGMA REAL,
        OI_1304_AMP REAL,
        OI_1304_AMP_IVAR REAL,
        OI_1304_FLUX REAL,
        OI_1304_FLUX_IVAR REAL,
        OI_1304_BOXFLUX REAL,
        OI_1304_VSHIFT REAL,
        OI_1304_SIGMA REAL,
        OI_1304_CONT REAL,
        OI_1304_CONT_IVAR REAL,
        OI_1304_EW REAL,
        OI_1304_EW_IVAR REAL,
        OI_1304_FLUX_LIMIT REAL,
        OI_1304_EW_LIMIT REAL,
        OI_1304_CHI2 REAL,
        OI_1304_NPIX INTEGER,
        SILIV_1396_AMP REAL,
        SILIV_1396_AMP_IVAR REAL,
        SILIV_1396_FLUX REAL,
        SILIV_1396_FLUX_IVAR REAL,
        SILIV_1396_BOXFLUX REAL,
        SILIV_1396_VSHIFT REAL,
        SILIV_1396_SIGMA REAL,
        SILIV_1396_CONT REAL,
        SILIV_1396_CONT_IVAR REAL,
        SILIV_1396_EW REAL,
        SILIV_1396_EW_IVAR REAL,
        SILIV_1396_FLUX_LIMIT REAL,
        SILIV_1396_EW_LIMIT REAL,
        SILIV_1396_CHI2 REAL,
        SILIV_1396_NPIX INTEGER,
        CIV_1549_AMP REAL,
        CIV_1549_AMP_IVAR REAL,
        CIV_1549_FLUX REAL,
        CIV_1549_FLUX_IVAR REAL,
        CIV_1549_BOXFLUX REAL,
        CIV_1549_VSHIFT REAL,
        CIV_1549_SIGMA REAL,
        CIV_1549_CONT REAL,
        CIV_1549_CONT_IVAR REAL,
        CIV_1549_EW REAL,
        CIV_1549_EW_IVAR REAL,
        CIV_1549_FLUX_LIMIT REAL,
        CIV_1549_EW_LIMIT REAL,
        CIV_1549_CHI2 REAL,
        CIV_1549_NPIX INTEGER,
        SILIII_1892_AMP REAL,
        SILIII_1892_AMP_IVAR REAL,
        SILIII_1892_FLUX REAL,
        SILIII_1892_FLUX_IVAR REAL,
        SILIII_1892_BOXFLUX REAL,
        SILIII_1892_VSHIFT REAL,
        SILIII_1892_SIGMA REAL,
        SILIII_1892_CONT REAL,
        SILIII_1892_CONT_IVAR REAL,
        SILIII_1892_EW REAL,
        SILIII_1892_EW_IVAR REAL,
        SILIII_1892_FLUX_LIMIT REAL,
        SILIII_1892_EW_LIMIT REAL,
        SILIII_1892_CHI2 REAL,
        SILIII_1892_NPIX INTEGER,
        CIII_1908_AMP REAL,
        CIII_1908_AMP_IVAR REAL,
        CIII_1908_FLUX REAL,
        CIII_1908_FLUX_IVAR REAL,
        CIII_1908_BOXFLUX REAL,
        CIII_1908_VSHIFT REAL,
        CIII_1908_SIGMA REAL,
        CIII_1908_CONT REAL,
        CIII_1908_CONT_IVAR REAL,
        CIII_1908_EW REAL,
        CIII_1908_EW_IVAR REAL,
        CIII_1908_FLUX_LIMIT REAL,
        CIII_1908_EW_LIMIT REAL,
        CIII_1908_CHI2 REAL,
        CIII_1908_NPIX INTEGER,
        MGII_2796_AMP REAL,
        MGII_2796_AMP_IVAR REAL,
        MGII_2796_FLUX REAL,
        MGII_2796_FLUX_IVAR REAL,
        MGII_2796_BOXFLUX REAL,
        MGII_2796_VSHIFT REAL,
        MGII_2796_SIGMA REAL,
        MGII_2796_CONT REAL,
        MGII_2796_CONT_IVAR REAL,
        MGII_2796_EW REAL,
        MGII_2796_EW_IVAR REAL,
        MGII_2796_FLUX_LIMIT REAL,
        MGII_2796_EW_LIMIT REAL,
        MGII_2796_CHI2 REAL,
        MGII_2796_NPIX INTEGER,
        MGII_2803_AMP REAL,
        MGII_2803_AMP_IVAR REAL,
        MGII_2803_FLUX REAL,
        MGII_2803_FLUX_IVAR REAL,
        MGII_2803_BOXFLUX REAL,
        MGII_2803_VSHIFT REAL,
        MGII_2803_SIGMA REAL,
        MGII_2803_CONT REAL,
        MGII_2803_CONT_IVAR REAL,
        MGII_2803_EW REAL,
        MGII_2803_EW_IVAR REAL,
        MGII_2803_FLUX_LIMIT REAL,
        MGII_2803_EW_LIMIT REAL,
        MGII_2803_CHI2 REAL,
        MGII_2803_NPIX INTEGER,
        NEV_3346_AMP REAL,
        NEV_3346_AMP_IVAR REAL,
        NEV_3346_FLUX REAL,
        NEV_3346_FLUX_IVAR REAL,
        NEV_3346_BOXFLUX REAL,
        NEV_3346_VSHIFT REAL,
        NEV_3346_SIGMA REAL,
        NEV_3346_CONT REAL,
        NEV_3346_CONT_IVAR REAL,
        NEV_3346_EW REAL,
        NEV_3346_EW_IVAR REAL,
        NEV_3346_FLUX_LIMIT REAL,
        NEV_3346_EW_LIMIT REAL,
        NEV_3346_CHI2 REAL,
        NEV_3346_NPIX INTEGER,
        NEV_3426_AMP REAL,
        NEV_3426_AMP_IVAR REAL,
        NEV_3426_FLUX REAL,
        NEV_3426_FLUX_IVAR REAL,
        NEV_3426_BOXFLUX REAL,
        NEV_3426_VSHIFT REAL,
        NEV_3426_SIGMA REAL,
        NEV_3426_CONT REAL,
        NEV_3426_CONT_IVAR REAL,
        NEV_3426_EW REAL,
        NEV_3426_EW_IVAR REAL,
        NEV_3426_FLUX_LIMIT REAL,
        NEV_3426_EW_LIMIT REAL,
        NEV_3426_CHI2 REAL,
        NEV_3426_NPIX INTEGER,
        OII_3726_AMP REAL,
        OII_3726_AMP_IVAR REAL,
        OII_3726_FLUX REAL,
        OII_3726_FLUX_IVAR REAL,
        OII_3726_BOXFLUX REAL,
        OII_3726_VSHIFT REAL,
        OII_3726_SIGMA REAL,
        OII_3726_CONT REAL,
        OII_3726_CONT_IVAR REAL,
        OII_3726_EW REAL,
        OII_3726_EW_IVAR REAL,
        OII_3726_FLUX_LIMIT REAL,
        OII_3726_EW_LIMIT REAL,
        OII_3726_CHI2 REAL,
        OII_3726_NPIX INTEGER,
        OII_3729_AMP REAL,
        OII_3729_AMP_IVAR REAL,
        OII_3729_FLUX REAL,
        OII_3729_FLUX_IVAR REAL,
        OII_3729_BOXFLUX REAL,
        OII_3729_VSHIFT REAL,
        OII_3729_SIGMA REAL,
        OII_3729_CONT REAL,
        OII_3729_CONT_IVAR REAL,
        OII_3729_EW REAL,
        OII_3729_EW_IVAR REAL,
        OII_3729_FLUX_LIMIT REAL,
        OII_3729_EW_LIMIT REAL,
        OII_3729_CHI2 REAL,
        OII_3729_NPIX INTEGER,
        NEIII_3869_AMP REAL,
        NEIII_3869_AMP_IVAR REAL,
        NEIII_3869_FLUX REAL,
        NEIII_3869_FLUX_IVAR REAL,
        NEIII_3869_BOXFLUX REAL,
        NEIII_3869_VSHIFT REAL,
        NEIII_3869_SIGMA REAL,
        NEIII_3869_CONT REAL,
        NEIII_3869_CONT_IVAR REAL,
        NEIII_3869_EW REAL,
        NEIII_3869_EW_IVAR REAL,
        NEIII_3869_FLUX_LIMIT REAL,
        NEIII_3869_EW_LIMIT REAL,
        NEIII_3869_CHI2 REAL,
        NEIII_3869_NPIX INTEGER,
        HEI_3889_AMP REAL,
        HEI_3889_AMP_IVAR REAL,
        HEI_3889_FLUX REAL,
        HEI_3889_FLUX_IVAR REAL,
        HEI_3889_BOXFLUX REAL,
        HEI_3889_VSHIFT REAL,
        HEI_3889_SIGMA REAL,
        HEI_3889_CONT REAL,
        HEI_3889_CONT_IVAR REAL,
        HEI_3889_EW REAL,
        HEI_3889_EW_IVAR REAL,
        HEI_3889_FLUX_LIMIT REAL,
        HEI_3889_EW_LIMIT REAL,
        HEI_3889_CHI2 REAL,
        HEI_3889_NPIX INTEGER,
        H6_AMP REAL,
        H6_AMP_IVAR REAL,
        H6_FLUX REAL,
        H6_FLUX_IVAR REAL,
        H6_BOXFLUX REAL,
        H6_VSHIFT REAL,
        H6_SIGMA REAL,
        H6_CONT REAL,
        H6_CONT_IVAR REAL,
        H6_EW REAL,
        H6_EW_IVAR REAL,
        H6_FLUX_LIMIT REAL,
        H6_EW_LIMIT REAL,
        H6_CHI2 REAL,
        H6_NPIX INTEGER,
        HEPSILON_AMP REAL,
        HEPSILON_AMP_IVAR REAL,
        HEPSILON_FLUX REAL,
        HEPSILON_FLUX_IVAR REAL,
        HEPSILON_BOXFLUX REAL,
        HEPSILON_VSHIFT REAL,
        HEPSILON_SIGMA REAL,
        HEPSILON_CONT REAL,
        HEPSILON_CONT_IVAR REAL,
        HEPSILON_EW REAL,
        HEPSILON_EW_IVAR REAL,
        HEPSILON_FLUX_LIMIT REAL,
        HEPSILON_EW_LIMIT REAL,
        HEPSILON_CHI2 REAL,
        HEPSILON_NPIX INTEGER,
        HDELTA_AMP REAL,
        HDELTA_AMP_IVAR REAL,
        HDELTA_FLUX REAL,
        HDELTA_FLUX_IVAR REAL,
        HDELTA_BOXFLUX REAL,
        HDELTA_VSHIFT REAL,
        HDELTA_SIGMA REAL,
        HDELTA_CONT REAL,
        HDELTA_CONT_IVAR REAL,
        HDELTA_EW REAL,
        HDELTA_EW_IVAR REAL,
        HDELTA_FLUX_LIMIT REAL,
        HDELTA_EW_LIMIT REAL,
        HDELTA_CHI2 REAL,
        HDELTA_NPIX INTEGER,
        HGAMMA_AMP REAL,
        HGAMMA_AMP_IVAR REAL,
        HGAMMA_FLUX REAL,
        HGAMMA_FLUX_IVAR REAL,
        HGAMMA_BOXFLUX REAL,
        HGAMMA_VSHIFT REAL,
        HGAMMA_SIGMA REAL,
        HGAMMA_CONT REAL,
        HGAMMA_CONT_IVAR REAL,
        HGAMMA_EW REAL,
        HGAMMA_EW_IVAR REAL,
        HGAMMA_FLUX_LIMIT REAL,
        HGAMMA_EW_LIMIT REAL,
        HGAMMA_CHI2 REAL,
        HGAMMA_NPIX INTEGER,
        OIII_4363_AMP REAL,
        OIII_4363_AMP_IVAR REAL,
        OIII_4363_FLUX REAL,
        OIII_4363_FLUX_IVAR REAL,
        OIII_4363_BOXFLUX REAL,
        OIII_4363_VSHIFT REAL,
        OIII_4363_SIGMA REAL,
        OIII_4363_CONT REAL,
        OIII_4363_CONT_IVAR REAL,
        OIII_4363_EW REAL,
        OIII_4363_EW_IVAR REAL,
        OIII_4363_FLUX_LIMIT REAL,
        OIII_4363_EW_LIMIT REAL,
        OIII_4363_CHI2 REAL,
        OIII_4363_NPIX INTEGER,
        HEI_4471_AMP REAL,
        HEI_4471_AMP_IVAR REAL,
        HEI_4471_FLUX REAL,
        HEI_4471_FLUX_IVAR REAL,
        HEI_4471_BOXFLUX REAL,
        HEI_4471_VSHIFT REAL,
        HEI_4471_SIGMA REAL,
        HEI_4471_CONT REAL,
        HEI_4471_CONT_IVAR REAL,
        HEI_4471_EW REAL,
        HEI_4471_EW_IVAR REAL,
        HEI_4471_FLUX_LIMIT REAL,
        HEI_4471_EW_LIMIT REAL,
        HEI_4471_CHI2 REAL,
        HEI_4471_NPIX INTEGER,
        HEII_4686_AMP REAL,
        HEII_4686_AMP_IVAR REAL,
        HEII_4686_FLUX REAL,
        HEII_4686_FLUX_IVAR REAL,
        HEII_4686_BOXFLUX REAL,
        HEII_4686_VSHIFT REAL,
        HEII_4686_SIGMA REAL,
        HEII_4686_CONT REAL,
        HEII_4686_CONT_IVAR REAL,
        HEII_4686_EW REAL,
        HEII_4686_EW_IVAR REAL,
        HEII_4686_FLUX_LIMIT REAL,
        HEII_4686_EW_LIMIT REAL,
        HEII_4686_CHI2 REAL,
        HEII_4686_NPIX INTEGER,
        HBETA_AMP REAL,
        HBETA_AMP_IVAR REAL,
        HBETA_FLUX REAL,
        HBETA_FLUX_IVAR REAL,
        HBETA_BOXFLUX REAL,
        HBETA_VSHIFT REAL,
        HBETA_SIGMA REAL,
        HBETA_CONT REAL,
        HBETA_CONT_IVAR REAL,
        HBETA_EW REAL,
        HBETA_EW_IVAR REAL,
        HBETA_FLUX_LIMIT REAL,
        HBETA_EW_LIMIT REAL,
        HBETA_CHI2 REAL,
        HBETA_NPIX INTEGER,
        OIII_4959_AMP REAL,
        OIII_4959_AMP_IVAR REAL,
        OIII_4959_FLUX REAL,
        OIII_4959_FLUX_IVAR REAL,
        OIII_4959_BOXFLUX REAL,
        OIII_4959_VSHIFT REAL,
        OIII_4959_SIGMA REAL,
        OIII_4959_CONT REAL,
        OIII_4959_CONT_IVAR REAL,
        OIII_4959_EW REAL,
        OIII_4959_EW_IVAR REAL,
        OIII_4959_FLUX_LIMIT REAL,
        OIII_4959_EW_LIMIT REAL,
        OIII_4959_CHI2 REAL,
        OIII_4959_NPIX INTEGER,
        OIII_5007_AMP REAL,
        OIII_5007_AMP_IVAR REAL,
        OIII_5007_FLUX REAL,
        OIII_5007_FLUX_IVAR REAL,
        OIII_5007_BOXFLUX REAL,
        OIII_5007_VSHIFT REAL,
        OIII_5007_SIGMA REAL,
        OIII_5007_CONT REAL,
        OIII_5007_CONT_IVAR REAL,
        OIII_5007_EW REAL,
        OIII_5007_EW_IVAR REAL,
        OIII_5007_FLUX_LIMIT REAL,
        OIII_5007_EW_LIMIT REAL,
        OIII_5007_CHI2 REAL,
        OIII_5007_NPIX INTEGER,
        NII_5755_AMP REAL,
        NII_5755_AMP_IVAR REAL,
        NII_5755_FLUX REAL,
        NII_5755_FLUX_IVAR REAL,
        NII_5755_BOXFLUX REAL,
        NII_5755_VSHIFT REAL,
        NII_5755_SIGMA REAL,
        NII_5755_CONT REAL,
        NII_5755_CONT_IVAR REAL,
        NII_5755_EW REAL,
        NII_5755_EW_IVAR REAL,
        NII_5755_FLUX_LIMIT REAL,
        NII_5755_EW_LIMIT REAL,
        NII_5755_CHI2 REAL,
        NII_5755_NPIX INTEGER,
        HEI_5876_AMP REAL,
        HEI_5876_AMP_IVAR REAL,
        HEI_5876_FLUX REAL,
        HEI_5876_FLUX_IVAR REAL,
        HEI_5876_BOXFLUX REAL,
        HEI_5876_VSHIFT REAL,
        HEI_5876_SIGMA REAL,
        HEI_5876_CONT REAL,
        HEI_5876_CONT_IVAR REAL,
        HEI_5876_EW REAL,
        HEI_5876_EW_IVAR REAL,
        HEI_5876_FLUX_LIMIT REAL,
        HEI_5876_EW_LIMIT REAL,
        HEI_5876_CHI2 REAL,
        HEI_5876_NPIX INTEGER,
        OI_6300_AMP REAL,
        OI_6300_AMP_IVAR REAL,
        OI_6300_FLUX REAL,
        OI_6300_FLUX_IVAR REAL,
        OI_6300_BOXFLUX REAL,
        OI_6300_VSHIFT REAL,
        OI_6300_SIGMA REAL,
        OI_6300_CONT REAL,
        OI_6300_CONT_IVAR REAL,
        OI_6300_EW REAL,
        OI_6300_EW_IVAR REAL,
        OI_6300_FLUX_LIMIT REAL,
        OI_6300_EW_LIMIT REAL,
        OI_6300_CHI2 REAL,
        OI_6300_NPIX INTEGER,
        NII_6548_AMP REAL,
        NII_6548_AMP_IVAR REAL,
        NII_6548_FLUX REAL,
        NII_6548_FLUX_IVAR REAL,
        NII_6548_BOXFLUX REAL,
        NII_6548_VSHIFT REAL,
        NII_6548_SIGMA REAL,
        NII_6548_CONT REAL,
        NII_6548_CONT_IVAR REAL,
        NII_6548_EW REAL,
        NII_6548_EW_IVAR REAL,
        NII_6548_FLUX_LIMIT REAL,
        NII_6548_EW_LIMIT REAL,
        NII_6548_CHI2 REAL,
        NII_6548_NPIX INTEGER,
        HALPHA_AMP REAL,
        HALPHA_AMP_IVAR REAL,
        HALPHA_FLUX REAL,
        HALPHA_FLUX_IVAR REAL,
        HALPHA_BOXFLUX REAL,
        HALPHA_VSHIFT REAL,
        HALPHA_SIGMA REAL,
        HALPHA_CONT REAL,
        HALPHA_CONT_IVAR REAL,
        HALPHA_EW REAL,
        HALPHA_EW_IVAR REAL,
        HALPHA_FLUX_LIMIT REAL,
        HALPHA_EW_LIMIT REAL,
        HALPHA_CHI2 REAL,
        HALPHA_NPIX INTEGER,
        NII_6584_AMP REAL,
        NII_6584_AMP_IVAR REAL,
        NII_6584_FLUX REAL,
        NII_6584_FLUX_IVAR REAL,
        NII_6584_BOXFLUX REAL,
        NII_6584_VSHIFT REAL,
        NII_6584_SIGMA REAL,
        NII_6584_CONT REAL,
        NII_6584_CONT_IVAR REAL,
        NII_6584_EW REAL,
        NII_6584_EW_IVAR REAL,
        NII_6584_FLUX_LIMIT REAL,
        NII_6584_EW_LIMIT REAL,
        NII_6584_CHI2 REAL,
        NII_6584_NPIX INTEGER,
        SII_6716_AMP REAL,
        SII_6716_AMP_IVAR REAL,
        SII_6716_FLUX REAL,
        SII_6716_FLUX_IVAR REAL,
        SII_6716_BOXFLUX REAL,
        SII_6716_VSHIFT REAL,
        SII_6716_SIGMA REAL,
        SII_6716_CONT REAL,
        SII_6716_CONT_IVAR REAL,
        SII_6716_EW REAL,
        SII_6716_EW_IVAR REAL,
        SII_6716_FLUX_LIMIT REAL,
        SII_6716_EW_LIMIT REAL,
        SII_6716_CHI2 REAL,
        SII_6716_NPIX INTEGER,
        SII_6731_AMP REAL,
        SII_6731_AMP_IVAR REAL,
        SII_6731_FLUX REAL,
        SII_6731_FLUX_IVAR REAL,
        SII_6731_BOXFLUX REAL,
        SII_6731_VSHIFT REAL,
        SII_6731_SIGMA REAL,
        SII_6731_CONT REAL,
        SII_6731_CONT_IVAR REAL,
        SII_6731_EW REAL,
        SII_6731_EW_IVAR REAL,
        SII_6731_FLUX_LIMIT REAL,
        SII_6731_EW_LIMIT REAL,
        SII_6731_CHI2 REAL,
        SII_6731_NPIX INTEGER,
        OII_7320_AMP REAL,
        OII_7320_AMP_IVAR REAL,
        OII_7320_FLUX REAL,
        OII_7320_FLUX_IVAR REAL,
        OII_7320_BOXFLUX REAL,
        OII_7320_VSHIFT REAL,
        OII_7320_SIGMA REAL,
        OII_7320_CONT REAL,
        OII_7320_CONT_IVAR REAL,
        OII_7320_EW REAL,
        OII_7320_EW_IVAR REAL,
        OII_7320_FLUX_LIMIT REAL,
        OII_7320_EW_LIMIT REAL,
        OII_7320_CHI2 REAL,
        OII_7320_NPIX INTEGER,
        OII_7330_AMP REAL,
        OII_7330_AMP_IVAR REAL,
        OII_7330_FLUX REAL,
        OII_7330_FLUX_IVAR REAL,
        OII_7330_BOXFLUX REAL,
        OII_7330_VSHIFT REAL,
        OII_7330_SIGMA REAL,
        OII_7330_CONT REAL,
        OII_7330_CONT_IVAR REAL,
        OII_7330_EW REAL,
        OII_7330_EW_IVAR REAL,
        OII_7330_FLUX_LIMIT REAL,
        OII_7330_EW_LIMIT REAL,
        OII_7330_CHI2 REAL,
        OII_7330_NPIX INTEGER,
        SIII_9069_AMP REAL,
        SIII_9069_AMP_IVAR REAL,
        SIII_9069_FLUX REAL,
        SIII_9069_FLUX_IVAR REAL,
        SIII_9069_BOXFLUX REAL,
        SIII_9069_VSHIFT REAL,
        SIII_9069_SIGMA REAL,
        SIII_9069_CONT REAL,
        SIII_9069_CONT_IVAR REAL,
        SIII_9069_EW REAL,
        SIII_9069_EW_IVAR REAL,
        SIII_9069_FLUX_LIMIT REAL,
        SIII_9069_EW_LIMIT REAL,
        SIII_9069_CHI2 REAL,
        SIII_9069_NPIX INTEGER,
        SIII_9532_AMP REAL,
        SIII_9532_AMP_IVAR REAL,
        SIII_9532_FLUX REAL,
        SIII_9532_FLUX_IVAR REAL,
        SIII_9532_BOXFLUX REAL,
        SIII_9532_VSHIFT REAL,
        SIII_9532_SIGMA REAL,
        SIII_9532_CONT REAL,
        SIII_9532_CONT_IVAR REAL,
        SIII_9532_EW REAL,
        SIII_9532_EW_IVAR REAL,
        SIII_9532_FLUX_LIMIT REAL,
        SIII_9532_EW_LIMIT REAL,
        SIII_9532_CHI2 REAL,
        SIII_9532_NPIX INTEGER
        CONTINUUM_CHI2 REAL,
        CONTINUUM_AGE REAL,
        CONTINUUM_AV REAL,
        CONTINUUM_AV_IVAR REAL,
        DN4000_MODEL REAL,
        FLUX_SYNTH_MODEL_G REAL,
        FLUX_SYNTH_MODEL_R REAL,
        FLUX_SYNTH_MODEL_Z REAL,
        FLUX_SYNTH_MODEL_W1 REAL,
        FLUX_SYNTH_MODEL_W2 REAL,
        KCORR_U REAL,
        ABSMAG_U REAL,
        ABSMAG_IVAR_U REAL,
        KCORR_B REAL,
        ABSMAG_B REAL,
        ABSMAG_IVAR_B REAL,
        KCORR_V REAL,
        ABSMAG_V REAL,
        ABSMAG_IVAR_V REAL,
        KCORR_SDSS_U REAL,
        ABSMAG_SDSS_U REAL,
        ABSMAG_IVAR_SDSS_U REAL,
        KCORR_SDSS_G REAL,
        ABSMAG_SDSS_G REAL,
        ABSMAG_IVAR_SDSS_G REAL,
        KCORR_SDSS_R REAL,
        ABSMAG_SDSS_R REAL,
        ABSMAG_IVAR_SDSS_R REAL,
        KCORR_SDSS_I REAL,
        ABSMAG_SDSS_I REAL,
        ABSMAG_IVAR_SDSS_I REAL,
        KCORR_SDSS_Z REAL,
        ABSMAG_SDSS_Z REAL,
        ABSMAG_IVAR_SDSS_Z REAL,
        KCORR_W1 REAL,
        ABSMAG_W1 REAL,
        ABSMAG_IVAR_W1 REAL,
        target  TEXT,
        survey  TEXT
    );
    """
                        ]

    @staticmethod
    def fill_table():   
        #dirs=["https://data.desi.lbl.gov/desi/spectro/fastspecfit/everest/catalogs/"]
        surveys=["sv3","main"]
        targets = ["bright","dark"]
        for survey in surveys:
            for target in targets:
                
                filename_spec = f"/cfs/cdirs/desi/spectro/fastspecfit/everest/catalogs/fastspec-everest-{survey}-{target}.fits"
                filename_photo = f"/cfs/cdirs/desi/spectro/fastspecfit/everest/catalogs/fastphot-everest-{survey}-{target}.fits"
                try:
                    dat1 = Table.read(filename_spec, format='fits',hdu=1)
                    names1 = [name for name in dat1.colnames if len(dat1[name].shape) <= 1]
                    names1=names1[5:]
                    dat1=dat1[names1]
                    
                    dat0 = Table.read(filename_spec, format='fits',hdu=2)
                    names0 = [name for name in dat0.colnames if len(dat0[name].shape) <= 1]
                    dat0=dat0[names0]

                    dat2 = Table.read(filename_photo, format='fits',hdu=1)
                    names2 = [name for name in dat2.colnames if len(dat2[name].shape) <= 1]
                    names2=names2[1:]
                    dat2=dat2[names2]
            
                    
                    dat=hstack([dat0, dat1, dat2],metadata_conflicts='silent')
                except:
                    print(f"{filename} not found")
                    continue
                dat.convert_bytestring_to_unicode()
                df = dat.to_pandas()
                df['target']=numpy.full(df.shape[0],target)
                df['survey']=numpy.full(df.shape[0],survey)
                df.columns= df.columns.str.lower()
                dtypesToSchema(df.dtypes)
                wef
                try:
                    df.to_sql('fastspecphoto_data',engine,index=False,if_exists='append')
                except:
                    dtypesToSchema(df.dtypes)
                    sys.exit()           

# ======================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="DESI database loader" )
    parser.add_argument( "-s", "--secrets-file", default="/global/homes/a/akim/secrets/desi_pg.txt",
                         help=( "File with \"username password\" for database "
                                "(default: /global/homes/a/akim/secrets/desi_pg.txt" ) )
    # parser.add_argument( "-a", "--action", default="daily-fill",
    #                      help=( "What to do (default: daily-fill) "
    #                             "Possibilities: daily-fill load-mosthosts create-all-missing-tables" ) )
    parser.add_argument( "-a", "--action", default=None,
                         help=( "What to do (default: None) "
                                "Possibilities: daily-fill load-mosthosts create-all-missing-tables" ) )
    parser.add_argument( "--really-do", default=False, action='store_true',
                         help=argparse.SUPPRESS )
    parser.add_argument( "-m", "--mosthosts-file", default="DESI_transients_w_tns.txt",
                         help="Name of CSV file for use with load-mosthosts (default: DESI_transients_w_tns.txt)" )
    parser.add_argument( "--dbhost", default="decatdb.lbl.gov",
                         help="Host machine of database (default: decatdb.lbl.gov)" )
    parser.add_argument( "--db", default="desi",
                         help="Name of database (default: desi)" )
    
    args = parser.parse_args()

    create_db_engine( args.secrets_file, dbhost=args.dbhost, dbname=args.db )

    if args.action == "daily-fill":
        # running as a cron job on cori10
        zbest_daily.fill_table()
        fibermap_daily.fill_table()
        exposure_tables_daily.fill_table()

    elif args.action == "load-mosthosts":
        if not args.really_do:
            print( "Running this will delete the current mosthosts.mosthosts\n"
                   "table and make a new one from the specified file.\n"
                   "If this is in fact what you want to do, rerun with --really-do" )
        else:
            mosthosts.create_table( overwrite=True )
            mosthosts.fill_table( args.mosthosts_file )
    elif args.action == "load-fastspecphoto":
        if not args.really_do:
            print( "Running this will delete the current fastspecphoto_data\n"
                   "table and make a new one from the specified file.\n"
                   "If this is in fact what you want to do, rerun with --really-do" )
        else:
#             fastspecphoto.create_table(overwrite=True )
            fastspecphoto.fill_table( )
    elif args.action == "create-all-missing-tables":
        # TODO : a better place to list all the classes with tables
        #  that need to be called here
        zbest_daily.create_table()
        dr9_pv.create_table()
        secondary.create_table()
        zcatalog_denali.create_table()
        everest_zpix_redshifts.create_table()
        everest_zpix_fibermap.create_table()
        # TODO : different values of hdu, coadd for ztile_everest
        ztile_everest.create_table()
        mtl.create_table()
        exposure_tables_daily.create_table()
        proposals_pv.create_table()
        fibermap_daily.create_table()
        mosthosts.create_table()
        fastspecphoto.create_table()
        
        
    else:
        raise ValueError( f'Unknown action {args.action}' )

else:
    # I put this here for backwards compatibility just in case Alex runs
    # this as a module inside other python scripts he has -- RKNOP
    # 2021-09-16
    create_db_engine( "/global/homes/a/akim/secrets/desi_pg.txt" )
    
