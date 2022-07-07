import os
import sys
import glob
import argparse
from astropy.io import ascii
import pandas
import psycopg2
import re
from astropy.table import Table, hstack
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
                        
                        """
        CREATE INDEX IF NOT EXISTS zpix_redshifts_q3c_ang2ipix_mean_fiber_idx
             ON everest.zpix_redshifts USING btree(public.q3c_ang2ipix(mean_fiber_ra, mean_fiber_dec))
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
    """,
                        """
        CREATE INDEX IF NOT EXISTS fibermap_daily_q3c_ang2ipix_idx
            ON fibermap_daily USING btree(public.q3c_ang2ipix(fiber_ra, fiber_dec))
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

class everest_daily_fibermap(tablebaseclass):
    tablename = "everest.daily_fibermap"
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS everest.daily_fibermap (
        targetid  BIGINT,
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
        release  SMALLINT,
        brickname  TEXT,
        brickid  INTEGER,
        brick_objid  INTEGER,
        morphtype  TEXT,
        ebv  REAL,
        flux_g  REAL,
        flux_r  REAL,
        flux_z  REAL,
        flux_w1  REAL,
        flux_w2  REAL,
        flux_ivar_g  REAL,
        flux_ivar_r  REAL,
        flux_ivar_z  REAL,
        flux_ivar_w1  REAL,
        flux_ivar_w2  REAL,
        fiberflux_g  REAL,
        fiberflux_r  REAL,
        fiberflux_z  REAL,
        fibertotflux_g  REAL,
        fibertotflux_r  REAL,
        fibertotflux_z  REAL,
        maskbits  SMALLINT,
        sersic  REAL,
        shape_r  REAL,
        shape_e1  REAL,
        shape_e2  REAL,
        ref_id  BIGINT,
        ref_cat  TEXT,
        gaia_phot_g_mean_mag  REAL,
        gaia_phot_bp_mean_mag  REAL,
        gaia_phot_rp_mean_mag  REAL,
        parallax  REAL,
        photsys  TEXT,
        priority_init  BIGINT,
        numobs_init  BIGINT,
        desi_target  BIGINT,
        bgs_target  BIGINT,
        mws_target  BIGINT,
        scnd_target  BIGINT,
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
        yyyymmdd  BIGINT
        PRIMARY KEY (targetid, yyyymmdd)
    )
    """,
                        """
        CREATE INDEX IF NOT EXISTS everest.daily_fibermap_q3c_ang2ipix_idx
            ON everst.daily_fibermap USING btree(public.q3c_ang2ipix(fiber_ra, fiber_dec))
    """
                         ]
  
    @staticmethod
    def fill_table(): 
        global engine
        firstdate='20210901'  # the last date that daily directory was created

        dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"
        # dir_root = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/"
        
        
        
        #find the last date
        dates_db=[]
        try:
            with engine.connect() as conn:
                for row in conn.execute(f"SELECT DISTINCT YYYYMMDD FROM everest.daily_fibermap"):
                    dates_db.append(row[0])
        except:
            pass
            
        dates = []
        for path in glob.glob(f'{dir_root}/*/202?????'):
            split = path.split('/')
            dates.append(split[-1])
            
        dates = numpy.unique(dates)
        dates=dates[dates >= firstdate]
        
        for date in dates:
            if int(date) not in dates_db:
                print(date)
                # Do things in terms of dates
                dfs=[]
                base='redrock'
                    
                for path in glob.glob(f'{dir_root}/*/{date}'):
                    split = path.split('/')
                    tile = split[-2]
                    if tile.isnumeric():
                        print(date,tile)
                        for i in range(10):
                            filename = f'{dir_root}/{tile}/{date}/{base}-{i}-{tile}-thru{date}.fits'
                            try:
                                dat = Table.read(filename, format='fits',hdu=2)
                            except:
                                print(f"{filename} not found")
                                continue

                            dat.convert_bytestring_to_unicode()
                            df = dat.to_pandas()
                            df['YYYYMMDD']=numpy.full(df.shape[0],int(date))
                            
                            df.columns= df.columns.str.lower()
                            dfs.append(df)
                if len(dfs)>0:
                    try:
                        dfs = pandas.concat(dfs, ignore_index=True, sort=False)
                        dfs.to_sql(f'daily_fibermap',engine,index=False,if_exists='append',schema='everest')
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

class fastspecfit_meta(tablebaseclass):
    tablename = "everest.fastspecfit_meta"
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS everest.fastspecfit_meta (  
        targetid  BIGINT,
        ra  DOUBLE PRECISION,
        dec  DOUBLE PRECISION,
        survey  TEXT,
        faprgrm  TEXT,
        hpxpixel  INTEGER,
        desi_target  BIGINT,
        bgs_target  BIGINT,
        mws_target  BIGINT,
        sv1_desi_target  BIGINT,
        sv1_bgs_target  BIGINT,
        sv1_mws_target  BIGINT,
        sv2_desi_target  BIGINT,
        sv2_bgs_target  BIGINT,
        sv2_mws_target  BIGINT,
        sv3_desi_target  BIGINT,
        sv3_bgs_target  BIGINT,
        sv3_mws_target  BIGINT,
        scnd_target  BIGINT,
        sv1_scnd_target  BIGINT,
        sv2_scnd_target  BIGINT,
        sv3_scnd_target  BIGINT,
        z  DOUBLE PRECISION,
        zwarn  BIGINT,
        deltachi2  DOUBLE PRECISION,
        spectype  TEXT,
        photsys  TEXT,
        mw_transmission_g  REAL,
        mw_transmission_r  REAL,
        mw_transmission_z  REAL,
        mw_transmission_w1  REAL,
        mw_transmission_w2  REAL,
        fiberflux_g  REAL,
        fiberflux_r  REAL,
        fiberflux_z  REAL,
        fibertotflux_g  REAL,
        fibertotflux_r  REAL,
        fibertotflux_z  REAL,
        flux_g  REAL,
        flux_r  REAL,
        flux_z  REAL,
        flux_w1  REAL,
        flux_w2  REAL,
        flux_ivar_g  REAL,
        flux_ivar_r  REAL,
        flux_ivar_z  REAL,
        flux_ivar_w1  REAL,
        flux_ivar_w2  REAL,
        PRIMARY KEY (targetid, survey)
        );
    """]


    @staticmethod
    def fill_table():   
        #dirs=["https://data.desi.lbl.gov/desi/spectro/fastspecfit/everest/catalogs/"]
        surveys=["sv3","main"]
        targets = ["bright","dark"]
        for survey in surveys:
            for target in targets:
                print(survey, target)
                filename_spec = f"/global/cfs/cdirs/desi/spectro/fastspecfit/everest/catalogs/fastspec-everest-{survey}-{target}.fits"
                try:              
                    dat = Table.read(filename_spec, format='fits',hdu=2)
                except:
                    print(f"{filename} not found")
                    continue
                dat.convert_bytestring_to_unicode()
                df = dat.to_pandas()
                df.columns= df.columns.str.lower()
                try:
                    df.to_sql('fastspecfit_meta',engine,index=False,if_exists='append',schema='everest')
                except:
                    dtypesToSchema(df.dtypes)
                    sys.exit()

                    
class fastspecfit_spec(tablebaseclass):
    tablename = "everest.fastspecfit_spec"
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS everest.fastspecfit_spec ( 
        targetid  BIGINT,
        continuum_z  DOUBLE PRECISION,
        continuum_chi2  REAL,
        continuum_age  REAL,
        continuum_av  REAL,
        continuum_av_ivar  REAL,
        continuum_vdisp  REAL,
        continuum_vdisp_ivar  REAL,
        continuum_snr_b  REAL,
        continuum_snr_r  REAL,
        continuum_snr_z  REAL,
        continuum_smoothcorr_b  REAL,
        continuum_smoothcorr_r  REAL,
        continuum_smoothcorr_z  REAL,
        dn4000  REAL,
        dn4000_ivar  REAL,
        dn4000_model  REAL,
        flux_synth_g  REAL,
        flux_synth_r  REAL,
        flux_synth_z  REAL,
        flux_synth_model_g  REAL,
        flux_synth_model_r  REAL,
        flux_synth_model_z  REAL,
        balmer_z  DOUBLE PRECISION,
        forbidden_z  DOUBLE PRECISION,
        broad_z  DOUBLE PRECISION,
        balmer_sigma  REAL,
        forbidden_sigma  REAL,
        broad_sigma  REAL,
        oi_1304_amp  REAL,
        oi_1304_amp_ivar  REAL,
        oi_1304_flux  REAL,
        oi_1304_flux_ivar  REAL,
        oi_1304_boxflux  REAL,
        oi_1304_vshift  REAL,
        oi_1304_sigma  REAL,
        oi_1304_cont  REAL,
        oi_1304_cont_ivar  REAL,
        oi_1304_ew  REAL,
        oi_1304_ew_ivar  REAL,
        oi_1304_flux_limit  REAL,
        oi_1304_ew_limit  REAL,
        oi_1304_chi2  REAL,
        oi_1304_npix  INTEGER,
        siliv_1396_amp  REAL,
        siliv_1396_amp_ivar  REAL,
        siliv_1396_flux  REAL,
        siliv_1396_flux_ivar  REAL,
        siliv_1396_boxflux  REAL,
        siliv_1396_vshift  REAL,
        siliv_1396_sigma  REAL,
        siliv_1396_cont  REAL,
        siliv_1396_cont_ivar  REAL,
        siliv_1396_ew  REAL,
        siliv_1396_ew_ivar  REAL,
        siliv_1396_flux_limit  REAL,
        siliv_1396_ew_limit  REAL,
        siliv_1396_chi2  REAL,
        siliv_1396_npix  INTEGER,
        civ_1549_amp  REAL,
        civ_1549_amp_ivar  REAL,
        civ_1549_flux  REAL,
        civ_1549_flux_ivar  REAL,
        civ_1549_boxflux  REAL,
        civ_1549_vshift  REAL,
        civ_1549_sigma  REAL,
        civ_1549_cont  REAL,
        civ_1549_cont_ivar  REAL,
        civ_1549_ew  REAL,
        civ_1549_ew_ivar  REAL,
        civ_1549_flux_limit  REAL,
        civ_1549_ew_limit  REAL,
        civ_1549_chi2  REAL,
        civ_1549_npix  INTEGER,
        siliii_1892_amp  REAL,
        siliii_1892_amp_ivar  REAL,
        siliii_1892_flux  REAL,
        siliii_1892_flux_ivar  REAL,
        siliii_1892_boxflux  REAL,
        siliii_1892_vshift  REAL,
        siliii_1892_sigma  REAL,
        siliii_1892_cont  REAL,
        siliii_1892_cont_ivar  REAL,
        siliii_1892_ew  REAL,
        siliii_1892_ew_ivar  REAL,
        siliii_1892_flux_limit  REAL,
        siliii_1892_ew_limit  REAL,
        siliii_1892_chi2  REAL,
        siliii_1892_npix  INTEGER,
        ciii_1908_amp  REAL,
        ciii_1908_amp_ivar  REAL,
        ciii_1908_flux  REAL,
        ciii_1908_flux_ivar  REAL,
        ciii_1908_boxflux  REAL,
        ciii_1908_vshift  REAL,
        ciii_1908_sigma  REAL,
        ciii_1908_cont  REAL,
        ciii_1908_cont_ivar  REAL,
        ciii_1908_ew  REAL,
        ciii_1908_ew_ivar  REAL,
        ciii_1908_flux_limit  REAL,
        ciii_1908_ew_limit  REAL,
        ciii_1908_chi2  REAL,
        ciii_1908_npix  INTEGER,
        mgii_2796_amp  REAL,
        mgii_2796_amp_ivar  REAL,
        mgii_2796_flux  REAL,
        mgii_2796_flux_ivar  REAL,
        mgii_2796_boxflux  REAL,
        mgii_2796_vshift  REAL,
        mgii_2796_sigma  REAL,
        mgii_2796_cont  REAL,
        mgii_2796_cont_ivar  REAL,
        mgii_2796_ew  REAL,
        mgii_2796_ew_ivar  REAL,
        mgii_2796_flux_limit  REAL,
        mgii_2796_ew_limit  REAL,
        mgii_2796_chi2  REAL,
        mgii_2796_npix  INTEGER,
        mgii_2803_amp  REAL,
        mgii_2803_amp_ivar  REAL,
        mgii_2803_flux  REAL,
        mgii_2803_flux_ivar  REAL,
        mgii_2803_boxflux  REAL,
        mgii_2803_vshift  REAL,
        mgii_2803_sigma  REAL,
        mgii_2803_cont  REAL,
        mgii_2803_cont_ivar  REAL,
        mgii_2803_ew  REAL,
        mgii_2803_ew_ivar  REAL,
        mgii_2803_flux_limit  REAL,
        mgii_2803_ew_limit  REAL,
        mgii_2803_chi2  REAL,
        mgii_2803_npix  INTEGER,
        nev_3346_amp  REAL,
        nev_3346_amp_ivar  REAL,
        nev_3346_flux  REAL,
        nev_3346_flux_ivar  REAL,
        nev_3346_boxflux  REAL,
        nev_3346_vshift  REAL,
        nev_3346_sigma  REAL,
        nev_3346_cont  REAL,
        nev_3346_cont_ivar  REAL,
        nev_3346_ew  REAL,
        nev_3346_ew_ivar  REAL,
        nev_3346_flux_limit  REAL,
        nev_3346_ew_limit  REAL,
        nev_3346_chi2  REAL,
        nev_3346_npix  INTEGER,
        nev_3426_amp  REAL,
        nev_3426_amp_ivar  REAL,
        nev_3426_flux  REAL,
        nev_3426_flux_ivar  REAL,
        nev_3426_boxflux  REAL,
        nev_3426_vshift  REAL,
        nev_3426_sigma  REAL,
        nev_3426_cont  REAL,
        nev_3426_cont_ivar  REAL,
        nev_3426_ew  REAL,
        nev_3426_ew_ivar  REAL,
        nev_3426_flux_limit  REAL,
        nev_3426_ew_limit  REAL,
        nev_3426_chi2  REAL,
        nev_3426_npix  INTEGER,
        oii_3726_amp  REAL,
        oii_3726_amp_ivar  REAL,
        oii_3726_flux  REAL,
        oii_3726_flux_ivar  REAL,
        oii_3726_boxflux  REAL,
        oii_3726_vshift  REAL,
        oii_3726_sigma  REAL,
        oii_3726_cont  REAL,
        oii_3726_cont_ivar  REAL,
        oii_3726_ew  REAL,
        oii_3726_ew_ivar  REAL,
        oii_3726_flux_limit  REAL,
        oii_3726_ew_limit  REAL,
        oii_3726_chi2  REAL,
        oii_3726_npix  INTEGER,
        oii_3729_amp  REAL,
        oii_3729_amp_ivar  REAL,
        oii_3729_flux  REAL,
        oii_3729_flux_ivar  REAL,
        oii_3729_boxflux  REAL,
        oii_3729_vshift  REAL,
        oii_3729_sigma  REAL,
        oii_3729_cont  REAL,
        oii_3729_cont_ivar  REAL,
        oii_3729_ew  REAL,
        oii_3729_ew_ivar  REAL,
        oii_3729_flux_limit  REAL,
        oii_3729_ew_limit  REAL,
        oii_3729_chi2  REAL,
        oii_3729_npix  INTEGER,
        neiii_3869_amp  REAL,
        neiii_3869_amp_ivar  REAL,
        neiii_3869_flux  REAL,
        neiii_3869_flux_ivar  REAL,
        neiii_3869_boxflux  REAL,
        neiii_3869_vshift  REAL,
        neiii_3869_sigma  REAL,
        neiii_3869_cont  REAL,
        neiii_3869_cont_ivar  REAL,
        neiii_3869_ew  REAL,
        neiii_3869_ew_ivar  REAL,
        neiii_3869_flux_limit  REAL,
        neiii_3869_ew_limit  REAL,
        neiii_3869_chi2  REAL,
        neiii_3869_npix  INTEGER,
        hei_3889_amp  REAL,
        hei_3889_amp_ivar  REAL,
        hei_3889_flux  REAL,
        hei_3889_flux_ivar  REAL,
        hei_3889_boxflux  REAL,
        hei_3889_vshift  REAL,
        hei_3889_sigma  REAL,
        hei_3889_cont  REAL,
        hei_3889_cont_ivar  REAL,
        hei_3889_ew  REAL,
        hei_3889_ew_ivar  REAL,
        hei_3889_flux_limit  REAL,
        hei_3889_ew_limit  REAL,
        hei_3889_chi2  REAL,
        hei_3889_npix  INTEGER,
        h6_amp  REAL,
        h6_amp_ivar  REAL,
        h6_flux  REAL,
        h6_flux_ivar  REAL,
        h6_boxflux  REAL,
        h6_vshift  REAL,
        h6_sigma  REAL,
        h6_cont  REAL,
        h6_cont_ivar  REAL,
        h6_ew  REAL,
        h6_ew_ivar  REAL,
        h6_flux_limit  REAL,
        h6_ew_limit  REAL,
        h6_chi2  REAL,
        h6_npix  INTEGER,
        hepsilon_amp  REAL,
        hepsilon_amp_ivar  REAL,
        hepsilon_flux  REAL,
        hepsilon_flux_ivar  REAL,
        hepsilon_boxflux  REAL,
        hepsilon_vshift  REAL,
        hepsilon_sigma  REAL,
        hepsilon_cont  REAL,
        hepsilon_cont_ivar  REAL,
        hepsilon_ew  REAL,
        hepsilon_ew_ivar  REAL,
        hepsilon_flux_limit  REAL,
        hepsilon_ew_limit  REAL,
        hepsilon_chi2  REAL,
        hepsilon_npix  INTEGER,
        hdelta_amp  REAL,
        hdelta_amp_ivar  REAL,
        hdelta_flux  REAL,
        hdelta_flux_ivar  REAL,
        hdelta_boxflux  REAL,
        hdelta_vshift  REAL,
        hdelta_sigma  REAL,
        hdelta_cont  REAL,
        hdelta_cont_ivar  REAL,
        hdelta_ew  REAL,
        hdelta_ew_ivar  REAL,
        hdelta_flux_limit  REAL,
        hdelta_ew_limit  REAL,
        hdelta_chi2  REAL,
        hdelta_npix  INTEGER,
        hgamma_amp  REAL,
        hgamma_amp_ivar  REAL,
        hgamma_flux  REAL,
        hgamma_flux_ivar  REAL,
        hgamma_boxflux  REAL,
        hgamma_vshift  REAL,
        hgamma_sigma  REAL,
        hgamma_cont  REAL,
        hgamma_cont_ivar  REAL,
        hgamma_ew  REAL,
        hgamma_ew_ivar  REAL,
        hgamma_flux_limit  REAL,
        hgamma_ew_limit  REAL,
        hgamma_chi2  REAL,
        hgamma_npix  INTEGER,
        oiii_4363_amp  REAL,
        oiii_4363_amp_ivar  REAL,
        oiii_4363_flux  REAL,
        oiii_4363_flux_ivar  REAL,
        oiii_4363_boxflux  REAL,
        oiii_4363_vshift  REAL,
        oiii_4363_sigma  REAL,
        oiii_4363_cont  REAL,
        oiii_4363_cont_ivar  REAL,
        oiii_4363_ew  REAL,
        oiii_4363_ew_ivar  REAL,
        oiii_4363_flux_limit  REAL,
        oiii_4363_ew_limit  REAL,
        oiii_4363_chi2  REAL,
        oiii_4363_npix  INTEGER,
        hei_4471_amp  REAL,
        hei_4471_amp_ivar  REAL,
        hei_4471_flux  REAL,
        hei_4471_flux_ivar  REAL,
        hei_4471_boxflux  REAL,
        hei_4471_vshift  REAL,
        hei_4471_sigma  REAL,
        hei_4471_cont  REAL,
        hei_4471_cont_ivar  REAL,
        hei_4471_ew  REAL,
        hei_4471_ew_ivar  REAL,
        hei_4471_flux_limit  REAL,
        hei_4471_ew_limit  REAL,
        hei_4471_chi2  REAL,
        hei_4471_npix  INTEGER,
        heii_4686_amp  REAL,
        heii_4686_amp_ivar  REAL,
        heii_4686_flux  REAL,
        heii_4686_flux_ivar  REAL,
        heii_4686_boxflux  REAL,
        heii_4686_vshift  REAL,
        heii_4686_sigma  REAL,
        heii_4686_cont  REAL,
        heii_4686_cont_ivar  REAL,
        heii_4686_ew  REAL,
        heii_4686_ew_ivar  REAL,
        heii_4686_flux_limit  REAL,
        heii_4686_ew_limit  REAL,
        heii_4686_chi2  REAL,
        heii_4686_npix  INTEGER,
        hbeta_amp  REAL,
        hbeta_amp_ivar  REAL,
        hbeta_flux  REAL,
        hbeta_flux_ivar  REAL,
        hbeta_boxflux  REAL,
        hbeta_vshift  REAL,
        hbeta_sigma  REAL,
        hbeta_cont  REAL,
        hbeta_cont_ivar  REAL,
        hbeta_ew  REAL,
        hbeta_ew_ivar  REAL,
        hbeta_flux_limit  REAL,
        hbeta_ew_limit  REAL,
        hbeta_chi2  REAL,
        hbeta_npix  INTEGER,
        oiii_4959_amp  REAL,
        oiii_4959_amp_ivar  REAL,
        oiii_4959_flux  REAL,
        oiii_4959_flux_ivar  REAL,
        oiii_4959_boxflux  REAL,
        oiii_4959_vshift  REAL,
        oiii_4959_sigma  REAL,
        oiii_4959_cont  REAL,
        oiii_4959_cont_ivar  REAL,
        oiii_4959_ew  REAL,
        oiii_4959_ew_ivar  REAL,
        oiii_4959_flux_limit  REAL,
        oiii_4959_ew_limit  REAL,
        oiii_4959_chi2  REAL,
        oiii_4959_npix  INTEGER,
        oiii_5007_amp  REAL,
        oiii_5007_amp_ivar  REAL,
        oiii_5007_flux  REAL,
        oiii_5007_flux_ivar  REAL,
        oiii_5007_boxflux  REAL,
        oiii_5007_vshift  REAL,
        oiii_5007_sigma  REAL,
        oiii_5007_cont  REAL,
        oiii_5007_cont_ivar  REAL,
        oiii_5007_ew  REAL,
        oiii_5007_ew_ivar  REAL,
        oiii_5007_flux_limit  REAL,
        oiii_5007_ew_limit  REAL,
        oiii_5007_chi2  REAL,
        oiii_5007_npix  INTEGER,
        nii_5755_amp  REAL,
        nii_5755_amp_ivar  REAL,
        nii_5755_flux  REAL,
        nii_5755_flux_ivar  REAL,
        nii_5755_boxflux  REAL,
        nii_5755_vshift  REAL,
        nii_5755_sigma  REAL,
        nii_5755_cont  REAL,
        nii_5755_cont_ivar  REAL,
        nii_5755_ew  REAL,
        nii_5755_ew_ivar  REAL,
        nii_5755_flux_limit  REAL,
        nii_5755_ew_limit  REAL,
        nii_5755_chi2  REAL,
        nii_5755_npix  INTEGER,
        hei_5876_amp  REAL,
        hei_5876_amp_ivar  REAL,
        hei_5876_flux  REAL,
        hei_5876_flux_ivar  REAL,
        hei_5876_boxflux  REAL,
        hei_5876_vshift  REAL,
        hei_5876_sigma  REAL,
        hei_5876_cont  REAL,
        hei_5876_cont_ivar  REAL,
        hei_5876_ew  REAL,
        hei_5876_ew_ivar  REAL,
        hei_5876_flux_limit  REAL,
        hei_5876_ew_limit  REAL,
        hei_5876_chi2  REAL,
        hei_5876_npix  INTEGER,
        oi_6300_amp  REAL,
        oi_6300_amp_ivar  REAL,
        oi_6300_flux  REAL,
        oi_6300_flux_ivar  REAL,
        oi_6300_boxflux  REAL,
        oi_6300_vshift  REAL,
        oi_6300_sigma  REAL,
        oi_6300_cont  REAL,
        oi_6300_cont_ivar  REAL,
        oi_6300_ew  REAL,
        oi_6300_ew_ivar  REAL,
        oi_6300_flux_limit  REAL,
        oi_6300_ew_limit  REAL,
        oi_6300_chi2  REAL,
        oi_6300_npix  INTEGER,
        nii_6548_amp  REAL,
        nii_6548_amp_ivar  REAL,
        nii_6548_flux  REAL,
        nii_6548_flux_ivar  REAL,
        nii_6548_boxflux  REAL,
        nii_6548_vshift  REAL,
        nii_6548_sigma  REAL,
        nii_6548_cont  REAL,
        nii_6548_cont_ivar  REAL,
        nii_6548_ew  REAL,
        nii_6548_ew_ivar  REAL,
        nii_6548_flux_limit  REAL,
        nii_6548_ew_limit  REAL,
        nii_6548_chi2  REAL,
        nii_6548_npix  INTEGER,
        halpha_amp  REAL,
        halpha_amp_ivar  REAL,
        halpha_flux  REAL,
        halpha_flux_ivar  REAL,
        halpha_boxflux  REAL,
        halpha_vshift  REAL,
        halpha_sigma  REAL,
        halpha_cont  REAL,
        halpha_cont_ivar  REAL,
        halpha_ew  REAL,
        halpha_ew_ivar  REAL,
        halpha_flux_limit  REAL,
        halpha_ew_limit  REAL,
        halpha_chi2  REAL,
        halpha_npix  INTEGER,
        nii_6584_amp  REAL,
        nii_6584_amp_ivar  REAL,
        nii_6584_flux  REAL,
        nii_6584_flux_ivar  REAL,
        nii_6584_boxflux  REAL,
        nii_6584_vshift  REAL,
        nii_6584_sigma  REAL,
        nii_6584_cont  REAL,
        nii_6584_cont_ivar  REAL,
        nii_6584_ew  REAL,
        nii_6584_ew_ivar  REAL,
        nii_6584_flux_limit  REAL,
        nii_6584_ew_limit  REAL,
        nii_6584_chi2  REAL,
        nii_6584_npix  INTEGER,
        sii_6716_amp  REAL,
        sii_6716_amp_ivar  REAL,
        sii_6716_flux  REAL,
        sii_6716_flux_ivar  REAL,
        sii_6716_boxflux  REAL,
        sii_6716_vshift  REAL,
        sii_6716_sigma  REAL,
        sii_6716_cont  REAL,
        sii_6716_cont_ivar  REAL,
        sii_6716_ew  REAL,
        sii_6716_ew_ivar  REAL,
        sii_6716_flux_limit  REAL,
        sii_6716_ew_limit  REAL,
        sii_6716_chi2  REAL,
        sii_6716_npix  INTEGER,
        sii_6731_amp  REAL,
        sii_6731_amp_ivar  REAL,
        sii_6731_flux  REAL,
        sii_6731_flux_ivar  REAL,
        sii_6731_boxflux  REAL,
        sii_6731_vshift  REAL,
        sii_6731_sigma  REAL,
        sii_6731_cont  REAL,
        sii_6731_cont_ivar  REAL,
        sii_6731_ew  REAL,
        sii_6731_ew_ivar  REAL,
        sii_6731_flux_limit  REAL,
        sii_6731_ew_limit  REAL,
        sii_6731_chi2  REAL,
        sii_6731_npix  INTEGER,
        oii_7320_amp  REAL,
        oii_7320_amp_ivar  REAL,
        oii_7320_flux  REAL,
        oii_7320_flux_ivar  REAL,
        oii_7320_boxflux  REAL,
        oii_7320_vshift  REAL,
        oii_7320_sigma  REAL,
        oii_7320_cont  REAL,
        oii_7320_cont_ivar  REAL,
        oii_7320_ew  REAL,
        oii_7320_ew_ivar  REAL,
        oii_7320_flux_limit  REAL,
        oii_7320_ew_limit  REAL,
        oii_7320_chi2  REAL,
        oii_7320_npix  INTEGER,
        oii_7330_amp  REAL,
        oii_7330_amp_ivar  REAL,
        oii_7330_flux  REAL,
        oii_7330_flux_ivar  REAL,
        oii_7330_boxflux  REAL,
        oii_7330_vshift  REAL,
        oii_7330_sigma  REAL,
        oii_7330_cont  REAL,
        oii_7330_cont_ivar  REAL,
        oii_7330_ew  REAL,
        oii_7330_ew_ivar  REAL,
        oii_7330_flux_limit  REAL,
        oii_7330_ew_limit  REAL,
        oii_7330_chi2  REAL,
        oii_7330_npix  INTEGER,
        siii_9069_amp  REAL,
        siii_9069_amp_ivar  REAL,
        siii_9069_flux  REAL,
        siii_9069_flux_ivar  REAL,
        siii_9069_boxflux  REAL,
        siii_9069_vshift  REAL,
        siii_9069_sigma  REAL,
        siii_9069_cont  REAL,
        siii_9069_cont_ivar  REAL,
        siii_9069_ew  REAL,
        siii_9069_ew_ivar  REAL,
        siii_9069_flux_limit  REAL,
        siii_9069_ew_limit  REAL,
        siii_9069_chi2  REAL,
        siii_9069_npix  INTEGER,
        siii_9532_amp  REAL,
        siii_9532_amp_ivar  REAL,
        siii_9532_flux  REAL,
        siii_9532_flux_ivar  REAL,
        siii_9532_boxflux  REAL,
        siii_9532_vshift  REAL,
        siii_9532_sigma  REAL,
        siii_9532_cont  REAL,
        siii_9532_cont_ivar  REAL,
        siii_9532_ew  REAL,
        siii_9532_ew_ivar  REAL,
        siii_9532_flux_limit  REAL,
        siii_9532_ew_limit  REAL,
        siii_9532_chi2  REAL,
        siii_9532_npix  INTEGER,
        continuum_coeff_0  DOUBLE PRECISION,
        continuum_coeff_1  DOUBLE PRECISION,
        continuum_coeff_2  DOUBLE PRECISION,
        continuum_coeff_3  DOUBLE PRECISION,
        continuum_coeff_4  DOUBLE PRECISION,
        continuum_coeff_5  DOUBLE PRECISION,
        continuum_coeff_6  DOUBLE PRECISION,
        continuum_coeff_7  DOUBLE PRECISION,
        continuum_coeff_8  DOUBLE PRECISION,
        continuum_coeff_9  DOUBLE PRECISION,
        target  TEXT,
        survey  TEXT,
        PRIMARY KEY (targetid, survey, target)
        );
    """]


    @staticmethod
    def fill_table():   
        surveys=["sv3","main"]
        targets = ["bright","dark"]
        for survey in surveys:
            for target in targets:
                print(survey, target)
                filename_spec = f"/global/cfs/cdirs/desi/spectro/fastspecfit/everest/catalogs/fastspec-everest-{survey}-{target}.fits"

                try:               
                    dat = Table.read(filename_spec, format='fits',hdu=1)
                except:
                    print(f"{filename} not found")
                    continue
                
                for icoeff in range(0,10):
                    dat[f'CONTINUUM_COEFF_{icoeff}']= dat['CONTINUUM_COEFF'][0:len(dat),icoeff]
                dat.remove_column('CONTINUUM_COEFF')
                    
                
                dat.convert_bytestring_to_unicode()
                df = dat.to_pandas()
                df['target']=numpy.full(df.shape[0],target)
                df['survey']=numpy.full(df.shape[0],survey)
                df.columns= df.columns.str.lower()
                try:
                    df.to_sql('fastspecfit_spec',engine,index=False,if_exists='append',schema='everest',chunksize=10000)
                except:
                    dtypesToSchema(df.dtypes)
                    sys.exit()

class fastspecfit_phot(tablebaseclass):
    tablename = "everest.fastspecfit_phot"
    createstatements = [ """
        CREATE TABLE IF NOT EXISTS everest.fastspecfit_phot ( 
        targetid  BIGINT,
        continuum_chi2  REAL,
        continuum_age  REAL,
        continuum_av  REAL,
        continuum_av_ivar  REAL,
        dn4000_model  REAL,
        flux_synth_model_g  REAL,
        flux_synth_model_r  REAL,
        flux_synth_model_z  REAL,
        flux_synth_model_w1  REAL,
        flux_synth_model_w2  REAL,
        kcorr_u  REAL,
        absmag_u  REAL,
        absmag_ivar_u  REAL,
        kcorr_b  REAL,
        absmag_b  REAL,
        absmag_ivar_b  REAL,
        kcorr_v  REAL,
        absmag_v  REAL,
        absmag_ivar_v  REAL,
        kcorr_sdss_u  REAL,
        absmag_sdss_u  REAL,
        absmag_ivar_sdss_u  REAL,
        kcorr_sdss_g  REAL,
        absmag_sdss_g  REAL,
        absmag_ivar_sdss_g  REAL,
        kcorr_sdss_r  REAL,
        absmag_sdss_r  REAL,
        absmag_ivar_sdss_r  REAL,
        kcorr_sdss_i  REAL,
        absmag_sdss_i  REAL,
        absmag_ivar_sdss_i  REAL,
        kcorr_sdss_z  REAL,
        absmag_sdss_z  REAL,
        absmag_ivar_sdss_z  REAL,
        kcorr_w1  REAL,
        absmag_w1  REAL,
        absmag_ivar_w1  REAL,
        continuum_coeff_0  DOUBLE PRECISION,
        continuum_coeff_1  DOUBLE PRECISION,
        continuum_coeff_2  DOUBLE PRECISION,
        continuum_coeff_3  DOUBLE PRECISION,
        continuum_coeff_4  DOUBLE PRECISION,
        continuum_coeff_5  DOUBLE PRECISION,
        continuum_coeff_6  DOUBLE PRECISION,
        continuum_coeff_7  DOUBLE PRECISION,
        continuum_coeff_8  DOUBLE PRECISION,
        continuum_coeff_9  DOUBLE PRECISION,
        target  TEXT,
        survey  TEXT,
        PRIMARY KEY (targetid, survey, target)
        );
    """]


    @staticmethod
    def fill_table():   
        surveys=["sv3","main"]
        targets = ["bright","dark"]
        for survey in surveys:
            for target in targets:
                print(survey, target)
                filename_photo = f"/global/cfs/cdirs/desi/spectro/fastspecfit/everest/catalogs/fastphot-everest-{survey}-{target}.fits"
                try:
                    dat = Table.read(filename_photo, format='fits',hdu=1)
                except:
                    print(f"{filename} not found")
                    continue
                for icoeff in range(0,10):
                    dat[f'CONTINUUM_COEFF_{icoeff}']= dat['CONTINUUM_COEFF'][0:len(dat),icoeff]
                dat.remove_column('CONTINUUM_COEFF')
                
                dat.convert_bytestring_to_unicode()
                df = dat.to_pandas()
                df['target']=numpy.full(df.shape[0],target)
                df['survey']=numpy.full(df.shape[0],survey)
                df.columns= df.columns.str.lower()
                try:
                    df.to_sql('fastspecfit_phot',engine,index=False,if_exists='append',schema='everest',chunksize=10000)
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
        everest_daily_fibermap.fill_table()

    elif args.action == "load-mosthosts":
        if not args.really_do:
            print( "Running this will delete the current mosthosts.mosthosts\n"
                   "table and make a new one from the specified file.\n"
                   "If this is in fact what you want to do, rerun with --really-do" )
        else:
            mosthosts.create_table( overwrite=True )
            mosthosts.fill_table( args.mosthosts_file )
    elif args.action == "load-fastspecfit":
        if not args.really_do:
            print( "Running this will delete the current fastspecphoto_data\n"
                   "table and make a new one from the specified file.\n"
                   "If this is in fact what you want to do, rerun with --really-do" )
        else:
            fastspecfit_meta.create_table(overwrite=True )
            fastspecfit_meta.fill_table( )
#             fastspecfit_spec.create_table(overwrite=True )
#             fastspecfit_spec.fill_table( )
#             fastspecfit_phot.create_table(overwrite=True )
#             fastspecfit_phot.fill_table( )
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
    
