#!/usr/bin/env python

import argparse
from timedomain.filters import *
from timedomain.iterators import *
from timedomain.sp_utils import *
from timedomain.fs_utils import *
import timedomain.config as config
import sys
from astropy.table import Table
from desiutil.log import get_logger, DEBUG
log = get_logger(DEBUG)


def main(args):
    """ Main entry point of the app """
    print("Start ", args)
    logic = getattr(sys.modules[__name__], args.logic)
    iterator_= getattr(sys.modules[__name__], args.iterator)

    prunelogic = args.logic[0:args.logic.find("Logic")]

    ### Get the tile and array from the arguments

    if args.obsdates_tilenumbers!=None:
        obsdates_tilenumbers_str = args.obsdates_tilenumbers
        obsdates_tilenumbers = np.chararray((len(obsdates_tilenumbers_str),2),itemsize=10,unicode=True)
        for i in range(len(obsdates_tilenumbers_str)):
            obsdates_tilenumbers[i,:]=obsdates_tilenumbers_str[i].split('|')
        log.info(obsdates_tilenumbers)
    date = []
    tile = []
    for obsdate,tile_number in obsdates_tilenumbers:
        if obsdate >= config.mindate:
            date.append(obsdate)
            tile.append(tile_number)

#     iterator = TileDate_SpectraPairs_Iterator(tile, date, subdir=args.subdir,trunk=args.trunk, verbose=True)
#     iterator = TileDate_TargetPairs_Iterator(tile, date, subdir=args.subdir,trunk=args.trunk, verbose=True)
    iterator = iterator_(tile, date, subdir=args.subdir,trunk=args.trunk, verbose=True)

    # make this none for results to appear in the notebook

    for (pspectra, meta) in iterator:
        pspectra0,pspectra1 = pspectra
        spdf = ["diff",logic.__name__,args.subdir,args.trunk,meta[0]['date']]

        # This should work in the far future
#         zf = os.path.join(fs_utils.redux,args.subdir,'tiles/cumulative',meta[0]['tile'])
#         zdirs=glob(zf+'/*/')
#         zdirs.sort()
#         zf = glob(zdirs[-1]+'/zbest-{}-{}-*.fits'.format(meta[0]['panel'],meta[0]['tile']))[0]

        zf = fitsfile(meta[0]['tile'],meta[0]['date'],meta[0]['panel'],subdir=args.subdir,trunk='zbest')
#         zf1 = fitsfile(meta[1]['tile'],meta[1]['date'],meta[1]['panel'],subdir='daily',trunk='zbest')

#         # For now consider missing zbest files as catastrophic errors
        if (zf is None):
            log.warning("Missing zbest")
            sys.exit(1)
            
        zbest = Table.read(zf, 'ZBEST')
        
        # which of these are real targets
        triggered, diff = logic.filter(pspectra0,pspectra1, zbest, norm=True,ston_cut=5)

#         # plot triggered objects
        if triggered.sum() > 0:

            wheretriggered = np.where(triggered)[0]

            for sig in wheretriggered.flat:

                targetid = diff.fibermap['TARGETID'].data[sig].astype('str')
                log.info("TARGETID {}".format(targetid))
#                 SkyPortal.nukeCandidate(targetid, "DESIDIFF {}".format(prunelogic))
                
                # save the DESI-inferred redshift
#                 zf = fitsfile(meta[0]['tile'],meta[0]['date'],meta[0]['panel'],subdir='daily',trunk='zbest')
#                 zbest = Table.read(zf, 'ZBEST')
                zin = numpy.where(zbest['TARGETID'] == diff.fibermap['TARGETID'].data[sig])[0]
                zin=zin[0]
 
                # save the candidate
                data_override = {
                  "origin": "DESIDIFF {}-{}".format(args.iterator, args.logic),
                  "altdata": {'logic':args.logic, 'iterator':args.iterator, 'subdir':args.subdir, 'trunk':args.trunk,
                              'tile':meta[0]['tile'], 'date0':meta[0]['date'],'date1':meta[1]['date'], 'panel':meta[0]['panel'],
                              'fibermap': {'Z':str(zbest['Z'][zin]), 'ZERR':str(zbest['ZERR'][zin]), 'SPECTYPE':zbest['SPECTYPE'][zin]}}
                }

                SkyPortal.postCandidate(sig, diff.fibermap, "DESIDIFF {}".format(prunelogic),data_override=data_override)
   
                # Annotate
                data_override = {
                    "group_ids": [SkyPortal.group_id('DESI'), SkyPortal.group_id('Wayne State')],
                    "data": {'Z':str(zbest['Z'][zin]), 'ZERR':str(zbest['ZERR'][zin]), 'SPECTYPE':zbest['SPECTYPE'][zin]}
                }

                SkyPortal.postAnnotation(sig, diff.fibermap,data_override=data_override)

                data = {
                    "redshift": str(zbest['Z'][zin])
                }
                SkyPortal.api('PATCH', '{}/api/sources/DESI{}'.format(SkyPortal.url,targetid),data=data)
            
                # save Spectra
                targetid = diff.fibermap['TARGETID'].data[sig].astype('str')
                data_override = {
                    "origin": "{}-{}".format(meta[0]['date'],meta[1]['date']),
                    "group_ids": [SkyPortal.group_id('DESI'), SkyPortal.group_id('Wayne State')],
                    "altdata": {'logic':args.logic, 'iterator':args.iterator, 'subdir':args.subdir, 'trunk':args.trunk,
                    'tile':meta[0]['tile'], 'date0':meta[0]['date'],'date1':meta[1]['date'], 'panel':meta[0]['panel']}
                }              
                SkyPortal.postSpectra(targetid, diff, data_override=data_override,coadd_camera=True)
                    
                data_override = {
                    "group_ids": [SkyPortal.group_id('DESI'), SkyPortal.group_id('Wayne State')],
                    "altdata": {'subdir':args.subdir, 'trunk':args.trunk,
                              'tile':meta[0]['tile'], 'date':meta[0]['date'], 'panel':meta[0]['panel']}
                }                   
                SkyPortal.postSpectra(targetid, pspectra0, data_override=data_override,coadd_camera=True)

                data_override = {
                    "group_ids": [SkyPortal.group_id('DESI'), SkyPortal.group_id('Wayne State')],
                    "altdata": {'subdir':args.subdir, 'trunk':args.trunk,
                              'tile':meta[1]['tile'], 'date':meta[1]['date'], 'panel':meta[1]['panel']}
                } 
                SkyPortal.postSpectra(targetid, pspectra1, data_override=data_override,coadd_camera=True)
                                    
                # get all photometry associated with this guy.  Note that the altdata is for the 3 filters.
                # SkyPortal will be fixed so that a scalar will be broadcast to arrays
                # this can be made more efficient.
                alldates = tileToDates(meta[0]['tile'], subdir=args.subdir)
                for d in alldates:
                    ffile = fs_utils.fitsfile(meta[0]['tile'], d, meta[0]['panel'], subdir=args.subdir,trunk=args.trunk)
                    if ffile is not None:
                        sp = read_spectra(ffile)[sig]
                        if sp.fibermap['FIBERSTATUS'][0] ==0:
                            data_override = {
                                "group_ids": [SkyPortal.group_id('DESI'), SkyPortal.group_id('Wayne State')],
                                "altdata": {'subdir':args.subdir, 'trunk':args.trunk,
                                    'tile':meta[0]['tile'], 'date':d,'panel': meta[0]['panel']}
                            }

                            # Modify this so that it gets all phases, not just the pair
                            SkyPortal.postPhotometry(targetid, sp,coadd_camera=True, data_override=data_override)
                logic.plotter(sig,pspectra0, pspectra1, diff, savepdf=spdf)
                
                wefw

    print("End")
                
if __name__ == "__main__":
    
    # ./diff.py 20201223 CVLogic Date_SpectraPairs_Iterator daily coadd
    # ./diff.py 20201223 CVLogic Date_TargetPairs_Iterator daily spectra
    
# date = "20201223"
# subdir = 'daily'
# trunk='coadd'
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("iterator", help="Required positional argument")
    parser.add_argument("logic", help="Required positional argument")    
    parser.add_argument("subdir", help="Required positional argument")
    parser.add_argument("trunk", help="Required positional argument")

    
    #If there are more than 1 obsdate, provide a 2D array
    parser.add_argument('-o', '--obsdates_tilenumbers', nargs='+', type=str,default=None,
                        help='str array with columns obsdate, tilenumber, separated by |')
    

    args = parser.parse_args()
    main(args)
