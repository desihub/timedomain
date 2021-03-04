#!/usr/bin/env python

import argparse
from timedomain.filters import *
from timedomain.iterators import *
from timedomain.sp_utils import *
import sys


__version__=0.1

def main(args):
    """ Main entry point of the app """
    print("Start ", args)
    logic = getattr(sys.modules[__name__], args.logic)
    
    ### Get the tile and array from the arguments

    if args.obsdates_tilenumbers!=None:
        obsdates_tilenumbers_str = args.obsdates_tilenumbers
        obsdates_tilenumbers = np.chararray((len(obsdates_tilenumbers_str),2),itemsize=10,unicode=True)
        for i in range(len(obsdates_tilenumbers_str)):
            obsdates_tilenumbers[i,:]=obsdates_tilenumbers_str[i].split('|')
        print(obsdates_tilenumbers)
    date = []
    tile = []
    for obsdate,tile_number in obsdates_tilenumbers:
        date.append(obsdate)
        tile.append(tile_number)  

    iterator = TileDate_TargetPairs_Iterator(tile, date)

    # make this none for results to appear in the notebook
#     spdf = ["diff",logic.__name__,args.subdir,args.trunk,args.date]

    for pspectra0,pspectra1 in iterator:

        # which of these are real targets
        triggered, diff = logic.filter(pspectra0,pspectra1, norm=True,ston_cut=5)

        # plot triggered objects
        if triggered.sum() > 0:

            wheretriggered = np.where(triggered)[0]

            for sig in wheretriggered.flat:
                SkyPortal.postCandidate(sig, diff.fibermap, "DESIDIFF_CV_daily")
                targetid = diff.fibermap['TARGETID'].data[sig].astype('str')
                SkyPortal.postSpectra(targetid, diff)
                SkyPortal.postSpectra(targetid, pspectra0)
                SkyPortal.postSpectra(targetid, pspectra1)
                logic.plotter(sig,pspectra0, pspectra1, diff, savepdf=spdf)

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
    parser.add_argument("logic", help="Required positional argument")    
    parser.add_argument("subdir", help="Required positional argument")
    parser.add_argument("trunk", help="Required positional argument")
    
    #If there are more than 1 obsdate, provide a 2D array
    parser.add_argument('-o', '--obsdates_tilenumbers', nargs='+', type=str,default=None,
                        help='str array with columns obsdate, tilenumber, separated by |')
    

    args = parser.parse_args()
    main(args)