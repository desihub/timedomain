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
    iterator = TileDate_TargetPairs_Iterator(tile, date)

    # make this none for results to appear in the notebook
#     spdf = ["diff",logic.__name__,args.subdir,args.trunk,args.date]
    for pspectra0,pspectra1 in iterator(args.date,subdir=args.subdir,trunk=args.trunk, verbose=True):
        print(pspectra0,pspectra1)

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