#!/usr/bin/env python

import argparse
from timedomain.filters import *
from timedomain.iterators import *
import sys


__version__=0.1

def main(args):
    """ Main entry point of the app """
    print("Start ", args)
    logic = getattr(sys.modules[__name__], args.logic)
    iterator = getattr(sys.modules[__name__], args.iterator)

    # make this none for results to appear in the notebook
    spdf = ["diff",logic.__name__,args.subdir,args.trunk,args.date]
    for pspectra0,pspectra1 in iterator(args.date,subdir=args.subdir,trunk=args.trunk, verbose=True):
        # which of these are real targets
        triggered, diff = logic.filter(pspectra0,pspectra1, norm=True,ston_cut=5)

        # plot triggered objects
        if triggered.sum() > 0:

            wheretriggered = np.where(triggered)[0]

            for sig in wheretriggered.flat:
#                 SkyPortal.postCandidate(sig, diff.fibermap)
#                 targetid = diff.fibermap['TARGETID'].data[sig].astype('str')
#                 SkyPortal.postSpectra(targetid, diff)
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
    parser.add_argument("date", help="Required positional argument")
    parser.add_argument("logic", help="Required positional argument")    
    parser.add_argument("iterator", help="Required positional argument")
    parser.add_argument("subdir", help="Required positional argument")
    parser.add_argument("trunk", help="Required positional argument")
    
    # Optional argument flag which defaults to False
#     parser.add_argument('-f', '--flag', action="store_true", default=False)

    # Optional argument which requires a parameter (eg. -d test)
#     parser.add_argument("-n", "--name", action="store", dest="name")
#     parser.add_argument('-i','--iargs', nargs='+', action="store", dest="iargs")


    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
#     parser.add_argument(
#         '-v',
#         '--verbose',
#         action='count',
#         default=0,
#         help="Verbosity (-v, -vv, etc)")

#     # Specify output of '--version'
#     parser.add_argument(
#         '--version',
#         action='version',
#         version='%(prog)s (version {version})'.format(version=__version__))

    args = parser.parse_args()
    main(args)