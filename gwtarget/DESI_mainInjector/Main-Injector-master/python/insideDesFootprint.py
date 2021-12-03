import numpy as np
import matplotlib.path

def insideFootprint (ra, dec) :
    ix = ra > 180
    ra[ix] = ra[ix]-360.
    footprint = getFootprint()
    ix = footprint.contains_points( zip(ra,dec) )
    return ix

def getFootprint() :
    ra, dec = getFootprintRaDec()
    footprint = desPath(ra,dec)
    return footprint

def getFootprintRaDec() :
    import os
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    footFile = gw_data_dir + "round19_v0.txt"
    #ra,dec = np.genfromtxt(footFile,unpack=True,skiprows=30)
    ra,dec = np.genfromtxt(footFile,unpack=True,comments="#")
    return ra,dec

def desPath(raDes, decDes) :
    footprint = matplotlib.path.Path(zip(raDes, decDes))
    return footprint

