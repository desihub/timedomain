import numpy as np
import matplotlib.path

def insideFootprint (ra, dec) :
    ix = ra > 180
    ra[ix] = ra[ix]-360.
    footprint = getFootprint()
    ix = footprint.contains_points( zip(ra,dec) )
    return ix

def getFootprint(blob=1) :
    ra, dec = getFootprintRaDec(blob)
    footprint = desPath(ra,dec)
    return footprint

def getFootprintRaDec(blob=1) :
    import os
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    footFile = gw_data_dir + "desi_footprint_{}.txt".format(blob)
    ra,dec = np.genfromtxt(footFile,unpack=True,comments="#")
    #ix, = np.where(ra> 180)
    #ra[ix] = ra[ix] - 360
    return ra,dec

def desPath(raDes, decDes) :
    footprint = matplotlib.path.Path(zip(raDes, decDes))
    return footprint

