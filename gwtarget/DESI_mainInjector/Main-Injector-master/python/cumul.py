# countour_sum:
#   find the number of pixels that contain some total
#   sum of vals (if val is normed to 1,
#   this is a percentage) 
# area
#   one can convert this to area
#
import healpy as hp
import numpy as np
def area(ra,dec,vals, threshold, nside, max_area) :
    npix, cumulative  = contour_sum(ra,dec,vals,threshold)
    if npix == ra.size :
        area = max_area
    else :
        area_per_pixel = hp.nside2pixarea(nside)*((360./2/np.pi)**2)
        area = npix*area_per_pixel
    return area, cumulative

def probability_covered(ra,dec,vals, area, nside, max_area) :
    area_per_pixel = hp.nside2pixarea(nside)*((360./2/np.pi)**2)
    npix= area/area_per_pixel
    ix = np.argsort(vals)[::-1] ;# descending order
    probability = (vals[ix][0:npix]).sum()
    return probability

def contour_sum(ra, dec, vals, threshold) :
    ix = np.argsort(vals)[::-1] ;# descending order
    ra_sort = ra[ix]
    dec_sort = dec[ix]
    vals_sort = vals[ix]

    cumulative = 0
    i = 0
    while cumulative < threshold and i < ra.size :
        cumulative +=  vals_sort[i]
        i += 1
    if i == ra.size :
        i -= 1
        print "\t contour_sim: ",
        print " npix set to max as sum < threshold, ",
        print "{:.2f}<{:.2f} ".format(vals.sum(),threshold)
    npixels =ra_sort[:i].size 
    return npixels,cumulative
