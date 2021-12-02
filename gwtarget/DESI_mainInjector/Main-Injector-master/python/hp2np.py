import healpy as hp
import numpy as np

license="""
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

#
# from a healpix map, 
# return  ra, dec, val
#   fluxConservation = False => averaging maps when changing resolution
#   fluxConservation = True => sums maps when changing resolution
#
def hp2np (hp_map_file, nan=True, degrade=False, fluxConservation=True, field=0, verbose=False) :
    hm = hp.read_map(hp_map_file, field=field, verbose=verbose)
    ra,dec,vals = map2np(hm, resolution=degrade, fluxConservation=fluxConservation, verbose=verbose)
    return ra,dec,vals

#
# healpix maps have implicit ra,dec-
# one is supposed to know from position along
# vector what the ra,dec is.
# Often this is inconvenient. 
#
# Give me, ra, dec, val.
def map2np (hp_map, resolution=False, fluxConservation=True, verbose=False) :
    nside = hp.npix2nside(len(hp_map))
    if verbose: print "\t map2np: \t res= ", nside
    if resolution :
        if fluxConservation :
            hp_map = hp.ud_grade(hp_map, resolution, power=-2)
        else :
            hp_map = hp.ud_grade(hp_map, resolution)
        nside = hp.npix2nside(len(hp_map))
        if verbose: print "map2np \t changed resolution to ", nside
    ix = range(0,hp_map.size)
    # pix2and wants the indicies numbers of the pixel to get the coord of
    theta,phi = hp.pix2ang(nside,ix)
    theta = theta*360./2./np.pi; 
    phi = phi*360./2./np.pi
    ra = phi
    dec = 90-theta
    ix=np.nonzero(ra > 180); ra[ix]=ra[ix]-360.
    return ra,dec,hp_map

#
# convert a non-hp map into a hp map
#
def convert(template_ra, template_dec, ra, dec, vals) :
    import scipy.spatial
    tree = scipy.spatial.KDTree( zip(ra,dec) )
    newmap = []
    for i in range(0,template_ra.size) :
        data = tree.query(np.array([template_ra[i], template_dec[i]]))
        index = data[1]
        newmap.append(vals[index])
    newmap = np.array(newmap)
    return newmap
        
def radec2hp(ra, dec, nsides=1024) :
# coordinate conversion
    phi = ra*2*np.pi/360.;
    theta = (90-dec)*2*np.pi/360.
    pix = hp.ang2pix(nsides,theta,phi)
# build an empty healpix map
    npix = hp.nside2npix(nsides)
    vector = np.empty(npix)
    vector[:] = np.NAN
# fill it
    unique_pixels  = np.unique(pix)
    for up in unique_pixels:
        ix = np.nonzero(pix==up)
        count = ix[0].size
        vector[up] = count
    return vector
