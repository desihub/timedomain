license="""
   Copyright (C) 2015 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style. (Note the rate at which people who stand on the
   summit of K2 to successfully make it down.)

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

import numpy as np
import decam2hp

#
# Example usage
#
#    obs.resetTime(mjd+time)
#    obs.limitMag("i", exposure=90)
#
#    sm=sourceProb.map(obs);
#    models_at_t = modelsAtTimeT (models, time)
#    abs_mag = models_at_t[0]
#    sm.modelAbsoluteMagnitude = abs_mag
#    sm.searchDistance = np.array([distance,])
#    sm.calculateProb()
#
#    raHex, decHex, idHex, hexVals, rank = cutAndHexalate (obs, sm, raHexen, decHexen)

#
# Given an obs object and a sm object, along with a map of hex centers,
# this routine returns the sum of the probability inside of DECam camera outlines
# for all hexes in the map, along with a rank of the hexes sorted from large to small.
#
# domain knowlege- the blanco cuts keep everything inside an hour angle range
def hexalateNHexes (obs, sm, nHexes, allskyDesHexes) :
    raHexen, decHexen, idHexen, hexVals, rank = cutAndHexalate(obs, sm, allskyDesHexes)
    raHexen = raHexen[0:nHexes]
    decHexen = decHexen[0:nHexes]
    idHexen = idHexen[0:nHexes]
    hexVals = hexVals[0:nHexes]
    rank = rank[0:nHexes]
    return raHexen, decHexen, idHexen, hexVals, rank

# be aware that while most of my ra,dec are in degrees,
# those in obs and sm are in radians
def cutAndHexalate (obs, sm, camera, hexFile) :
    #allskyDesHexes="../data/all-sky-hexCenters-"+camera+".txt"
    raHexen, decHexen, idHexen = getHexCenters(hexFile)
    raHexen, decHexen, idHexen, hexVals, rank = cutAndHexalateOnRaDec(
        obs, sm, raHexen, decHexen, idHexen, camera)
    return raHexen, decHexen, idHexen, hexVals, rank


def cutAndHexalateOnRaDec (obs, sm, raHexen, decHexen, idHexen, tree, camera, cutProbs=False) :
    #cutProbs calls without Overlap
    verbose = False
    obsHourAngle = obs.ha*360./(2*np.pi)
    obsRa        = obs.ra*360./(2*np.pi)
    obsDec       = obs.dec*360./(2*np.pi)
    # based on blanco horizen limits  (may not be needed with tree)                                                                                                                          
    ix = (abs(obsHourAngle) <= 83. ) & (obsDec < 43.)
    ix2 = decHexen < 43.

    probabilities = obs.map*sm.probMap
    if verbose  :
        print "\t cutAndHexalate probabilities sum",probabilities.sum()

    hexVals = np.zeros(raHexen.size)
    if not cutProbs:
        hexVals[ix2] = decam2hp.hexalateMap(obsRa,obsDec, probabilities, tree,
                                            raHexen[ix2], decHexen[ix2], camera, verbose=False)
    else:
        hexVals[ix2] = decam2hp.hexalateMapWithoutOverlap(obsRa,obsDec, probabilities, tree,
                                            raHexen[ix2], decHexen[ix2], camera, verbose=False)
    if verbose  :
        print "hexVals max", hexVals.max()
    rank=np.argsort(hexVals);
    rank = rank[::-1];# sort from large to small by flipping natural argsort order                                                                                                   
 
    return raHexen, decHexen, idHexen, hexVals, rank


def getHexCenters (hexFile) :
    #allskyDesHexes="../data/all-sky-hexCenters-"+camera+".txt"
    ra,dec = np.genfromtxt(hexFile, unpack=True, \
        usecols=(0,1),comments="#")
    hex_id = getHexId(ra,dec)
    return ra,dec,hex_id

#
# a modified DES convention: just rounded ra,dec with + or - sign
# of type "24-31", "179+10"
def getHexId(ra, dec) :
    intra = np.round(ra)
    intdec = np.abs(np.round(dec));
    sign=np.full(ra.size, "+",dtype="str");
    ix = dec < 0; sign[ix]="-";
    id = np.array([],dtype="str")
    for i in range(0,ra.size) :
        name = str(np.int(intra[i])) + sign[i] + str(np.int(intdec[i]))
        id = np.append(id, name)
    return id


