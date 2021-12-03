import numpy as np
import os
import rotate
from scipy.interpolate import interp1d

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

# McBryde Thomas V, the flat polar quartic
#  http://books.google.com/books?id=8LZeu8RxOIsC&pg=PA151&lpg=PA151&dq=McBryde-Thomas+Flat-Polar+Quartic+equations&source=bl&ots=Amq75MbxBD&sig=5YJQoCpr3FDXxz2k_W58T0B94I8&hl=en&sa=X&ei=uakSUtzUEqH94APUy4Aw&ved=0CDsQ6AEwAg#v=onepage&q=McBryde-Thomas%20Flat-Polar%20Quartic%20equations&f=false
# which is 153, Cartographic Science: A Compendium of Map Projections, with Derivations By Donald Fenna

#
# alpha is the degrees to rotate ra
# beta is the degrees to rotate dec
#
def mcbryde (ra,dec, southUp=0, test=0, alpha=0, beta=0, isLine=False) :
    global psi_spline
    if test: rai = ra; deci=dec
    pi = np.pi
    R = 100.

    #alpha = -39
    #alpha = -30 ;# for a projection on a nice part of the des
    #alpha = 0
    #alpha = -80
    ra, dec = mtCoord (ra,dec, alpha, beta) 
    # the alpha, beta roations play havoc on lines but not maps
    if isLine :
        ix = np.argsort(ra)
        ra = ra[ix]; dec=dec[ix]

    try :
        psi = psi_spline(dec)
    except :
        readMT()
        psi = psi_spline(np.array(dec,dtype='float'))
    ra  =  ra*2*pi/360.
    dec = dec*2*pi/360.
    psi = psi*2*np.pi/360.

    psi = psi/2.
    #x = ra * np.sqrt(6./7.) * ( 2*np.cos(2*psi) -1)
    #y = (9./np.sqrt(7.))*np.sin(psi)

    a = np.sqrt(2 + np.sqrt(2))
    x = 1 + 2*np.cos(2*psi)/np.cos(psi)
    x = ra*x/(np.sqrt(3.)*a)
    y = 2*np.sqrt(3)*np.sin(psi)/a

    x = x*R
    y = y*R


    if southUp : 
        y = -y
    if type(x) is np.array   and len(x) == 1 : x = x[0]; y = y[0]
    if type(x) is np.ndarray and len(x) == 1 : x = x[0]; y = y[0]
    return x,y

    if test :
        print "psi= ", psi*2*360/(2*np.pi)
        print "ra,dec= ", rai, deci
        print "x,y= ", round(x,3), round(y,3)
        print "-180,  0 and 180,  0"
        alpha = 0; ra = -180; psi=alpha/2.; psi = psi*2*np.pi/360.; ra = ra*2*np.pi/360.
        #print np.cos(psi),  np.cos(2*psi), np.sin(psi)
        #       1.0          1.0                 0
        #q1 = 0.31246*R*ra*(1 + 2* np.cos(psi*2)/np.cos(psi))
        #r1 = 1.87476*R*np.sin(psi)
        q1 = 0.31246*R*ra*3
        r1 = 0.
        alpha = 0; ra = 180; psi=alpha/2.; psi = psi*2*np.pi/360.; ra = ra*2*np.pi/360.
        #print np.cos(psi),  np.cos(2*psi), np.sin(psi)
        #q1 = 0.31246*R*ra*(1 + 2* np.cos(psi*2)/np.cos(psi))
        #r1 = 1.87476*R*np.sin(psi)
        q2 = 31.246*ra*3
        r2 = 0
        print round(q1,3), round(r1,3), round(q2,3), round(r2,3)
    
        print "-180,-90 and 180,-90"
        alpha = -90; ra = -180; psi=alpha/2.; psi = psi*2*np.pi/360.; ra = ra*2*np.pi/360.
        #print np.cos(psi),  np.cos(2*psi), np.sin(psi)
        #       0.707         0.0                 -0.707
        #q1 = 0.31246*R*ra*(1 + 2* np.cos(psi*2)/np.cos(psi))
        #r1 = 1.87476*R*np.sin(psi)
        q1 = 0.31246*R*ra*(1 + 2* 0)
        r1 = 1.87476*R*-0.707
        alpha = -90; ra = 180; psi=alpha/2.; psi = psi*2*np.pi/360.; ra = ra*2*np.pi/360.
        #print np.cos(psi),  np.cos(2*psi), np.sin(psi)
        #q2 = 0.31246*R*ra*(1 + 2* np.cos(psi*2)/np.cos(psi))
        #r2 = 1.87476*R*np.sin(psi)
        q2 = 31.246*ra
        r2 =-132.545
        print round(q1,3), round(r1,3), round(q2,3), round(r2,3)


        alpha = -75; rao = 120; 
        psi=alpha; 
        psi = psi_spline(alpha)
        psi=psi/2.; 
        psi = psi*2*np.pi/360.; ra = rao*2*np.pi/360.
        q1 = 0.31246*R*ra*(1 + 2* np.cos(psi*2)/np.cos(psi))
        r1 = 1.87476*R*np.sin(psi)
        print("alt x,y for {:d} {:d}= {:8.3f} {:8.3f}".format(rao,alpha, round(q1,3), round(r1,3)))

        alpha = -85; rao = 120; 
        psi=alpha; 
        psi = psi_spline(alpha)
        psi=psi/2.; 
        psi = psi*2*np.pi/360.; ra = rao*2*np.pi/360.
        q1 = 0.31246*R*ra*(1 + 2* np.cos(psi*2)/np.cos(psi))
        r1 = 1.87476*R*np.sin(psi)
        print("alt x,y for {:d} {:d}= {:8.3f} {:8.3f}".format(rao,alpha, round(q1,3), round(r1,3)))
        print("cal x,y for {:d} {:d}= {:8.3f} {:8.3f}".format(rai,deci, round(x,3), round(y,3)))


def mtCoord (ra,dec, alpha, beta) :
    lon = ra; lat = dec
    x,y,z = rotate.sphericalToCartesian(lon,lat)
    x,y,z = rotate.rotateAboutZaxis(x,y,z, alpha)
    x,y,z = rotate.rotateAboutYaxis(x,y,z, beta)
    lon,lat,r = rotate.cartesianToSpherical(x,y,z)
    try :
        index = np.nonzero(lon > 180)
        lon[index] = lon[index]-360.
        index = np.nonzero(lon < -180)
        lon[index] = lon[index]+360.
    except :
        print "mtCoord go Boom- but caught, distrust lon limits"
        pass
    return lon,lat

def solveMT () :
    data_dir = os.environ["DESGW_DATA_DIR"]
    file = data_dir + "mcbrydethomas-psi.dat"
    fd = open(file,"w")
    for i in range(-90,91) :
        dec = i*2*np.pi/360.
        a = np.sin(dec)*(2+np.sqrt(2))/2.
        close = 10000; enuff = -1
        for j in range(-9000,9000) :
            psi = (j/100.)*2*np.pi/360.
            b = np.sin(psi) + np.sin(psi/2.)
            delta = np.abs(b-a)
            if delta < close :
                enuff = j/100.
                close = delta
                bigb = b
        print i, enuff, round(close,4), round(a,2), round(bigb,2)
        fd.write("{:d} {:6.3f}\n".format(i,enuff))
    fd.close()
def readMT() :
    global psi_spline
    data_dir = os.environ["DESGW_DATA_DIR"]
    file = data_dir + "mcbrydethomas-psi.dat"
    dec,psi = np.genfromtxt(file,unpack=True)
    psispline =interp1d(np.array(dec,dtype='float'),np.array(psi,dtype='float'),fill_value="extrapolate",bounds_error=False)
    psi_spline = psispline
    return psispline, dec, psi



