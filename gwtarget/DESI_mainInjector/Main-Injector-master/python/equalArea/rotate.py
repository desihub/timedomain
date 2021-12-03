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

def rotate (x, y, angle) :
    angle = angle*2*np.pi/360.
    rotMatrix = np.array([[np.cos(angle), -np.sin(angle)], 
                       [np.sin(angle),  np.cos(angle)]], dtype=np.float64)
    newx = x*rotMatrix[0,0] + y*rotMatrix[0,1]
    newy = x*rotMatrix[1,0] + y*rotMatrix[1,1]
    return newx,newy
    
#
# x axis points at ra = 0
# y axis points at ra = 90
# z axis points at dec = 90
#
def sphericalToCartesian (ra, dec, r=1) :
    ra = ra*2*np.pi/360.
    #dec = dec*2*np.pi/360.
    # npd = north polar distance = co-latitude
#    npd = 90+dec
#    npd = npd*2*np.pi/360.
#    x = r * np.cos(ra)*np.sin(npd)
#    y = r * np.sin(ra)*np.sin(npd)
#    z = r * np.cos(npd)
    dec = dec*2*np.pi/360.
    x = r * np.cos(ra)*np.cos(dec)
    y = r * np.sin(ra)*np.cos(dec)
    z = r * np.sin(dec)
    return x,y,z


def cartesianToSpherical (x, y, z) :
    r = np.sqrt( x**2 + y**2 + z**2)
# the accuracy isn't good enough?
    #index =  np.nonzero(abs(r) < 1e-16)
    #index1 = np.nonzero(abs(r) >= 1e-16)
    #npd = np.zeros(len(r))
    #npd[index] = np.arccos(z)
    #npd[index1] = np.arccos(z/r)
    npd = np.arccos(z/r)
    ra  = np.arctan2(y,x)
    ra  =   ra*360/(2*np.pi)
    dec =   90- npd*360/(2*np.pi) 
    try :
        index = np.nonzero(ra < 0); ra[index] = ra[index]+360;
        index = np.nonzero(ra >= 360); ra[index] = ra[index]-360;
        index = np.nonzero(dec > 90); dec[index] = 90-(dec[index]-90);
        index = np.nonzero(dec < -90); dec[index] = -90+(dec[index]+90)
    except :
        pass
    return ra, dec, r

def rotateAboutXaxis (x, y, z, alpha, verbose = 0) :
    if verbose : print "\t x axis rotation of ", alpha, "given  ", x[0], y[0], z[0]
    alpha = alpha*2*np.pi/360.
    xp = x
    yp = y*np.cos(alpha) - z*np.sin(alpha)
    zp = y*np.sin(alpha) + z*np.cos(alpha)
    return xp,yp,zp
def rotateAboutYaxis (x, y, z, alpha, verbose = 0) :
    if verbose : print "\t y axis rotation of ", alpha, "given  ", x[0], y[0], z[0]
    alpha = alpha*2*np.pi/360.
    # correct
    xp = z*np.sin(alpha) + x*np.cos(alpha)
    yp = y
    zp = z*np.cos(alpha) - x*np.sin(alpha)
    # but at one point functional
    xp = x*np.cos(alpha) - z*np.sin(alpha)
    yp = y
    zp = x*np.sin(alpha) + z*np.cos(alpha)
    return xp,yp,zp
def rotateAboutZaxis (x, y, z, alpha, verbose = 0) :
    if verbose : print "\t z axis rotation of ", alpha, "given  ", x[0], y[0], z[0]
    alpha = alpha*2*np.pi/360.
    xp = x*np.cos(alpha) - y*np.sin(alpha)
    yp = x*np.sin(alpha) + y*np.cos(alpha)
    zp = z
    return xp,yp,zp

def getEulerAngles (ra, dec) :
    alpha = 0
    beta = dec
    gamma = ra
    return alpha, beta, gamma

def test_euler (raCen, decCen) :
    alpha, beta, gamma = getEulerAngles (raCen, decCen)
    print "rotate about axes by Euler angles: ", alpha, beta, gamma

    ra = [-30, -10, -10, 0,  0, 30]
    dec = [ 0,   1,  -1, 3, -3,  0]
    ra = np.array(ra); dec = np.array(dec)
    print "ra,dec pairs: ",
    for i in range(0,len(ra)) :
        print "{:6.2f} {:6.2f}          ".format(ra[i], dec[i]),
    print ""
    x,y,z = sphericalToCartesian(ra,dec,1)

    i = 1
    print "{:d} start              : {:10.5f}  {:10.5f}  {:10.5f}".format(i,float(x[i]),float(y[i]),float(z[i]))
    x,y,z =rotateAboutZaxis(x,y,z,alpha)
    print "Z-rot done {:5.1f} deg : {:10.5f}  {:10.5f}  {:10.5f}".format(alpha, float(x[i]),float(y[i]),float(z[i]))

    x,y,z = rotateAboutYaxis(x,y,z,beta)
    print "Y-rot done {:5.1f} deg : {:10.5f}  {:10.5f}  {:10.5f}".format(beta, float(x[i]),float(y[i]),float(z[i]))

    x,y,z = rotateAboutZaxis(x,y,z,gamma)
    print "Z-rot done {:5.1f} deg : {:10.5f}  {:10.5f}  {:10.5f}".format(gamma, float(x[i]),float(y[i]),float(z[i]))

    ra,dec,r = cartesianToSpherical (x, y, z)
    ra = np.array(ra); dec = np.array(dec)
    print "ra,dec pairs out: ",
    for i in range(0,len(ra)) :
        print "{:6.2f} {:6.2f}          ".format(ra[i], dec[i]),
    print ""

# definintions:
#   zxz rotation conventions, following the solid body
#   the ellipse is centered at the z=0, x=0, y=1 position, ra=90, dec=0
#   the target is ra,dec; our aim is to place the y-axis onto it, after the three rotations
#       line of nodes: this will be at ra+- 90 degrees
#           thus alpha = x-x' = 0-(ra-90) = 90-ra
#           then beta = rotate about x'axis by declination, 
#           and gamma = 0, as the X axis is coincident with the nodes
#           
def getEulerAngles2 (ra, dec) :
    alpha = -90
    beta = -dec
    gamma =  ra+90
    alpha = 0
    beta = dec
    gamma =  ra
    return alpha, beta, gamma

def test_euler_mk2 (raCen, decCen) :
    alpha, beta, gamma = getEulerAngles2 (raCen, decCen)
    print "rotate about axes by Euler angles: ", alpha, beta, gamma

    ra =  [60, 80, 80, 90, 90, 90.]
    dec = [ 0,  1, -1,  3, -3, 0.]
    ra = np.array(ra); dec = np.array(dec)
    ra = ra-90
    print "ra,dec pairs: ",
    for i in range(0,len(ra)) :
        print "{:6.2f} {:6.2f}          ".format(ra[i], dec[i]),
    print ""
    x,y,z = sphericalToCartesian(ra,dec,1)

    i = 1
    print "{:d} start              : {:10.5f}  {:10.5f}  {:10.5f}".format(i,float(x[i]),float(y[i]),float(z[i]))
    x,y,z =rotateAboutZaxis(x,y,z,alpha)
    print "Z-rot done {:5.1f} deg : {:10.5f}  {:10.5f}  {:10.5f}".format(alpha, float(x[i]),float(y[i]),float(z[i]))

    x,y,z = rotateAboutXaxis(x,y,z,beta)
    print "X-rot done {:5.1f} deg : {:10.5f}  {:10.5f}  {:10.5f}".format(beta, float(x[i]),float(y[i]),float(z[i]))

    x,y,z = rotateAboutZaxis(x,y,z,gamma)
    print "Z-rot done {:5.1f} deg : {:10.5f}  {:10.5f}  {:10.5f}".format(gamma, float(x[i]),float(y[i]),float(z[i]))

    ra,dec,r = cartesianToSpherical (x, y, z)
    ra = np.array(ra); dec = np.array(dec)
    print "ra,dec pairs out: ",
    for i in range(0,len(ra)) :
        print "{:6.2f} {:6.2f}          ".format(ra[i], dec[i]),
    print ""





#
# a combination of cunning and euler rotation;
#   a solid body rotation of an great circle ellipse
#   centered on ra, dec=0, rotated up to the dec
#   of interest
#
def rotateAboutEasternPole (raCen, decCen, ra, dec) :
    x, y, z = sphToCartesian(raCen, ra, dec)
    x, y, z = rotateAboutXaxis (x, y, z, decCen)
    ra, dec = cartesianToSph (raCen, x, y, z) 
    return ra,dec

# redefine ra to be 90 away from this ra
def sphToCartesian(ra0, ra, dec, r=1) :
    ra = (ra-(ra0-90))*2*np.pi/360.
    dec = dec*2*np.pi/360.
    x = r * np.cos(ra)*np.cos(dec)
    y = r * np.sin(ra)*np.cos(dec)
    z = r * np.sin(dec)
    return x,y,z

def cartesianToSph (ra0, x, y, z) :
    r = np.sqrt( x**2 + y**2 + z**2)
    npd = np.arccos(z/r)
    ra  = np.arctan2(y,x)
    ra  =   (ra0-90) + ra*360/(2*np.pi)
    dec =   90- npd*360/(2*np.pi) 
    index = np.nonzero(ra < 0); ra[index] = ra[index]+360;
    index = np.nonzero(ra >= 360); ra[index] = ra[index]-360;
    index = np.nonzero(dec > 90); dec[index] = 90-(dec[index]-90);
    index = np.nonzero(dec < -90); dec[index] = -90+(dec[index]+90)
    return ra, dec
