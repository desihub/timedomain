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

def dustTransmission(filter, ebv) :
    A = dustA(filter)
    dustTransmission = 10**(-0.4*A*(ebv - ebv.min()))
    return dustTransmission

def dustA(filter) :
    if filter == "u" : A = 4.544
    elif filter == "g" : A = 3.764
    elif filter == "r" : A = 2.765
    elif filter == "i" : A =  2.483
    elif filter == "z" : A =  1.935
    elif filter == "y" : A =  1.515
    else :
        raise Exception ("no such filter {}".format(filter))
    return A

def loadDust(dir = "/home/s1/annis/daedalean/desgw-map/data/", 
        file = "plank-ebv-HFI_CompMap_ThermaDustModel.fits") :
    print "\t loading dust map ",dir+file
    hpFile = True
    if not hpFile :
        import pyfits
        f=pyfits.open(dir+file)
        t = f[1].data;
        ra = t.field("ra")
        dec = t.field("dec")
        ebv = t.field("ebv")
    else :
        import hp2np
        ra,dec,ebv= hp2np.hp2np(dir+file)

    ixra=np.nonzero(ra > 180); ra[ixra]=ra[ixra]-360.
    ra = ra*2*np.pi/360.
    dec = dec*2*np.pi/360.
    return ra, dec, ebv

