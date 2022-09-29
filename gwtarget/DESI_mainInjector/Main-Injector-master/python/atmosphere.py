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

def transmission(airmass, filter, refAirmass=1.3) :
    extCoeff = extinctionModel (filter) 
    transmission = 10**(-0.4 * extCoeff*(airmass-refAirmass))
    return transmission

def dirtTransmission(zd) :
    deg_to_rad = 2*np.pi/360.
    ix =np.nonzero(zd > 90*deg_to_rad)
    transmission = 0.0*zd+1.0
    transmission[ix] = 10.0**-100
    return transmission

def lunarDirtTransmission(moon_sep) :
    deg_to_rad = 2*np.pi/360.
    ix =np.nonzero(moon_sep < 0.5*deg_to_rad)
    transmission = 0.0*moon_sep+1.0
    transmission[ix] = 10.0**-100
    return transmission

def extinctionModel (filter) :
    if filter == "u" : k=0.58
    elif filter == "g" : k=0.18
    elif filter == "r" : k=0.09
    elif filter == "i" : k=0.08
    elif filter == "z" : k=0.08
    elif filter == "y" : k=0.08
    else : raise Exception ("no such filter {}".format(filter))
    return k

