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

# airmass seeing dependence &
#  komologorov seeing wavelength dep
def seeingWithAirmassAndLambda(airmass, filter, seeingAtZenith=1.0) :
    
    seeing = seeingAtZenith*airmass**(3./5.)

    refWavelength = 775 # i-band
    wavelength=filterEffWavelength(filter)
    komol = (wavelength/refWavelength)**-0.2
    seeing = seeing*komol
    return seeing

def filterEffWavelength(filter) :
    if filter == "u" : effLam = 380.
    elif filter == "g" : effLam = 475.
    elif filter == "r" : effLam = 635.
    elif filter == "i" : effLam = 775.
    elif filter == "z" : effLam = 925.
    elif filter == "y" : effLam = 1000.
    else : raise Exception ("no such filter {}".format(filter))
    return effLam

#
