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

def skyFiducial(filter) :
    zenith_distance = 40. ;# degrees = 1.3 airmasses
    sky = sky_brightness_model(filter, zenith_distance) 
    return sky

# main routine
def sky_brightness_at_time(filter, zenith_distance, moon_zd, moon_sep, moon_phase) :
    import atmosphere
    extinction = atmosphere.extinctionModel(filter)
    
    sky_mag = total_sky_brightness(filter, zenith_distance, extinction, 
        moon_phase, moon_sep, moon_zd)
    return sky_mag

# first of many helper routines
def total_sky_brightness (filter, zenith_distance, extinction, phase, moon_sep, moon_zd) :
    sky = sky_brightness_model(filter, zenith_distance)
    moon = moon_brightness_model (filter, 
        extinction, phase, moon_sep, moon_zd,zenith_distance)
    # these are additive
    linear_sky =  10**(-0.4*sky  +9)
    linear_moon = 10**(-0.4*moon +9)
    total = sky + -2.5*np.log10( (linear_sky+linear_moon)/linear_sky)
    return total

# zenith_distance is in radian
def sky_brightness_model(filter, zenith_distance) :
    import atmosphere
    extinction = atmosphere.extinctionModel(filter)
    
    if   filter == "g" : sky = 22.0;  slope = 0.0114
    elif filter == "r" : sky = 20.9;  slope = 0.0131
    elif filter == "i" : sky = 19.96; slope = 0.0081
    elif filter == "z" : sky = 18.63; slope = 0.0081
    elif filter == "y" : sky = 18.0;  slope = 0.0150
    elif filter == "Y" : sky = 18.0;  slope = 0.0150
    else : raise Exception ("no such filter {}".format(filter)) 

    # those values are defined such that delta skybrightness = 0 at 1.3 airmass
    # so remove that
    sky = sky + slope*(zenith_distance-(40.*2*np.pi/360.))
    # and add it tback in
    sky = ks_sky_brightness (sky, extinction, zenith_distance)
    return sky

#
# Let us implement the sky brightness model of Krisciunas and Schaefer 1991
#
# the increase of sky brightness with zenith distance is 
# independent of color, though extinction is color-dependent
# returns magnitude
def ks_sky_brightness (sky, k, zenith_distance) :
    airmass = scattering_airmass(zenith_distance)
    sky_zd = -2.5*np.log10( airmass) + k*(airmass - 1)
    mag = sky + sky_zd
    return mag

#
# moony sky brightness
# K&S eq 15
def moon_brightness_model ( filter, extinction, phase, moon_sep, moon_zd, obj_zd) :
    test = 0
    lunar_airmass = scattering_airmass(moon_zd)
    obj_airmass   = scattering_airmass(obj_zd)
    if test : 
        print "\t airmass of moon & object : ", lunar_airmass, obj_airmass
    lunar_airmass = extinction*lunar_airmass
    obj_airmass = -2.5*np.log10( 1-10**(-0.4*extinction*obj_airmass))

    scattering = scattering_function(moon_sep)
    if moon_zd > 93 : phase = 0.0 
    lunar_mag = moon_mag(phase, filter)
    sky = scattering + lunar_mag + lunar_airmass + obj_airmass
    return sky

# K&S equation 3
def scattering_airmass(zenith_distance) :
    rad_per_deg=(2*np.pi/360.)
    sin_zd = np.sin( zenith_distance )
    airmass = 1/(np.sqrt( 1 - 0.96*sin_zd**2 ))
    # the function is symmetrical about zd=90, not what we want
    if airmass.size > 1 :
        ix=np.nonzero(zenith_distance > 90*rad_per_deg)
        if ix[0].size > 0 :
            airmass[ix] = airmass.max()
    return airmass

# K&S eq 21
# moon_sep in radians
def scattering_function (moon_sep) :
    rad_per_deg=(2*np.pi/360.)
    term_1a = 10**5.36
    term_1b = 1.06 + np.cos(moon_sep )**2
    term_2 = 10**(6.15-moon_sep/(40.*rad_per_deg))
    scattering = -2.5*np.log10(term_1a * term_1b + term_2)
    return scattering

# K&S eq 20
# phase in %ilumination ala skycalc
#   but note it goes from 0 to 1
def moon_mag (phase, filter) :
    vmag = 3.84 + 0.026*phase + 4.0e-9 * phase**4

    # vmag is calculated as -2.5 * log(lumens)
    vmag  = 10**(-0.4*vmag)

# B = 34.08 exp(20.7233 - 0.92104 V)
# ln B = ln 34.08 + 20.7233 - 0.92104V
# 0.92104 V = ln 34.08 + 20.7233 - ln B
# V = (1/0.92104) * (lnv 34.08 + 20.7233 - ln B )
# V = 1.0958 * (3.52871 + 20.7233 - ln B)
# V = 26.5754 - 1.0958 ln B
    vmag = 26.5754 - 1.0958*np.log(vmag) ;# from K&S eq 1

    # assume moon colors are those of the sun
    delgV = 0.34 ;# ~ from sdss.org filter transforms page 
    delgr = 0.44 ;# ~ from sdss.org filter transforms page 
    delri = 0.11 ;# ~ from sdss.org filter transforms page 
    deliz = 0.03 ;# ~ from sdss.org filter transforms page 
    delzy = 0.00 ;# assume
    if   filter == "V" : mag = 0
    elif filter == "g" : mag = delgV
    elif filter == "r" : mag = delgV - delgr
    elif filter == "i" : mag = delgV - delgr - delri
    elif filter == "z" : mag = delgV - delgr - delri - deliz
    elif filter == "y" : mag = delgV - delgr - delri - deliz - delzy
    else : raise Exception("no such filter {}".format(filter))
    mag = mag+vmag
    return mag

