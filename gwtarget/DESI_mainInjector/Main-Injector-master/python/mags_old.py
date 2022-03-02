import numpy as np
import healpy as hp
import hp2np
import os
from pyslalib import slalib
import atmosphere  
import telescope
import dustModel  
import seeingModel
import skyModel
from equalArea import mcbryde

license="""
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style. (Note the rate at which people who stand on the
   summit of K2 successfully make it down.)

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

#
# get pyslalib, python over industrial grade spherical astronomy code slalib
# https://github.com/scottransom/pyslalib
# % wget  https://github.com/scottransom/pyslalib/archive/master.zip
# % unzip master.zip
# % cd pyslalib-master
# % make
# % python setup.py install --home=$WORK_DIR/python-home-stash/
# % python test/test_slalib.py
# % PYTHONPATH=$PYTHONPATH:$WORK_DIR/python-home-stash/ ;export PYTHONPATH

class observed(object):
    """
    Constuct the observed map.
    This produces the limiting magnitude map for a 30 sec i-band
    exposure using DECam/Blanco over the whole sky given a time in MJD.

    Internally all coordinates, time and space, are in radians, except MJD.
    """
    def __init__(self, 
            ra, dec, values, mjd, alpha=0., 
            degradeRes=True, # change map resolution
            doMaps=True,  # don't think of ra,dec,vals as a healpy map; don't do map work
            camera="decam",
            verbose=True   # be very loud or no
        ) :
        self.verbose      = verbose
        data_dir          = os.environ["DESGW_DATA_DIR"]
        ctio_lat          = -30.16527778
        ctio_lon          = -70.8125
        ctio_height       = 2215.
        kpno_lat          = 31.9600784
        kpno_lon          = -111.598169
        kpno_height       =  2067.
        if camera == "decam" :
            obs_lat = ctio_lat
            obs_lon = ctio_lon
            obs_height = ctio_height
        elif camera == "desi" :
            obs_lat = kpno_lat
            obs_lon = kpno_lon
            obs_height = kpno_height

        self.date         = slalib.sla_djcl(mjd)

        self.degToRad     = 0.0174532925199
        self.lat          = obs_lat*self.degToRad
        self.lon          = obs_lon*self.degToRad
        self.height       = obs_height

        self.mjd          = mjd
        self.ra           = ra*self.degToRad
        self.dec          = dec*self.degToRad
        self.map          = values
        if doMaps :
            self.nside        = hp.npix2nside(values.size)
        else :
            self.nside        = 0

        # observational astronomy
        self.lst          = self.mjdToLST(self.mjd, self.lon)
        self.ha,self.zd   = self.equatorialToObservational (self.ra, self.dec, self.lst, self.lat) 
        self.airmass      = self.airmassModel(self.zd) 
        self.sunData      = self.getSolarPosition (self.mjd, self.lon, self.lat, self.lst)
        self.sunZD        = self.sunData[2]
        self.moonData     = self.getLunarPosition (self.mjd, self.lon, self.lat, self.lst)
        self.moonZD       = self.moonData[2]
        self.moonSep      = self.getLunarSeparation(self.ra, self.dec, self.moonData)
        self.moonPhase    = self.getLunarPhase()
        self.moonRa       = self.moonData[0]
        self.moonDec      = self.moonData[1]

        # telescope limits
        self.limits       = 0.0

        # galactic dust
        dust_dir          = data_dir
        dust_file         = "plank-ebv-HFI_CompMap_ThermaDustModel.fits"
        ra,dec,ebv        = dustModel.loadDust(dust_dir, dust_file)
        if degradeRes :
            ra,dec,ebv    = hp2np.map2np (ebv, resolution=self.nside, fluxConservation=False)
        if not doMaps :
            ra      = self.ra/self.degToRad
            dec     = self.dec/self.degToRad
            ebv = hp.get_interp_val(ebv, ra,dec, lonlat=True)
            ra      = self.ra
            dec     = self.dec
        self.dust_ra      = ra
        self.dust_dec     = dec
        self.ebv          = ebv

        # stellar density in the form of a probability of recognizing object map
        print("\t loading inverse stellar density = probability of recognition map")
        ra,dec,precog     = hp2np.hp2np(data_dir + "precognize.fits")
        if degradeRes :
            ra,dec,precog = hp2np.map2np (precog, resolution=self.nside, fluxConservation=False)
        if not doMaps :
            ra      = self.ra/self.degToRad
            dec     = self.dec/self.degToRad
            precog = hp.get_interp_val(precog, ra,dec, lonlat=True)
            ra      = self.ra
            dec     = self.dec
        self.pra          = ra
        self.pdec         = dec
        self.precog       = precog

        # limiting mag
        self.maglim       = ""
        self.maglimall    = ""

        # where to put the central meridan in the  mcbryde-thomas projection
        # - the ra which one wants on the meridian
        self.mcbryde_alpha = alpha
        # set the figure axes so that the aspect ratio is set by data
        # for the right Mcbryde plot, which projects onto a flat x,y page
        # plt.axes().set_aspect('equal')

    # change to a new time
    def resetTime (self, mjd) :
        self.mjd          = mjd
        self.lst          = self.mjdToLST(self.mjd, self.lon)
        self.ha,self.zd   = self.equatorialToObservational (self.ra, self.dec, self.lst, self.lat) 
        self.airmass      = self.airmassModel(self.zd) 
        self.sunData      = self.getSolarPosition (self.mjd, self.lon, self.lat, self.lst)
        self.sunZD        = self.sunData[2]
        self.moonData     = self.getLunarPosition (self.mjd, self.lon, self.lat, self.lst)
        self.moonZD       = self.moonData[2]
        self.moonSep      = self.getLunarSeparation(self.ra, self.dec, self.moonData)
        self.moonRa       = self.moonData[0]
        self.moonDec      = self.moonData[1]

    def limitMag (self, filter, exposure=30, verbose=True) :
        #print "\t JTA limiting magnitude", self.mjd, self.mjdToLST(self.mjd, self.lon)
        #print self.ha[int(self.ha.size/2.)]
        mjd       = self.mjd
        zd        = self.zd
        sun_zd    = self.sunZD
        moon_zd   = self.moonZD
        moon_phase= self.moonPhase
        moon_sep  = self.moonSep
        airmass   = self.airmass
        ebv       = self.ebv
        alpha     = self.mcbryde_alpha
        sunIsUp = self.sunBrightnessModel (sun_zd )
        if sunIsUp :
            m = np.copy(ebv)*0.1+-10.0
            mglobal = np.copy(ebv)*0.1+-10.0
            print("\t ... the sun is up")
        else :
            dust    = self.dustTransmission(filter, ebv)
            atmo    = self.atmosphereTransmission(zd, airmass, filter, moon_sep)
            seeing  = self.seeing(airmass, filter, seeingAtZenith=0.9)
            sky     = self.skyBrightness(zd, moon_zd, moon_sep, moon_phase, filter) 
            skyFid  = self.skyBrightnessFid(filter)

            telescope = self.telescopeLimits(self.ha, self.dec)
            self.limits =telescope

            if filter == "g" : m_zp = 23.3
            elif filter == "r" :m_zp = 23.4
            elif filter == "i" :m_zp = 22.9
            elif filter == "z" :m_zp = 22.5
            elif filter == "y" :m_zp = 20.6
            else : raise Exception("no such filter")
            
            SN = telescope*dust*atmo*(1./seeing)*np.sqrt(exposure/30.)
            ix = np.nonzero( SN <= 0)
            SN[ix] = 1.0e-12
            m = m_zp + 2.5*np.log10( SN ) + 0.5*(sky - skyFid)
            # arbitarily limit limiting mag to that of moon
            ix = np.nonzero(m < 0)
            m[ix] =  0
            ix, = np.where(telescope == 0)
            m[ix] = -11.0

            SNglobal = dust*atmo*(1./seeing)*np.sqrt(exposure/30.)
            ix = np.nonzero( SNglobal <= 0)
            SNglobal[ix] = 1.0e-12
            mglobal = m_zp + 2.5*np.log10( SNglobal ) + 0.5*(sky - skyFid)
            # arbitarily limit limiting mag to that of moon
            ix = np.nonzero(mglobal < 0 )
            mglobal[ix] = 0


        # project into equal area map
        x,y = mcbryde.mcbryde(self.ra/self.degToRad, self.dec/self.degToRad, alpha=alpha)
        self.x = x
        self.y = y
        hx,hy = mcbryde.mcbryde(self.ha/self.degToRad, self.dec/self.degToRad, alpha=alpha)
        self.hx = hx
        self.hy = hy
        self.maglim = m
        self.maglimall = mglobal

    def dustTransmission(self, filter, ebv) :
        if self.verbose: print("\t ... dust")
        dustTransmission = dustModel.dustTransmission(filter, ebv)
        return dustTransmission

    def atmosphereTransmission(self, zd, airmass, filter, moon_sep, refAirmass=1.3) :
        if self.verbose: print("\t ... atmosphere")
        atransmission = atmosphere.transmission(airmass, filter, refAirmass)
        if self.verbose: print("\t ... earth")
        dtransmission = atmosphere.dirtTransmission(zd)
        if self.verbose: print("\t ... moon")
        mtransmission = atmosphere.lunarDirtTransmission(moon_sep)
        transmission = atransmission*dtransmission*mtransmission
        return transmission

    def telescopeLimits (self, ha, dec) :
        if self.verbose: print("\t ... telescope limits")
        limits = telescope.blancoLimits(ha, dec)
        return limits 

    def seeing(self, airmass, wavelength=775., seeingAtZenith=1.0) :
        if self.verbose: print("\t ... seeing")
        seeing = seeingModel.seeingWithAirmassAndLambda(airmass, wavelength, seeingAtZenith)
        return seeing

    def skyBrightness(self, zd, moon_zd, moon_sep, moon_phase, filter) :
        if self.verbose: print("\t ... sky brightness")
        sky = skyModel.sky_brightness_at_time( filter, zd, moon_zd, moon_sep, moon_phase) 
        return sky

    def skyBrightnessFid(self, filter) :
        skyfiducial = skyModel.skyFiducial(filter)
        return skyfiducial

    #
    # First, let's turn RA Dec into local sideral time, the start of all observations
    #
    # equation of Equinoxes is an ~1 second effect
    # "east longitude", where positive numbers increase as one moves east
    #   recall that lat, long is in radians
    def mjdToLST (self, mjd, eastLongitude) :
        if self.verbose: print("\t MJD to LST")
        gmst        = slalib.sla_gmst(mjd)
        eqEquinoxes = slalib.sla_eqeqx(mjd)
        lst         = gmst + eqEquinoxes + eastLongitude
        return lst 

    #
    # Now we turn LST and dec into Hour Angle and zenith distance
    #
    # Technically this is apparent ha,dec; i.e. we ignore diurnal abberation
    #   recall that lat, long is in radians
    # The best way to view this output is as hexbin(ra,dec,zd) or ha
    def equatorialToObservational (self, ra, dec, lst, latitude) :
        if self.verbose: print("\t LST to HA,zenith distance")
        ha = lst - ra
        zd = self.zenithDistance(ha, dec, latitude)
        ix = np.nonzero(ha > 180.*2*np.pi/360.)
        ha[ix] = ha[ix] - 360*2*np.pi/360.
        return ha,zd

    #
    # Calculating zenith distance
    #
    # sin(HC) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(LHA)
    def zenithDistance(self, ha, dec, latitude) :
        degToRad = 2.*np.pi/360.
        sinAltRad = np.sin(latitude)*np.sin(dec) + np.cos(latitude)*np.cos(dec)*np.cos(ha)
        altRad = np.arcsin(sinAltRad)
        zenithDist = 90*degToRad - altRad
        return zenithDist

    #
    # Calculating airmass
    #
    # from Young 1994 "Air mass and refraction", Applied Optics 33:1108-1110
    # max error = 0.0037 airmasses at zd=90 degrees. (!)
    def airmassModel(self, zd) :
        if self.verbose: print("\t zenith distance to airmass")
        airmass = zd*0+45 # a hack to set the values of airmass under the horizon
        ix = np.nonzero(zd <= 89.9999*2*np.pi/360.)
        coszd = np.cos(zd)
        coszdsq = coszd*coszd
        numerator = 1.002432*coszdsq + 0.148386*coszd + 0.0096467
        denominator = coszdsq*coszd + 0.149864*coszdsq + 0.0102963*coszd + 0.000303978
        airmass[ix] = numerator[ix]/denominator[ix]
        return airmass
        

    #
    # check lunar position
    #
    def getLunarPosition (self, mjd, eastLongitude, latitude, lst ) :
        ra, dec, diam = slalib.sla_rdplan(mjd, 3, eastLongitude, latitude)
        ha = lst - ra
        zd = self.zenithDistance(ha, dec, latitude)
        return ra,dec,zd

    #
    # returns moon phase in degrees 
        # K&S want moon phase as angle in degrees, 
        #   where 0 = full, 90 equals half, and  180 = new
    #
    def getLunarPhase (self) :
        moon_ra, moon_dec = self.moonData[0], self.moonData[1]
        sun_ra, sun_dec   = self.sunData[0], self.sunData[1]
        # moon is full when elongation = 180, new when elongation = 0
        moon_elongation = self.gc_separation(sun_ra, sun_dec, moon_ra, moon_dec)
        # K&S want moon phase as angle in degrees, 
        #   where 0 = full, 90 equals half, and  180 = new
        phase = (180*(2*np.pi/360.) - moon_elongation)*360./2/np.pi
        return phase

    def getLunarSeparation(self, ra, dec, moonData) :
        moon_ra, moon_dec = self.moonData[0], self.moonData[1]
        moon_sep = self.gc_separation(ra, dec, moon_ra, moon_dec)
        return moon_sep

    #
    # check solar position
    #
    def getSolarPosition (self, mjd, eastLongitude, latitude, lst ) :
        ra, dec, diam = slalib.sla_rdplan(mjd, 0, eastLongitude, latitude)
        ha = lst - ra
        zd = self.zenithDistance(ha, dec, latitude)
        return ra,dec,zd
    #
    # check to see if sun near the horizen
    #
    def sunBrightnessModel (self, sunZD) :
        twilight = 100.*2*np.pi/360. 
        if sunZD <= twilight :
            bright = True
        else :
            bright = False
        return bright
    
    #
    # great circle separations
    #
    def gc_separation(self, ra1, dec1, ra2, dec2) :
        delDec = dec1-dec2
        delRa = ra1-ra2
        dhav = self.haversine(delDec)
        rhav = self.haversine(delRa)
        hav = dhav + np.cos(dec1)*np.cos(dec2)*rhav
        gc_distance = self.ahaversine(hav)
        return gc_distance

    def haversine(self, theta) :
        hav = np.sin(theta/2.)**2
        return hav
    def ahaversine(self, x) :
        ahav = 2*np.arcsin(np.sqrt(x))
        return ahav

def findNightDuration(mjd, camera="decam") :
    ctio_lat          = -30.16527778
    ctio_lon          = -70.8125
    ctio_height       = 2215.
    kpno_lat          = 31.9600784
    kpno_lon          = -111.598169
    kpno_height       =  2067.
    if camera == "decam" :
        obs_lat = ctio_lat
        obs_lon = ctio_lon
        obs_height = ctio_height
    elif camera == "desi" :
        obs_lat = kpno_lat
        obs_lon = kpno_lon
        obs_height = kpno_height

    degToRad = 2.*np.pi/360.
    lat          = obs_lat*degToRad
    lon          = obs_lon*degToRad
    height       = obs_height

    imjd = np.int(mjd)
    start_mjd = imjd - 6./24.  ;# before sunset at CTIO

    sunset = ""
    sunrise = ""
    # check every minute
    for i in np.arange(0, 1., 1./(24.*60) ) :
        mjd = start_mjd + i
        gmst        = slalib.sla_gmst(mjd)
        eqEquinoxes = slalib.sla_eqeqx(mjd)
        lst         = gmst + eqEquinoxes + lon

        sunra, sundec, diam = slalib.sla_rdplan(mjd, 0, lon, lat)
        sunha = lst - sunra
        sinAltRad = np.sin(lat)*np.sin(sundec) + \
            np.cos(lat)*np.cos(sundec)*np.cos(sunha)
        altRad = np.arcsin(sinAltRad)
        zenithDist = 90*degToRad - altRad

        twilight = 100.*2*np.pi/360. 
        if zenithDist <= twilight :
            bright = True
        else :
            bright = False
        if sunset == "" and bright == False :
            sunset = mjd
        if sunset != "" and sunrise == "" and bright == True :
            sunrise = mjd
    duration = sunrise-sunset
    return duration, sunset, sunrise


