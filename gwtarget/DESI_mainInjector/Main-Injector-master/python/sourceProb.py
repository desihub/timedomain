import numpy as np
import os
import scipy.stats
import mags
import hp2np
import warnings

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

class map(object):
    """
    """
    def __init__(self, observation, type="bright", apparent_mag_source=21.5) :
        """
        """
        data_dir = os.environ["DESGW_DATA_DIR"] 

        self.limitingMag = observation.maglim
        self.limits      = observation.limits
        tryApparent = True
        if type == "bright" :
            if tryApparent :
                self.lumModel = "apparent"
                self.modelAbsoluteMagnitude = apparent_mag_source
                self.modelAbsoluteMagnitudeSigma = 0
            else : 
                # this can be used to scale the magnitudes from the source model
                # this can be done as the source abs magnitude for calculations
                # is taken from self.absMagMean, self.absMagSigma
                # whereas we're setting the absolute scale here with
                # self.modelAbsoluteMagnitude, self.modelAbsoluteMagnitudeSigma
                # Currently this is done in mapsAtTimeT.probabilityMaps
    
                self.lumModel = "kn-gaussian"
                # LIGO O1 and O2
                #self.modelAbsoluteMagnitude = -11.1
                # LIGO O3
                self.modelAbsoluteMagnitude = -15.5
                self.modelAbsoluteMagnitudeSigma = 1.0
        elif type == "dark" :
            if tryApparent :
                self.lumModel = "apparent"
                self.modelAbsoluteMagnitude = 21.5
                self.modelAbsoluteMagnitudeSigma = 0
            else : 
                # fixed luminosity
                self.lumModel = "bh-gaussian"
                self.modelAbsoluteMagnitude = -25.0
                self.modelAbsoluteMagnitudeSigma = 1.0
        else :
            raise Exception (
                "only trigger types known are bright and dark, not {}".format(type))
        self.absMagMean = self.modelAbsoluteMagnitude 
        self.absMagSigma = self.modelAbsoluteMagnitudeSigma 

        # get the P_recognition
        self.pra = observation.pra
        self.pdec = observation.pdec
        self.precog = observation.precog

        # keep the answer
        self.probMap = ""

        # for plotting contours, keep the intermediate data around
        self.xi, self.yi, self.zi = ["","",""]
        self.lastCtype = ""
        

    def calculateProb(self, ligo, ligo_distance, ligo_distance_sigma, verbose = True) :
        import scipy.integrate
        warnings.filterwarnings("ignore")
        # bookkeeping for plotting
        self.zi= ""

        #
        telescopeLimits = self.limits
        limitingMag = self.limitingMag 

        # check if sun is up
        if limitingMag.mean() < -9 :
            self.probMap = np.zeros(ligo.size)
            return 0

        # realwork
        ligo_spatial = ligo
        ligo_d_mean = ligo_distance
        ligo_d_var = ligo_distance_sigma**2

# we are going to change variables to distance modulus (magnitude),
# then, assume error is gaussian in DM and in absolute magnitude
# so that we can add guassians to get a gaussian.
        absMag_mean = self.absMagMean
        absMag_var = self.absMagSigma**2

        test_sum = ligo_d_var.sum()
        if self.lumModel == "apparent" :
            # implementing this in Feb 2020
            prob_map = np.zeros(ligo_spatial.size)
            apparent_mag = absMag_mean
            ix, = np.where(apparent_mag < limitingMag)
            prob_map[ix] = 1.0
            ix, = np.where((apparent_mag >= limitingMag) & (apparent_mag < limitingMag+0.5))
            prob_map[ix] = 0.5

        else :
# we can't afford to do every pixel. 
# for res= 128, a cut at > 1e-8 brings the number of pixels
# down from 196,608 to 10^4
            ix = (ligo_spatial > 1e-8) & (ligo_d_mean > 0) & (limitingMag > 0)
            distance_mod = 5*np.log10(ligo_d_mean[ix]*1e6/10.)
            # error propgation on converting from distance to distance modulus
            dm_var = ligo_d_var[ix]*(5./np.log(10)/ligo_d_mean[ix])
            # assume gaussians, so adding gaussians to get a gaussian
            ap_mag_mean = distance_mod + absMag_mean
            ap_mag_var = absMag_var + dm_var
    
            ic = 0
            prob_map = np.copy(ligo_spatial)*0.0
            for pix  in np.nonzero(ix)[0] :
                mag = ap_mag_mean[ic]
                var = ap_mag_var[ic]
# Now we weight by r^2 dr , thus assuming equal probability per unit volume
# We will also drop constants, equivalent to assuming we are going
# to renormalize by integrating from 0 to infinity
            # normalize
                norm_min, norm_max = 7.0, 25.0 # magnitudes
                ans = scipy.integrate.quad(self.probability_density, 
                    norm_min, norm_max, args=(mag, var ) )
                norm = ans[0]; norm_error = ans[1]

                # calculate prob
                mag_min, mag_max = 7.0, limitingMag[pix]
                ans = scipy.integrate.quad(self.probability_density, 
                    mag_min, mag_max, args=(mag, var ) )
                prob = ans[0]; prob_error = ans[1]
                prob = prob/norm 
                prob_map[pix] = prob
                ic += 1
            print "\t probMap: made for source absMag {:.1f}".format(absMag_mean)

# probability of recognition propto star density
        prob_map = prob_map * self.precog
# finally, eliminate every where the telescope cannot point
        prob_map = prob_map * telescopeLimits

        self.probMap = prob_map
        return 1 

#
# This is a gaussian in apparent magnitude, weighted by r^2 (itself
# trasformed into a distance modulus (i.e. magnitude)
#
# Now we weight by r^2 dr , thus assuming equal probability per unit volume
# probability_density = np.exp( (m - ap_mag_mean)/2/ap_mag_var) * \
#                       np.exp(0.6*np.log(10.)*(m - ap_mag_mean)) 
    def probability_density(self,m, mag_mean, mag_var ) :
        #print m, mag_mean, mag_var
        pd= np.exp( -(m-mag_mean)/2/mag_var) * \
                       np.exp(0.6*np.log(10.)*(m -mag_mean))
        return pd




    # one can choose to plot the ligo contours
    # or the ligo*obs_prob contours
    #   type="ligo"   type="ls"
    #           houranlge chooses HA instead of RA projected into mcbryde coords
    def plotLigoContours(self, obs, type="ligo", whiteLine="False", hourangle=False) :
        import matplotlib
        import matplotlib.pyplot as plt
        con_levels=10
        if self.zi == "" or self.lastCtype != type:
            print "\t calculating contours for type = ",type
            if hourangle == False :
                xmin = obs.x.min(); xmax = obs.x.max()
                ymin = obs.y.min(); ymax = obs.y.max()
            else :
                xmin = obs.hx.min(); xmax = obs.hx.max()
                ymin = obs.hy.min(); ymax = obs.y.max()
            xi=np.linspace(xmin, xmax, 500)
            yi=np.linspace(ymin, ymax, 500)
            if type == "ligo" :
                probMap = obs.map
            if type == "ls" :
                probMap = obs.map*self.probMap
            if hourangle == False :
                x = obs.x; y = obs.y
            else :
                x = obs.hx; y = obs.hy
            zi = matplotlib.mlab.griddata(x,y,probMap,xi,yi)
            self.xi, self.yi, self.zi = xi,yi,zi
            self.lastCtype = type
        if not whiteLine :
            plt.contour(self.xi,self.yi,self.zi,con_levels,linewidths=3,colors="k")
        plt.contour(self.xi,self.yi,self.zi,con_levels,linewidths=0.66,colors="w")
        
