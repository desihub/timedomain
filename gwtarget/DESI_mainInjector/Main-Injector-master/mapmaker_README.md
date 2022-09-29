Code to calculate the probability of seeing
a EM counterpart to a LIGO event starting from
the LIGO probability maps.

At present, calculates the limiting magnitude
of a 30 second i-band Blanco image, all sky.

Try:
./setup.sh
>>> import matplotlib.pyplot as plt; plt.ion() ; figure = plt.figure(1)
>>> import numpy as np; import os
>>> import mags

# get ligo map
>>> import hp2np
>>> import healpy as hp
>> skymap = "lalinference.fits.gz"; resolution = 128
>>> ra,dec,ligo = hp2np.hp2np(skymap)
>>> j,j,ligo_dist = hp2np.hp2np(skymap, degrade=resolution, field=1)
>>> j,j,ligo_dist_var =hp2np.hp2np(skymap, degrade=resolution, field=2)
>>> j,j,ligo_dist_norm = hp2np.hp2np(skymap, degrade=resolution, field=3)

# compute limiting mag
>>> mjd = 55435.60552; obs = mags.observed(ra,dec,ligo, mjd); omags = obs.limitMag("i")
# since this turns out to be during the day, change the time.
>>> obs.resetTime(mjd+(10./24.)); obs.limitMag("i")
# perhaps make the exposures longer than 30 sec, say 90 seconds
>>> obs.limitMag("i",exposure=90.)

# plot as ra,dec map (all coordinates are stored as radians)
>>> plt.clf(); plt.hexbin( obs.ra*360./2/np.pi, obs.dec*360./2/np.pi, obs.maglim, vmin=15); plt.colorbar()

# plot in an equal area projection- x,y are projected into a mcbryde-thomas projection
>>> plt.axes().set_aspect('equal')
>>> plt.clf(); plt.hexbin( obs.x, obs.y, obs.maglim, vmin=15); plt.colorbar()

# work with the source probabilty
>>> import sourceProb
>>> sm=sourceProb.map(obs,type="NS");  sm.calculateProb(ligo, ligo_dist, ligo_dist_sigma)
>>> plt.clf();plt.hexbin(obs.x,obs.y,sm.probMap); plt.colorbar()
>>> plt.clf();plt.hexbin(obs.x,obs.y,sm.probMap*ligo); plt.colorbar()

#  show the ligo probability contours on the limiting magnitude map
>>> plt.hexbin( obs.x, obs.y, obs.maglim, vmin=15); plt.colorbar()
>>> sm.plotLigoContours(obs)
#  show the ligo probability contours on the source detection probability map
>>> plt.hexbin( obs.x, obs.y, sm.probMap); plt.colorbar()
>>> sm.plotLigoContours(obs)
#  show the ligo*source probability contours on the limiting magnitude map
>>> plt.hexbin( obs.x, obs.y, obs.maglim, vmin=17); plt.colorbar()
>>> sm.plotLigoContours(obs, type="ls")

# do it better
>>> from equalArea import mcplot
>>> plt.clf(); mcplot.plot(obs.x,obs.y,obs.maglim,300,vmin=17); sm.plotLigoContours(obs,type="ls")
>>> plt.savefig("thefig.pdf", bbox_inches="tight")

# get the limiting distance map, assuming a fixed absolute magnitude
>>> sm=sourceProb.map(obs, type="NS");  sm.calculateProb()
>>> plt.clf(); plt.hexbin( obs.x, obs.y, sm.prob); plt.colorbar()

# lower the resolution of the maps for speed
>>> ra,dec,ligo = hp2np.hp2np(os.environ['DESGW_DATA_DIR'] + "/bayestar.fits.gz")
>>> ra2,dec2,ligo2 =hp2np.map2np(ligo,resolution=32)
>>> obs = mags.observed(ra2,dec2,ligo2 mjd, degradeRes=True);
>>> obs.resetTime(mjd+(10./24.));
>>> obs.limitMag("i", exposure=180);
>>> sm=sourceProb.map(obs, type="NS");sm.calculateProb();

# plot using hour angle instead of RA- area always in center
>>> mcplot.plot(obs.hx,obs.hy,sm.probMap,cmap="gray"); 
>>> sm.plotLigoContours(obs,type="ligo",hourangle=True);

More to come.
Jim Annis
Dec 5, 2014

=====================================================
# Now let us industrialize this.

# First for a simulation, sim 13681
>>> import getHexObservations; import mapsAtTimeT; import allMaps; import modelRead
#  === we have to get the simulation information
>>> simNum = 13681
>>> data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims/"
>>> sims, mjds, distances, ligoMap = allMaps.selectSingleSim(simNum, data_dir)
#  === let's reset the time  of the burst to make it useable
>>> mjds[0]=57295. - 0.5
#  === now we are ready to do our standard work
>>> ligo = hp.read_map(ligoMap)
>>> ra,dec,ligo = hp2np.map2np(ligo,256, fluxConservation=True)
>>> obs = mags.observed(ra,dec,ligo, mjds[0])
>>> obs.limitMag("i",exposure=180.)
# ==== get the neutron star explosion models
>>> models = modelRead.getModels()
# ==== calculate maps during a full night of observing
>>> probs,times = mapsAtTimeT.oneDayOfTotalProbability(obs,mjds[0],distances[0],models)
>>> mapsAtTimeT.probabilityMapSaver (obs, sims[0], mjds[0], distances[0], models, times, probs, "jack/")
# ==== figure out what to observe during that night
>>> n_slots, map_zero = getHexObservations.contemplateTheDivisionsOfTime(
    probs, times, hoursAvailable=6)
>>> best_slot = getHexObservations.now( n_slots, \
         mapDirectory="jack/", simNumber=13681, mapZero=map_zero)
# ==== show what we learned
>>> plt.clf();getHexObservations.observingPlot(figure,13681,5,"jack/",12)
#  recall that "figure" is from "figure = plt.figure(1)"


# And let's try an example that works for Dylan Brout's MainInjector
# have a devel version of maininjector around, say at:
# /data/des41.a/data/annis/maininjector
# === the recycler object knows about skymap, mjd, trigger_id, and outfolder
# === it has exposure_length as a config variable and distance from the event
>>> import numpy as np; import os; import getHexObservations
>>> out_dir = "jack/";  data_dir = out_dir
>>> web_dir = "jack/web/"
>>> skymap = "/data/des41.a/data/desgw/maininjector/test-triggers/M168556/M168556_bayestar.fits.gz"
>>> trigger_id = "/M168556" ; trigger_mjd = 55472.8245; distance = 91.24
>>> hours_available=6.; exposure_length=180.; 
# ==== calculate maps during a full night of observing
>>> probs,times = getHexObservations.prepare(skymap, trigger_mjd, trigger_id, \
    data_dir, exposure_length=exposure_length, distance=distance)
>>> n_slots, map_zero = getHexObservations.contemplateTheDivisionsOfTime(
    probs, times, hoursAvailable=hours_available); 
# ==== figure out what to observe during that night
>>> best_slot  = getHexObservations.now( \
    n_slots, mapDirectory=data_dir, simNumber=trigger_id, mapZero=map_zero)
>>> getHexObservations.makeObservingPlots( 
            n_slots, trigger_id, best_slot, data_dir)

# ======== debug hexalation
>>> obs, trigger_id, mjd, distance, models, times, probs,data_dir = \
         getHexObservations.prepare(skymap, trigger_mjd, trigger_id, \
        data_dir, exposure_length=exposure_length, distance=distance)
>>> gw_data_dir = os.environ["DESGW_DATA_DIR"];hexFile = gw_data_dir + "all-sky-hexCenters.txt"
>>> print probs,times
>>> obs,sm = mapsAtTimeT.probabilityMaps( obs, mjd, times[10], distance, models)
>>> raHexen, decHexen, hexVals, rank = hexalate.cutAndHexalate (\
        obs, sm, hexFile)

# ======== real event
>>> import numpy as np; import os; import getHexObservations
>>> out_dir = "jack2/";  data_dir = out_dir
>>> skymap = "/data/des41.a/data/desgw/maininjector_devel/real-triggers/G184098/G184098_bayestar.fits"
>>> skymap = "/data/des41.a/data/desgw/maininjector_devel/real-triggers/G184098/skyprobcc_cWB_complete.fits"
>>> trigger_id = "G184098" ; trigger_mjd = 57279.4102; distance = 10.
>>> hours_available=3.; exposure_length=180.; 
# ==== calculate maps during a full night of observing
>>> probs,times = getHexObservations.prepare(skymap, trigger_mjd, trigger_id, \
    data_dir, exposure_length=exposure_length, distance=distance)
>>> n_slots, map_zero = getHexObservations.contemplateTheDivisionsOfTime(
    probs, times, hoursAvailable=hours_available); 
# ==== figure out what to observe during that night
>>> best_slot  = getHexObservations.now( \
    n_slots, mapDirectory=data_dir, simNumber=trigger_id, mapZero=map_zero,
    skipJson=True)
>>> getHexObservations.makeObservingPlots( 
            n_slots, trigger_id, best_slot, data_dir)

# test out the new start_mjd and the rationalization of slot duration etc
#  30 minute per slot, 5 hexes/slot
# thursday mid night 57283 -0.5 call it 57282.6
>>> probs,times = getHexObservations.prepare(skymap, trigger_mjd, trigger_id, \
     data_dir, exposure_length=exposure_length, distance=distance, \
    start_mjd=57285.0, deltaTime =0.02083)
>>> n_slots, map_zero = getHexObservations.contemplateTheDivisionsOfTime(
    probs, times, hoursAvailable=hours_available, slotDuration=30.); 
>>> best_slot  = getHexObservations.now( \
    n_slots, mapDirectory=data_dir, simNumber=trigger_id, \
    mapZero=map_zero, maxHexesPerSlot=5, skipJson=True)
>>> getHexObservations.makeObservingPlots( 
            n_slots, trigger_id, best_slot, data_dir)
>>> reload(getHexObservations); reload(mapsAtTimeT)

>>> rahex,dechex = np.genfromtxt(hexFile, unpack=True)
>>> ra2, dec2 = np.genfromtxt("python/jack2/G184098-ra-dec-prob-mjd-slot.txt", unpack=True)
>>> ra_lmc, dec_lmc = np.genfromtxt("python/jack2/lmc.radec.txt", unpack=True)
>>> tra, tdec, tval = getHexObservations.tester(ra2, dec2, rahex, dechex, hexVals)
>> print  " in 3 slots", tval.sum()
>>> tra, tdec, tval = getHexObservations.tester(ra_lmc, dec_lmc, ra, dec, vals)
>> print  " in lmc ", tval.sum()
>>> ix = (ra2< 140); print ra2.size, ra2[ix].size
>>> tra, tdec, tval = getHexObservations.tester(ra2[ix], dec2[ix], ra, dec, vals)
>> print  " in 2.4 slots", tval.sum()

>>> ras, decs = np.genfromtxt("jack2/lasilla.txt", unpack=True, usecols=(1,2))
>>> 
>>> night 3 : 57285.0  data_dir="night3/"; start_mjd=57285.0 
>>> night 4 : 57286.0  data_dir="night4/"; start_mjd=57286.0 
>>> probs,times = getHexObservations.prepare(skymap, trigger_mjd, trigger_id,data_dir, exposure_length=exposure_length, distance=distance,start_mjd=57286.0, deltaTime =0.02083,onlyHexesAlreadyDone="../records/G184098-ra-dec-prob-mjd-slot.txt")
>>> n_slots, map_zero = getHexObservations.contemplateTheDivisionsOfTime(probs, times, hoursAvailable=hours_available, slotDuration=30.);
>>> best_slot  = getHexObservations.now(n_slots, mapDirectory=data_dir, simNumber=trigger_id,mapZero=map_zero, maxHexesPerSlot=5, skipJson=False)

More to come.
Jim Annis
Sept 11, 2015
