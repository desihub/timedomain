import numpy as np
import healpy as hp
import os
import decam2hp
import hexalate
import kasen_modelspace

import sourceProb
import modelRead

#
# mapsAtTimeT
#
# This file contains routines to compute probability maps
#   for a given burst time and given observation time 
#   in days since the burst.
#
# probabilityMaps
#   Use the mags and sourceProb objects to calculate
#   the limiting magnitude and probability maps at a given time
#   in a given filter and exposure time
#   The models needed as input are found via:
#       models = modelRead.getModels()
#   and then models are ready to roll
#
# totalProbability
#   Return the total probability of in the map at a given time
#   in a given filter and exposure. 
#       This is a thin wrapper about probabilityMaps
#
# manyDaysOfTotalProbability 
#   Return a list of times and total probabilites
#   at deltaTime starting at startOfDays and ending at endOfDays
#   using totalProbabilityAtTimeT
#       The deltaTime parameter sets the duration of a "slot",
#       in the langauge of getHexObservations
#
# oneDayOfTotalProbability 
#   A thin wrapper about manyDaysOfTotalProbability to do one 24 hour period
#
# probabilityMapSaver 
#   Given a list of times and total probabilites,
#   use probabilityMapsAtTimeT (if the probability is high enough)
#   to re-make the maps and save them to disk
#

# Over 11 days
#   calculate sourceProb.map and given a model, sm.calculateProb()
#   sum the probability in the map, append that number to an array,
#   return
#
#   On the subject of deltaTime:
#    deltaTime = 1./24.  ;# once per hour
#    deltaTime = 0.0223  ;# 32 minutes, = 1 hex (i,z,z,i) 180s+30s = 8 minutes 
#       so 4 hexes per 32 minute slot. 
#       More convenient for survey strategy planning (as 8 doesnt go into 60)
#   deltaTime = 0.223*2 ;# twice as long slots, (8 hexes/slot)


def oneDayOfTotalProbability (obs, deltaTime, start_mjd, 
        probTimeFile, gw_map_trigger, gw_map_strategy, gw_map_control) :

    burst_mjd      = gw_map_trigger.burst_mjd
    trigger_type   = gw_map_trigger.trigger_type
    spatial        = gw_map_trigger.ligo
    distance       = gw_map_trigger.ligo_dist
    distance_sig   = gw_map_trigger.ligo_dist_sig
    filter         = gw_map_strategy.working_filter
    exposure       = gw_map_strategy.summed_exposure_time
    kasen_fraction = gw_map_strategy.kasen_fraction 
    data_dir       = gw_map_control.datadir
    simple_distance = gw_map_trigger.distance
    simple_dist_err = gw_map_trigger.diststd


    # the work.
    start_of_days=0
    end_of_days=1
    #print "JTA debugging  ====="
    #end_of_days=0.1
    totalProbs,times,isDark = manyDaysOfTotalProbability(
        obs, burst_mjd, start_mjd, spatial, distance, distance_sig, 
        startOfDays=start_of_days, endOfDays=end_of_days, deltaTime=deltaTime, 
        probTimeFile=probTimeFile, trigger_type=trigger_type, 
        filter=filter, exposure=exposure, 
        kasen_fraction=kasen_fraction, data_dir = data_dir,
        simple_distance=simple_distance, simple_dist_err=simple_dist_err)

    return totalProbs,times, isDark

def manyDaysOfTotalProbability (
        obs, burst_mjd, start_mjd, spatial, distance, distance_sig, 
        startOfDays=0, endOfDays=11, deltaTime=0.0223, 
        probTimeFile="probTime.txt", trigger_type="bright",
        filter="i", exposure=90, 
        kasen_fraction=50, data_dir = "./", 
        simple_distance=100, simple_dist_err=30,
        verbose=True) :
    times = []
    totalProbs = []

    # in the language of getHexObservations:
    #   each slot is 32 minutes, each slot can hold 4 hexes
    dark = False
    isDark = []
    for time_since_start in np.arange(startOfDays,endOfDays,deltaTime) :
        print "================================== ",
        print "hours since Time Zero: {:5.1f}".format(time_since_start*24.),
        totalProb, sunIsUp = totalProbability(obs, burst_mjd, start_mjd, time_since_start, \
            spatial, distance, distance_sig, \
            filter=filter, exposure=exposure, \
            trigger_type=trigger_type, \
            kasen_fraction=kasen_fraction, data_dir = data_dir,
            simple_distance=simple_distance, simple_dist_err=simple_dist_err)
        if sunIsUp: 
            print "\t ... the sun is up"
        else:
            print ""
        times.append(time_since_start)
        totalProbs.append(totalProb)
        if not dark and not sunIsUp:
            dark = True ;# model sunset
        if dark and sunIsUp:
            dark = False ;# model sunrise
        isDark.append(dark)
    print "================================== "
    print "================================== ",
    print "full 24 hours calculated"
    print "================================== "
    totalProbs =np.array(totalProbs)
    times = np.array(times)
    isDark = np.array(isDark).astype(bool)
    # JTA successful half nights
    # an attempt to deal with half nights
    # if successful, these two booleans should go in yaml
    # somehow and be passed down to here
    darkCount = np.nonzero(isDark)[0]

    # informational
    ix, = np.where(totalProbs > 0)
    if verbose: print "total all-sky summed probability of detection (where > 0):"
    if verbose: print totalProbs[ix]
    if verbose: print "days since sunset:"
    if verbose: print times[ix]
    #if verbose: print "isDark:"
    #if verbose: print isDark.astype(int)
    #if verbose: print "dark slots count:",darkCount.size
    
#    print "===== times with total prob > 10**-2"
#    ix = totalProbs > 10**-2; 
#    if np.nonzero(ix)[0].size == 0 :
#        totalProbs = np.array([0,])
#        times = np.array([times[0],])
#    else :
#        totalProbs = totalProbs[ix]
#        times = times[ix]
#    print "total all-sky summed probability of detection (list1) and daysSinceBurst (list2)"
#    print totalProbs,"\n",times

    data = np.array([totalProbs, times, isDark.astype(int)]).T
    np.savetxt(probTimeFile, data, "%f %f %d")
    return totalProbs,times, isDark

#==============================================================
#
# core computation
#
def totalProbability(obs, burst_mjd, start_mjd, daysSinceStart, \
        spatial, distance, distance_sig, 
        filter="i", exposure=180, trigger_type="bright", 
        kasen_fraction=50, data_dir="./",
        simple_distance=100, simple_dist_err=30) :

    obs,sm,sunIsUp = probabilityMaps(obs, burst_mjd, start_mjd, daysSinceStart, \
        spatial, distance, distance_sig,
        filter, exposure, trigger_type=trigger_type, 
        kasen_fraction=kasen_fraction, data_dir = data_dir,
        simple_distance=simple_distance, simple_dist_err=simple_dist_err)
    if sunIsUp:
        totalProb = 0.0
    else :
        totalProb = (obs.map * sm.probMap).sum()
    return totalProb, sunIsUp

# drive the probability map calculations. In the end, distance only is used here
def probabilityMaps(obs, burst_mjd, start_mjd, daysSinceStart, \
        spatial, distance, distance_sig, 
        filter="i", exposure=180, trigger_type="bright", 
        kasen_fraction=50, data_dir="./", 
        simple_distance=100, simple_dist_err=30, verbose=True) :

    obs.resetTime(start_mjd+daysSinceStart)

    sunIsUp = obs.sunBrightnessModel(obs.sunZD)
    if sunIsUp: return obs, "sm", sunIsUp

    obs.limitMag(filter, exposure=exposure)
    if trigger_type == "bright" :
        # as of O3, let's not use source probability
        #   # we may need to rescale the light curve in the models
        #   models_at_t = modelRead.modelsAtTimeT (models, daysSinceBurst)
        #   model_scale, model_scale_time = models["scale"]
        #   AbsMag = sm.modelAbsoluteMagnitude
        #   new_scale = AbsMag - model_scale
        #   abs_mag = models_at_t[0] + new_scale
        #   sm.absMagMean = abs_mag
        # as of O3, let's not use source probability
        # turns out that in the case of distance_sig.sum()=0,
        # then the sm.calculateProb assumes a source ap mag of 20,
        # and builds a source probality map where 1 if limit_mag > 20.
        pass

    daysSinceBurst = start_mjd + daysSinceStart - burst_mjd
    #print "daysSinceBurst:",daysSinceBurst, burst_mjd
    #print "daysSinceBurst uses daysSinceStart: ", daysSinceStart,  start_mjd

    apparent_mag = kasen_modelspace.run_ap_mag_for_kasen_models (filter,
         simple_distance, simple_dist_err, daysSinceBurst,
         kasen_fraction, data_dir, doPlots=False, fast=True)

    sm=sourceProb.map(obs, type=trigger_type, apparent_mag_source=apparent_mag)
    result = sm.calculateProb(spatial, distance, distance_sig, verbose=verbose)
    if not result:
        sunIsUp = 1
    return obs,sm, sunIsUp

#==============================================================
#
# very  expensive hexelation of maps
#   thus it doublesa a routine that saves 
#   maps and hexelated maps to disk.
#
# for each time in the times,
# calculate the probabilities, 
#   and hexalate
#   then save the maps
#       times,probabilities are the output of manyDaysOfTotalProbability
#
# a onlyHexesAlreadyDone can point at a file that contains
# ra,decs that are already done, and these will replace the
# all sky hexes
# 
def probabilityMapSaver (obs, times, probabilities,
        gw_map_trigger, gw_map_strategy, gw_map_control,
        performHexalatationCalculation=True) :
    burst_mjd                  = gw_map_trigger.burst_mjd
    start_mjd                  = gw_map_trigger.start_mjd
    trigger_id                 = gw_map_trigger.trigger_id
    trigger_type               = gw_map_trigger.trigger_type
    ligo                       = gw_map_trigger.ligo
    distance                   = gw_map_trigger.ligo_dist
    distance_sig               = gw_map_trigger.ligo_dist_sig
    debug                      = gw_map_control.debug
    reject_hexes               = gw_map_control.reject_hexes
    data_dir                   = gw_map_control.datadir
    camera                     = gw_map_strategy.camera
    filter                     = gw_map_strategy.working_filter
    exposure                   = gw_map_strategy.summed_exposure_time
    kasen_fraction             = gw_map_strategy.kasen_fraction 
    simple_distance            = gw_map_trigger.distance
    simple_dist_err            = gw_map_trigger.diststd
    #onlyHexesAlreadyDone  = gw_map_control.this_tiling
    onlyHexesAlreadyDone  = []

    # one reads the tiling 9 hex centers as that is our default position
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    hexFile = gw_data_dir + "all-sky-hexCenters-"+camera+".txt"
    keep_flag = performHexalatationCalculation
    
    prob_slots = np.percentile(probabilities, 95)
    print("95th percentile of cumulative prob covered {}\n".format(prob_slots))
    ix, = np.where(probabilities > 0)
    # since we are now doing every slot in a night, counter 0 is at sunset
    counter = -1
    # let's keep track of when the sun is down
    made_maps_list  = np.array([])
    for time,prob  in zip(times, probabilities) :
        counter += 1
        performHexalatationCalculation = keep_flag
        if prob <= 0 : 
            performHexalatationCalculation = False
# terrible hack
        if debug:
            if prob <= prob_slots : 
                performHexalatationCalculation = False

        daysSinceStart = time
        #print "daysSinceStart:", daysSinceStart, time, start_mjd, burst_mjd
        obs,sm, isDark = \
            probabilityMaps( obs, burst_mjd, start_mjd, daysSinceStart, 
            ligo, distance, distance_sig,
            filter, exposure, trigger_type=trigger_type, 
            kasen_fraction=kasen_fraction, data_dir = data_dir,
            simple_distance=simple_distance, simple_dist_err=simple_dist_err)
        if obs.sunBrightnessModel (obs.sunZD ): 
            #print "\t\tThe sun is up. Continue", time, prob
            #print ""
            continue

        # Legend:
        # obs.ra, obs.dec, obs.map  = ligo map
        # obs.hx, obs.hy   = mcbryde projection of houar angle, dec
        # obs.maglim = limiting mag
        # sm.prob = limiting mag convolve abmag convolve volume
        # sm.probMap = total prob map
        # hexRa,hexDec,hexVals
        nameStem = os.path.join(data_dir, str(trigger_id) + "-{}".format(str(counter)))
        print "\t Writing map files as {} for time {:.3f} and prob {:.2e}".format(nameStem,time,prob)
        made_maps_list = np.append(made_maps_list, counter)

        name = nameStem + "-ra.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.ra)
        name = nameStem + "-dec.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.dec)
        name = nameStem + "-ha.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.ha)
        name = nameStem + "-hx.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.hx)
        name = nameStem + "-hy.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.hy)
        name = nameStem + "-x.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.x)
        name = nameStem + "-y.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.y)
        name = nameStem + "-map.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.map)
        name = nameStem + "-maglim.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.maglim)
        name = nameStem + "-maglim-global.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.maglimall)
# do we need this?
        #name = nameStem + "-prob.hp"
        #if os.path.exists(name): os.remove(name)
        #hp.write_map(name, sm.prob)
        name = nameStem + "-probMap.hp"
        if os.path.exists(name): os.remove(name)
        #sm.probMap = sm.probMap*gal
        try:
            hp.write_map(name, sm.probMap)
        except:
            # there is no prob map in the sm, so prob is zero everywhere
            hp.write_map(name, obs.map*0.0)

        treedata = decam2hp.buildtree(obs.ra*360./2/np.pi,obs.dec*360./2/np.pi,\
            nsides=hp.get_nside(obs.ra), recompute=True)
        tree = treedata[2]
        if performHexalatationCalculation :
            raHexen, decHexen, idHexen = hexalate.getHexCenters(hexFile)
            if len(onlyHexesAlreadyDone) > 0 :
                do_these = np.in1d(idHexen, onlyHexesAlreadyDone)
                do_these = np.nonzero(do_these)
                raHexen, decHexen, idHexen = \
                    raHexen[do_these], decHexen[do_these], idHexen[do_these]
            if len(reject_hexes) > 0 :
                dont_do_these = np.in1d(idHexen, reject_hexes)
                dont_do_these = np.nonzero(np.invert(dont_do_these))
                raHexen, decHexen, idHexen = \
                    raHexen[do_these], decHexen[do_these], idHexen[do_these]

        #    import pickle
        #    fd=open("saveme.pickle","wb")
        #    pickle.dump([obs, sm, raHexen, decHexen, idHexen, tree, camera], fd)
        #    fd.close()
            raHexen, decHexen, idHexen, hexVals, rank = \
                hexalate.cutAndHexalateOnRaDec ( obs, sm, raHexen, decHexen, idHexen, tree, camera, cutProbs=False)

            # where rank is to be understood as the indicies of the
            # ranked hexes in order; i.e., they have nothing to do with
            # raHexen, decHexen, hexVals except as a sorting key
            name = nameStem + "-hexVals.txt"
            if os.path.exists(name): os.remove(name)
            
            f = open(name,'w')
            for j in range(0,raHexen.size) :
                f.write("{:.6f}, {:.5f}, {:s}, {:.4e}, {:d}, {:.4f}\n".format(
                    raHexen[j],decHexen[j],idHexen[j],hexVals[j],rank[j],
                    (np.asfarray(rank*0.)+(start_mjd+time))[j]))
            f.close()
    
            #######################################################################################
            ##### Brout: new, we have to run again where we dont double up on probability in the ##
            ##### hexes, but this cannot be used for prioritization. ##############################
            raHexencut, decHexencut, idHexencut, hexValscut, rankcut = \
                hexalate.cutAndHexalateOnRaDec( obs, sm, raHexen, decHexen, idHexen, tree, camera, cutProbs=True)
            name = nameStem + "-hexVals-cutOverlappingProb.txt"
            if os.path.exists(name): os.remove(name)
            f = open(name,'w')
            for j in range(0,raHexencut.size) :
                f.write("{:.6f}, {:.5f}, {:s}, {:.4e}, {:d}, {:.4f}\n".format(
                    raHexencut[j],decHexencut[j],idHexencut[j],hexValscut[j],rankcut[j],
                    (np.asfarray(rankcut*0.)+(start_mjd+time))[j]))
            f.close()
    return made_maps_list

# pass apparent mag from yaml file to sourceProb


# Get the saved maps for each day and hour.
def readMaps (data_dir, simNumber, slot) :
    name = os.path.join(data_dir, str(simNumber) + "-{}".format(str(slot)))

    ra=hp.read_map(name+"-ra.hp");
    dec=hp.read_map(name+"-dec.hp");
    ha=hp.read_map(name+"-ha.hp");
    map=hp.read_map(name+"-map.hp");
    maglim=hp.read_map(name+"-maglim.hp");
    prob=hp.read_map(name+"-prob.hp");
    probMap=hp.read_map(name+"-probMap.hp");
    hx=hp.read_map(name+"-hx.hp");
    hy=hp.read_map(name+"-hy.hp");
    x=hp.read_map(name+"-x.hp");
    y=hp.read_map(name+"-y.hp");
    return ra, dec, map, maglim, prob, probMap, x,y, hx,hy

