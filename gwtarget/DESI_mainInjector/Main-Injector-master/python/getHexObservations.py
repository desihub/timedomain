import matplotlib
matplotlib.use("Agg"); # matplotlib.use("TkAgg") important for off line image generation
import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import hp2np
import decam2hp
import jsonMaker
import obsSlots
import gw_map_configure
import kasen_modelspace
import pickle
import glob
import shutil
import warnings
#
#       the routine mapsAtTimeT.oneDayOfTotalProbability
#           breaks the night into slots of 32 minutes duration
#       the routine mapsAtTimeT.probabilityMapSaver
#           makes maps for each of these slots
#       the routine getHexObservations.observing
#           once told how many slots can be used
#           returns "hoursObserving" containing ra,dec,prob,time
#               of the hexes to observe each slot
#       the routine obsSlots.observingStats
#           places the ra,dec,id,prob,times per slot onto single  lists
#           tells sum(prob) on the list
#       the routine getHexObservations.observingPlot
#           needs to be told simNumber, "slot", data_dir, nslots
#           and then will plot the sky limitMag+ligo contour
#           and all hexes (in the night) on the map for that slot
#

#========================================================================
#
# main routines: these seven are called in Dillon's recycler.py
#   make_maps
#   make_divisions_of_time
#   make_jsons
#   economics
#   makeObservingPlots      
#   nothingToObserveShowSomething
#   readMaps   
#
#========================================================================

# ==== prep the observation by calculating necessary maps
#       distance is a problem- the question is whether we want a horizon
#       distance or the estimated distance from LIGO?
#       >>> distance = 75. ;# Mpc, as an estimated horizon distance
#
#    deltaTime = 0.0223  => 32 minute slots  for izzi 90sec exp, 4 hex/slot
#    deltaTime = 0.0446  => 62 minute slots  for izzi 90sec exp, 8 hex/slot
#    deltaTime = 0.0417  => 60 minute slots  for izz 90sec exp, 10 hex/slot
#    deltaTime = 0.0208  => 30 minute slots  for izz 90sec exp,  5 hex/slot
#       The first is flexible to observing conditions,
#       while the second is x2 faster as there are less maps to hexalate
#       The third shows a different tiling/filter complement
#
#   exposure_length is only used in the computation of the limiting mag map
#
def make_maps(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) :
    import mapsAtTimeT
    import mags
    import modelRead
    import healpy as hp
    import os

    print "=================================================="
    print "                 make_maps "
    print "=================================================="

    skymap               = gw_map_trigger.skymap
    trigger_id           = gw_map_trigger.trigger_id
    trigger_type         = gw_map_trigger.trigger_type
    distance             = gw_map_trigger.distance
    dist_err             = gw_map_trigger.diststd
    days_since_burst     = gw_map_trigger.days_since_burst
    burst_mjd            = gw_map_trigger.burst_mjd
    start_mjd            = gw_map_trigger.start_mjd

    camera               = gw_map_strategy.camera
    exposure_list        = gw_map_strategy.exposure_list
    filter_list          = gw_map_strategy.filter_list
    tiling_list          = gw_map_strategy.tiling_list
    maxHexesPerSlot      = gw_map_strategy.maxHexesPerSlot
    overhead             = gw_map_strategy.overhead 
    kasen_fraction       = gw_map_strategy.kasen_fraction
    use_teff             = gw_map_strategy.use_teff
    summed_exposure_time = gw_map_strategy.summed_exposure_time 
    working_filter       = gw_map_strategy.working_filter  

    resolution           = gw_map_control.resolution
    this_tiling          = gw_map_control.this_tiling
    reject_hexes         = gw_map_control.reject_hexes
    debug                = gw_map_control.debug
    data_dir             = gw_map_control.datadir
    snarf_mi_maps        = gw_map_control.snarf_mi_maps
    mi_map_dir           = gw_map_control.mi_map_dir

    print "burst_mjd = {:.2f} at a distance of {:.1f} Mpc".format(
        burst_mjd, distance)

    print "day of burst = {:.2f}" .format(start_mjd),
    print ", calculation starting at  sunset --- with an additional {:.1f} days added\n".format(
         days_since_burst)
    start_mjd = start_mjd + days_since_burst

    if snarf_mi_maps:
        if mi_map_dir == "/data/des41.a/data/desgw/O3FULL/Main-Injector/OUTPUT/O3REAL/":
            mi_map_dir = mi_map_dir + trigger_id 
        # get everything we need from a different directory
        print "copying from {}/ to {}/".format(mi_map_dir, data_dir)
        os.system("cp {}/*hp {}/".format(mi_map_dir, data_dir))
        os.system("cp {}/*probabilityTimeCache*txt {}/".format(mi_map_dir, data_dir))
        os.system("cp {}/*-hexVals.txt {}/".format(mi_map_dir, data_dir))
        os.system("cp {}/*-hexVals-cutOverlappingProb.txt {}/".format(mi_map_dir, data_dir))
        return

    print "\t cleaning up old files",
    if os.path.exists(data_dir+"/json") :
        shutil.rmtree("json/")
    files = glob.glob(data_dir+"/*png"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*jpg"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*json"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*hp"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*txt"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*pickle"); 
    for f in files: os.remove(f)
    print "\t done cleaning up"

# we need to understand the night, how long it is, etc
    print "\t obs slots starting"
    answers = obsSlots.slotCalculations( start_mjd, exposure_list, tiling_list, 
        overhead, hexesPerSlot=maxHexesPerSlot, camera=camera) 
    hoursPerNight = answers["hoursPerNight"] ;# in minutes
    slotDuration = answers["slotDuration"] ;# in minutes
    deltaTime = slotDuration/(60.*24.) ;# in days
    print "\t obs slots done: hours in night: {:.1f}, slot duration: {:.1f} minutes, ".format(
        hoursPerNight, slotDuration),
    print("hexes per slot: {}\n".format(maxHexesPerSlot))

    #  ==== run Rob's evaluate the kasen sim universe code
    print "\t examining Kasen universe KN models for coverage"
    night_dur,sunset,sunrise = mags.findNightDuration(start_mjd, camera)
    midnight_since_burst = 24*(sunset+night_dur/2. - burst_mjd) 
    apparent_mag = kasen_modelspace.run_ap_mag_for_kasen_models (
        working_filter,
        distance, dist_err, 
        midnight_since_burst,
        kasen_fraction, data_dir,
        fast = False)
    print "\t\t at a time halfway through night, {:.2f} days after merger:".format(midnight_since_burst/24.)
    print "\t\t {}% requires observations at {} <= {:5.2f}\n".format(kasen_fraction, filter_list[0], apparent_mag)

        
    # ==== get the neutron star explosion models
    #models = modelRead.getModels()

    # === prep the maps
    ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)
    ligo_dist, ligo_dist_sig, ligo_dist_norm  = \
        distance*np.ones(ra.size), np.zeros(ra.size), np.zeros(ra.size)
    warnings.filterwarnings("error")
    try :
        junk,junk,ligo_dist =hp2np.hp2np(skymap, degrade=resolution, field=1)
        junk,junk,ligo_dist_sig =hp2np.hp2np(skymap, degrade=resolution, field=2)
        junk,junk,ligo_dist_norm =hp2np.hp2np(skymap, degrade=resolution, field=3)
    except RuntimeWarning:
        pass
    except:
        print "\t !!!!!!!! ------- no distance information in skymap ------ !!!!!!!!"

    # details of the observations
    obs = mags.observed(ra,dec,ligo, start_mjd, verbose=False)
    obs.limitMag(working_filter,exposure=summed_exposure_time)
    print "finished setting up exposure calculation"

    # ==== calculate maps during a full night of observing
    # essentially, give me the total prob, source + ligo for a 1 filter pass
    probabilityTimesCache = os.path.join(data_dir,\
        "probabilityTimesCache_"+str(trigger_id)+".txt")
    probs,times,isdark = mapsAtTimeT.oneDayOfTotalProbability(
        obs, deltaTime, start_mjd, probabilityTimesCache,
        gw_map_trigger, gw_map_strategy, gw_map_control)

    made_maps_list = mapsAtTimeT.probabilityMapSaver (obs, times, probs, 
        gw_map_trigger, gw_map_strategy, gw_map_control)

    gw_map_results.probability_per_slot = probs
    gw_map_results.time_of_slot         = times
    gw_map_results.slotDuration         = slotDuration
    gw_map_results.hoursPerNight        = hoursPerNight
    gw_map_results.isdark               = isdark
    gw_map_results.made_maps_list       = made_maps_list
    gw_map_results.moonRa               = obs.moonRa*360./2/np.pi
    gw_map_results.moonDec              = obs.moonDec*360./2/np.pi
    # where 0 = full, 90 equals half, and  180 = new
    gw_map_results.moonIllumination     = (180-obs.moonPhase)/180.
    pickle.dump(made_maps_list, open("made_maps.pickle","wb"))
    pickle.dump(gw_map_results, open("gw_map_results.pickle","wb"))

    return 

# ==== figure out what to observe
#
# Another fast routine
#   basic for this routine is how many hexes per slot
#
#   this_tiling and reject_hexes are inverse in thier effect:
#   only hexes within 36" of the ra,decs of the centers of this_tiling
#   are to be observed, while hexes wihtin within 36" of the ra,decs
#   of the centers of reject_hexes are skipped
#
#   if a doneFile is given, which should be
#   something with ra,dec in columns 1,2, like 
#       G184098-ra-dec-prob-mjd-slot.txt
#   then those ra,decs are interpreted as hex centers,
#   and hexes in the hexVals maps within 36" of the ra,decs
#   are removed- they are done, the question is what to observe next
#
def make_hexes( gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, 
    start_slot = -1, do_nslots=-1) :

    print "=================================================="
    print "                 make_hexes "
    print "=================================================="

    trigger_id = gw_map_trigger.trigger_id
    trigger_type  = gw_map_trigger.trigger_type
    maxHexesPerSlot  = gw_map_strategy.maxHexesPerSlot
    exposure_list  = gw_map_strategy.exposure_list
    filter_list  = gw_map_strategy.filter_list
    hoursAvailable  = gw_map_strategy.hoursAvailable
    data_dir = gw_map_control.datadir
    if start_slot != -1 or do_nslots != -1 :
        gw_map_control.start_slot = start_slot
        gw_map_control.do_nslots = do_nslots

    probs = gw_map_results.probability_per_slot 
    times = gw_map_results.time_of_slot 
    slotDuration = gw_map_results.slotDuration   
    hoursPerNight = gw_map_results.hoursPerNight 
    if (hoursPerNight == False)  :
        #reuse_results(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) 
        gw_map_results = reuse_results(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) 
        probs = gw_map_results.probability_per_slot 
        times = gw_map_results.time_of_slot 
        hoursPerNight = gw_map_results.hoursPerNight 
        slotDuration = gw_map_results.slotDuration   

# I believe as of Jan 2020, n_slots is related primarily to slotDuration and the length of the night
    n_slots, mapZero  = make_divisions_of_time (
        probs, times, hoursPerNight, hoursAvailable, slotDuration) 
    print "\t tonight is {:.2f} hours long with {} slots".format(hoursPerNight, n_slots)
        
    # if the number of slots is zero, nothing to observe or plot
    if n_slots == 0: 
        return 0
    # compute the observing schedule
    print "=============>>>>  observing"
    if start_slot != -1 or do_nslots != -1 :
        print "\t tonight we will use {} slots starting at {}".format(do_nslots, start_slot)
    #print "hoursObserving=obsSlots.observing("
    #print trigger_id,n_slots, data_dir, "mapZero=",mapZero, "maxHexesPerSlot =", maxHexesPerSlot,
    #print "do_nslots =", do_nslots, "start_slot=",start_slot
    #print ")"
    #raise Exception("here")
    hoursObserving=obsSlots.observing(
        trigger_id,n_slots, data_dir, mapZero=mapZero,
        maxHexesPerSlot = maxHexesPerSlot, do_nslots = do_nslots, start_slot=start_slot)

    # print stats to screen
    print "\n=============>>>>  observingStats from *ra-dec-id-* file"
    # save results to the record -- here is made the ra-dec-id- file
    writeObservingRecord(hoursObserving,   data_dir, gw_map_trigger, gw_map_control, gw_map_strategy)
    #ra,dec,id,prob,mjd,slotNumbers,islots = obsSlots.observingStats(hoursObserving, mapZero, do_nslots, start_slot)
    ra,dec,id,prob,mjd,slotNumbers,islots = obsSlots.observingStatsFromRaDecFile( trigger_id, data_dir, 
        hoursObserving, mapZero, do_nslots, start_slot)
    maxProb_slot = obsSlots.maxProbabilitySlot(prob,slotNumbers)
    hoursObserving["maxSlot"] = maxProb_slot
    pickle.dump(hoursObserving, open("{}-hoursObserving.pickle".format(trigger_id),"w"))
    # if the length of ra is one and value zero, nothing to observe or plot
    if ra.size == 1 and prob.sum() == 0: 
        print '\t ZERO PROBABILITY observed!'
        sum_ligo_prob = 0.
    else :
        # shall we measure the total ligo probability covered?
        # Lets find out how well we did in covering Ligo probability
        try:
            sum_ligo_prob = how_well_did_we_do( gw_map_trigger, gw_map_strategy, gw_map_control)
        except:
            print '\t ZERO PROBABILITY observed!'
            sum_ligo_prob = 0.
        
        if do_nslots == -1 :
            cumulPlot(trigger_id, data_dir) 
            cumulPlot2(trigger_id, data_dir) 

    gw_map_results.n_slots = n_slots
    gw_map_results.first_slot = mapZero
    gw_map_results.best_slot = maxProb_slot
    gw_map_results.slot_numbers = np.unique(slotNumbers)
    gw_map_results.sum_ligo_prob = sum_ligo_prob
    gw_map_results.n_hexes = ra.size
    if ra.size == 1 and prob.sum() == 0: 
        gw_map_results.n_hexes = 0
    pickle.dump(gw_map_results, open("gw_map_results.pickle","wb"))

    return 

# Make the json files
# you will have to add a
# import gw_helper to your code and use the function
# gw_helper.setup_tilings()

def make_jsons(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) :
    import shutil
    from gw_helper import gw_helper
    print "=================================================="
    print "                 make_jsons "
    print "=================================================="

    skymap        = gw_map_trigger.skymap
    trigger_id    = gw_map_trigger.trigger_id
    trigger_type  = gw_map_trigger.trigger_type
    exposure_list = gw_map_strategy.exposure_list
    filter_list   = gw_map_strategy.filter_list
    tiling_list   = gw_map_strategy.tiling_list
    propid        = gw_map_strategy.propid
    data_dir      = gw_map_control.datadir

    print "\t cleaning up old jsons"
    if os.path.exists(data_dir+"/json") :
        shutil.rmtree("json/")
    files = glob.glob(data_dir+"/*json"); 
    for f in files: os.remove(f)
    os.mkdir(data_dir+"/json")

    hoursPerNight = gw_map_results.hoursPerNight
    if (hoursPerNight == False)  :
        gw_map_results = reuse_results(
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) 
    n_hexes       = gw_map_results.n_hexes

    if n_hexes == 0 :
        print "\t no hexes to turn into json files"
        return 

    ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)

    list_of_jsons = turnObservingRecordIntoJSONs(
        ra,dec,id,prob,mjd,slotNum, trigger_id,
        exposure_list, filter_list, tiling_list, trigger_type, skymap, data_dir, propid) 

    for tiling in tiling_list[1:] :
        if tiling == 0: ra_shift, dec_shift =  0.00, 0.00
        if tiling == 1: ra_shift, dec_shift =  0.02, 0.00
        if tiling == 2: ra_shift, dec_shift =  0.00, 0.02
        if tiling == 3: ra_shift, dec_shift = -0.02, 0.00
        if tiling == 4: ra_shift, dec_shift =  0.00,-0.02
        print("\t tiling {} ra,dec shift by {},{}".format(tiling, ra_shift, dec_shift))
        gw_helper.setup_tilings(
            list_of_jsons,
            ra_shift = [ra_shift], dec_shift = [dec_shift],
            tilings = tiling,
            prefix='tiling_{}'.format(tiling),
            overwrite=False,
        )
    for file in list_of_jsons:
        dir = os.path.dirname(file)
        file = os.path.basename(file)
        shutil.copyfile(dir+"/"+file, "json/tiling_{}_".format(tiling_list[0]) + file)

#
# ====== there are possibilities. Show them.
#
def makeGifs (gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) :
    print "\n=================================================="
    print "                 make_gifs "
    print "=================================================="
    hoursPerNight = gw_map_results.hoursPerNight
    if (hoursPerNight == False)  :
        gw_map_results = reuse_results(
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) 
    n_hexes       = gw_map_results.n_hexes


    # make gif centered on ra=0,dec=0, all sky
    if  gw_map_control.allSky == True :
        n_plots = makeObservingPlots(
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, allSky=True)

    # make gif centered on hexes
    if  gw_map_control.centeredSky== True :
        print "n_hexes=",n_hexes
        if n_hexes != 0:
            n_plots = makeObservingPlots(
                gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, allSky=False)
        else :
            print "\t skipping centered gif as no hexes observed"

def makeObservingPlots( gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, allSky=True) :

    trigger_id      = gw_map_trigger.trigger_id
    camera          = gw_map_strategy.camera
    data_dir        = gw_map_control.datadir
    start_slot      = gw_map_control.start_slot
    do_nslots       = gw_map_control.do_nslots
    gif_resolution  = gw_map_control.gif_resolution

    # we are doing a quick run, avoiding most calculations as they are already done and on disk
    if (gw_map_results.hoursPerNight == False)  :
        #reuse_results(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, get_slots=True)
        gw_map_results = reuse_results(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) 
    n_slots         = gw_map_results.n_slots
    n_hexes         = gw_map_results.n_hexes
    first_slot      = gw_map_results.first_slot
    best_slot       = gw_map_results.best_slot
    probs           = gw_map_results.probability_per_slot 
    times           = gw_map_results.time_of_slot 
    moonRa          = gw_map_results.moonRa
    moonDec         = gw_map_results.moonDec
    moonIllumination = gw_map_results.moonIllumination

    made_maps_list = (gw_map_results.made_maps_list).astype(int)

    print "\n================ >>>>>>>>>>>>>>>>>>>>> =================== "
    print "makeObservingPlots(",n_slots, trigger_id, best_slot,data_dir," )"
    print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    print "We're going to do {} slots with best slot={}".format(made_maps_list.size, best_slot)
    figure = plt.figure(1,figsize=(8.5*1.618,8.5))

    # first, make the probability versus something plot
    ra,dec,id,prob,slotMjd,slotNumbers,dist = obsSlots.readObservingRecord( trigger_id, data_dir)

    # now make the hex observation plots
    counter = 1   ;# already made one
    #for i in np.unique(slotNumbers) :

    label = ""
    if allSky == False: label="centered-"

    for i in made_maps_list:
        
        observingPlot(figure,trigger_id,i,data_dir, n_slots, n_hexes, camera, 
            allSky=allSky, gif_resolution=gif_resolution, 
            moonRa=moonRa, moonDec=moonDec, moonIllumination=moonIllumination)
        name = str(trigger_id)+"-{}observingPlot-{}.png".format(label,i)
        plt.savefig(os.path.join(data_dir,name))
        counter += 1
        counter+= equalAreaPlot(figure,i,trigger_id,data_dir)

    string = "$(ls -v {}observingPlot*)  {}_{}animate.gif".format(data_dir+'/'+trigger_id+'-'+label, data_dir+'/'+trigger_id, label)
    print string
    os.system("convert  -delay 40 -loop 0  " + string)

    # return the number of plots made
    return counter

def how_well_did_we_do(gw_map_trigger, gw_map_strategy, gw_map_control) :
    skymap = gw_map_trigger.skymap
    trigger_id = gw_map_trigger.trigger_id
    camera = gw_map_strategy.camera
    resolution = gw_map_control.resolution
    data_dir = gw_map_control.datadir

    ra = gw_map_trigger.ligo_ra
    dec = gw_map_trigger.ligo_dec
    ligo = gw_map_trigger.ligo
    name = os.path.join(data_dir, str(trigger_id) + "-ra-dec-id-prob-mjd-slot-dist.txt")
    raH, decH = np.genfromtxt(name, unpack=True, usecols=(0,1))
    treedata = decam2hp.buildtree(ra, dec, resolution, recompute=True) 
    tree = treedata[2] 
    sum = decam2hp.hexalateMapWithoutOverlap(ra,dec,ligo,tree, raH,decH,camera, verbose=False) 
    print "\nTotal Ligo probability covered by hexes observed: {:.3f}%".format(sum.sum()*100.)
    print "   (from decam2hp.sumHexalatedMap of -ra-dec-id file)\n"
    return sum.sum()


# ========== do simple calculations on how to divide the night
#
#   one of the questions is how many hours to devote to observations
#       hoursAvailable,  another is the slot duration
#
def make_divisions_of_time (
        probs, times, hoursPerNight= 10., hoursAvailable=6, slotDuration=30.) :
    if hoursAvailable > hoursPerNight:
        hoursAvailable = hoursPerNight
        
    # if the number of slots is zero, nothing to observe or plot
    if np.size(times) == 0 : return 0,0
    if probs.sum() < 1e-9 : return 0,0
    verbose = 0
    n_slots = obsSlots.findNSlots(hoursAvailable,slotDuration=slotDuration)
    n_maps = times.size
    if verbose: print n_slots, n_maps
    if n_maps == n_slots : 
        mapZero = 0
    elif n_maps < n_slots : 
        mapZero = 0
        n_slots = n_maps
    elif n_maps > n_slots :
        mapZero = obsSlots.findStartMap ( probs, times, n_slots )
    else :
        raise Exception ("no possible way to get here")
    print "\n=============>>>>  make_divisions_of_time:"
    print "\t n_maps = {}, n_slots = {}, mapZero = {}, prob_max = {:.6}".format(
        n_maps, n_slots, mapZero, probs.max())
    return n_slots, mapZero

# make_hexes computes slot information, so get_slots = false
# make_observingPlots needs slot information, so get_slots = true
def reuse_results(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, get_slots=False) :
        gw_map_results = pickle.load(open("gw_map_results.pickle","rb"))
        return gw_map_results

# delete this when comfortable
def old_reuse_results(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, get_slots=False) :
        trigger_id      = gw_map_trigger.trigger_id
        data_dir        = gw_map_control.datadir
        camera           = gw_map_strategy.camera
        probabilityTimesCache = os.path.join(data_dir,\
        "probabilityTimesCache_"+str(trigger_id)+".txt")
        print "=============>>>> Reuse results via reuse_results"
        #print "\t using probabilities, times, and maps from", probabilityTimesCache
        if os.stat(probabilityTimesCache).st_size == 0 :
            probs, times = np.array([0,]),np.array([0,])
            print "\t got nothing- we're skipping this one"
        else :
            data = np.genfromtxt(probabilityTimesCache, unpack=True)
            probs, times,isdark = data[0],data[1],data[2]
        gw_map_results.probability_per_slot = probs
        gw_map_results.time_of_slot = times
        gw_map_results.isdark = isdark
        gw_map_results.made_maps_list = pickle.load(open("made_maps.pickle","rb"))

        exposure_list   = gw_map_strategy.exposure_list
        burst_mjd       = gw_map_trigger.burst_mjd 
        start_mjd       = gw_map_trigger.start_mjd 
        maxHexesPerSlot = gw_map_strategy.maxHexesPerSlot
        overhead        = gw_map_strategy.overhead
        answers = obsSlots.slotCalculations(start_mjd, exposure_list, overhead,
            maxHexesPerSlot, camera=camera)
        hoursPerNight = answers["hoursPerNight"] ;# in minutes
        slotDuration = answers["slotDuration"] ;# in minutesk
        gw_map_results.slotDuration = slotDuration
        gw_map_results.hoursPerNight = hoursPerNight

        if get_slots:
            hoursObserving = pickle.load(open("{}-hoursObserving.pickle".format(data_dir+'/'+trigger_id),"r"))
            ra,dec,id,prob,mjd,slotNum,islot = obsSlots.slotsObservingToNpArrays(hoursObserving)
            gw_map_results.n_slots    = hoursObserving["nslots"]
            gw_map_results.first_slot = hoursObserving["mapZero"]
            gw_map_results.best_slot  = hoursObserving["maxSlot"]
            gw_map_results.n_hexes    = ra.size

# ===== The economics analysis
#
#   area_left is th enumber of hexes we have left to observe this season
#   days_left is the number of days left in the season
#   rate is the effective rate of triggers
#       p_gw is that for which the table cumul_table_pgw50.txt was  made.
#
def economics (simNumber, best_slot, mapDirectory, 
        area_left=200., days_left=60., rate=1/30.,  p_gw = 0.10) :
    import healpy as hp
    import cumul
    import des_optimization
    import os
    gw_data_dir = os.environ["DESGW_DATA_DIR"]
    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(mapDirectory, simNumber, best_slot)

    area_bar_p,area_bar = np.genfromtxt(
        gw_data_dir+"/area_bar_table.txt",unpack=True)
    avge_cumu_area,avge_cumu = np.genfromtxt(
        gw_data_dir+"/cumul_table_pgw10.txt",unpack=True)

    obsProb = ligo*prob
    nsides = hp.get_nside(obsProb)
    # max area viewable by Blanco at one time is 11734. sq-degrees
    max_area=11734.
    area, cum_prob  = cumul.area(ra,dec,obsProb, p_gw, nsides, max_area=max_area)
    area_to_cover_p_gw = area
    #print avge_cumu_area
    #print area
    ix = np.searchsorted(avge_cumu_area, area)
    if ix >= avge_cumu_area.size :
        fraction_of_sims_better_than_this_trigger = 1.0
    else :
        fraction_of_sims_better_than_this_trigger = avge_cumu[ix]

    prob, N_max = des_optimization.evaluate_average_event(
        area_left, days_left, rate, avge_cumu, avge_cumu_area, area_bar, area_bar_p)

    if fraction_of_sims_better_than_this_trigger < 1./N_max :
        area, cum_prob = cumul.area(ra,dec,obsProb, prob, nsides, max_area=max_area)
        if area>area_left:
            print "\t maxing out area: \t {:.3f} -> ".format( cum_prob),
            cum_prob = cumul.probability_covered(ra,dec,obsProb, area_left, nsides, max_area=max_area)
            print "{:.3f}".format(cum_prob)
            area=area_left
    else :
        print "\t ignore event"
        area = 0
        prob = 0

    probability_covered = cum_prob
    quality = fraction_of_sims_better_than_this_trigger 
    return probability_covered, area, area_to_cover_p_gw, quality

#
# no, no, no, we actually can see something: lets see the best plots
#
#   raMap, decMap, ligoMap, maglimMap, probMap, haMap, xMap,yMap, hxMap,hyMap = readMaps(
#   ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = readMaps(
def readMaps(mapDir, simNumber, slot) :
    import healpy as hp
    # get the maps for a reasonable slot
    name = os.path.join(mapDir, str(simNumber) + "-"+str(slot))
    #print "\t reading ",name+"-ra.hp  & etc"
    raMap     =hp.read_map(name+"-ra.hp", verbose=False);
    decMap    =hp.read_map(name+"-dec.hp", verbose=False);
    haMap     =hp.read_map(name+"-ha.hp", verbose=False);
    xMap      =hp.read_map(name+"-x.hp", verbose=False);
    yMap      =hp.read_map(name+"-y.hp", verbose=False);
    hxMap     =hp.read_map(name+"-hx.hp", verbose=False);
    hyMap     =hp.read_map(name+"-hy.hp", verbose=False);
    ligoMap   =hp.read_map(name+"-map.hp", verbose=False);
    maglimMap =hp.read_map(name+"-maglim.hp", verbose=False);
    probMap   =hp.read_map(name+"-probMap.hp", verbose=False);
    haMap=haMap/(2*np.pi/360.)
    raMap=raMap/(2*np.pi/360.)
    decMap=decMap/(2*np.pi/360.)
    return raMap, decMap, ligoMap, maglimMap, probMap, \
        haMap, xMap, yMap, hxMap, hyMap

#========================================================================
# 
# support routines
# 
#========================================================================
# for economics analysis
def time_cost_per_hex (nvisits, overhead, exposure_length) :
    tot_exptime = (np.array(overhead)+np.array(exposure_length)).sum
    time_cost_per_hex = nvisits * tot_exptime #sec
    return time_cost_per_hex
    
# for economics analysis
def area_left (area_per_hex, time_budget, time_cost_per_hex) :
    area_per_hex * (time_budget * 3600)/(time_cost_per_hex)
    return area_per_hex
#time_cost_per_hex = nvisits * nexposures * (overhead + exposure_length) #sec
#area_left =  area_per_hex * (time_budget * 3600)/(time_cost_per_hex)

# place holder for the code brought from desisurvey...
def hoursPerNight (mjd, camera="decam") :
    import mags
    night,sunset,sunrise = mags.findNightDuration(mjd, camera)
    night = night*24.
    return night


#===================================================================
#
    
def turnObservingRecordIntoJSONs(
        ra,dec,id,prob,mjd,slotNumbers, trigger_id,
        exposure_list, filter_list, tiling_list, trigger_type, skymap, mapDirectory, propid) :
    seqtot =  ra.size
    seqzero,seqnum = 0,0

    # write slot json files
    list_of_jsons = []
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        slotMJD = mjd[ix][0]  ;# just get first mjd in this slot
        tmpname, name = jsonUTCName(slot, slotMJD, trigger_id, mapDirectory)
        tmpname = mapDirectory + tmpname
        name = mapDirectory + name
        jsonMaker.writeJson(ra[ix],dec[ix],id[ix],
            seqzero, seqnum, seqtot, exposureList= exposure_list, 
            filterList= filter_list, tilingList = tiling_list, trigger_id = trigger_id, 
            trigger_type=trigger_type, propid=propid, skymap=skymap, jsonFilename=tmpname)
        desJson(tmpname, name, mapDirectory, slotMJD) 
        seqzero += ra[ix].size
        list_of_jsons.append(mapDirectory+name)
    return list_of_jsons
        
# verbose can be 0, 1=info, 2=debug
def desJson(tmpname, name, data_dir, start_time, verbose = 1) :
    import os
    import json
    import gwwide
    warnings.filterwarnings("ignore")

    fd = open(tmpname,"r")
    gw_queue = json.load(fd); fd.close()
    # seconds to leave before hitting the Blanco limits
    time_buffer= 300.
    fixed_queue = gwwide.gwwide([],gw_queue, start_time, time_buffer, sort=True)
    fd = open(name,"w")
    print "make ericJson {}".format(name)
    json.dump(fixed_queue, fd, indent=4); fd.close()
    os.remove(tmpname)

def jsonUTCName (slot, mjd, simNumber, mapDirectory) :
    time = utcFromMjd(mjd, alt=True)
    tmpname, name = jsonName(slot, time, simNumber, mapDirectory)
    return tmpname, name

def utcFromMjd (mjd, alt=True) :
    year,month,day,hour,minute = utc_time_from_mjd(mjd)
    time = "UT-{}-{:02d}-{:02d} {:02d}:{:02d}:00".format(year,month,day,hour,minute)
    if alt:
        time = "UT-{}-{:02d}-{:02d}_{:02d}:{:02d}:00".format(year,month,day,hour,minute)
    return time

def utc_time_from_mjd (mjd) :
    from pyslalib import slalib
    date = slalib.sla_djcl(mjd)
    year = np.int(date[0])
    month= np.int(date[1])
    day = np.int(date[2])
    hour = np.int(date[3]*24.)
    minute = np.int( (date[3]*24.-hour)*60.  )
    return year, month, day, hour, minute

def jsonName (slot, utcString, simNumber, mapDirectory) :
    slot = "-{}-".format(np.int(slot))
    #tmpname = os.path.join(mapDirectory, str(simNumber) + slot + utcString + "-tmp.json")
    #name = os.path.join(mapDirectory, str(simNumber) + slot + utcString + ".json")
    tmpname = str(simNumber) + slot + utcString + "-tmp.json"
    name = str(simNumber) + slot + utcString + ".json"
    return tmpname, name

def jsonFromRaDecFile(radecfile, nslots, slotZero, 
        hexesPerSlot, simNumber, mjdList, trigger_id, trigger_type, data_dir) :
    ra,dec = np.genfromtxt(radecfile, unpack=True,usecols=(0,1),comments="#",propid="propid")

    seqtot =  ra.size
    seqzero = 0

    # instead, just reorder the ra,dec before feeding to this routine
    #ix = np.argsort(ra)

    counter = 0
    slot = slotZero
    slotRa = np.array([])
    slotDec = np.array([])
    for i in range(0,ra.size) :
        slotRa = np.append(slotRa, ra[i])
        slotDec = np.append(slotDec, dec[i])
        counter += 1
        if counter == hexesPerSlot :
            tmpname, name = jsonName(slot, mjdList[slot-slotZero], 
                simNumber,data_dir)
            jsonMaker.writeJson(slotRa,slotDec, 
                simNumber, seqzero+(hexesPerSlot*(slot-slotZero)), 
                seqtot, trigger_id, trigger_type, jsonFilename=tmpname, propid=propid)
            desJson(tmpname, name, data_dir) 
            counter = 0
            slot += 1
            slotRa = np.array([])
            slotDec = np.array([])
    if counter > 0 :
        tmpname, name = jsonName(slot, mjdList[slot-slotZero], 
            simNumber,data_dir)
        jsonMaker.writeJson(slotRa,slotDec, 
            simNumber, seqzero+(hexesPerSlot*(slot-slotZero)), 
            seqtot, trigger_id, trigger_type, jsonFilename=tmpname, propid=propid)
        desJson(tmpname, name, data_dir) 
        
# ==================================
# plotting 
# ==================================

def cumulPlot(trigger_id, data_dir) :
    from scipy.interpolate import interp1d
    verbose = False

    ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)
    u_slotNum = np.unique(slotNum)
    u_dur = np.unique(mjd)
    duration = 0.0000001
    if len(u_dur) > 1:
        duration = (u_dur[1]-u_dur[0])*24*60.

    new_x, new_y, new_mjd, new_hour = np.array([]), np.array([]), np.array([]), np.array([])
    for u in u_slotNum :
        ix = slotNum == u
        y = prob[ix].sum()
        new_x = np.append(new_x, u)
        new_y = np.append(new_y, y)
        new_mjd = np.append(new_mjd, mjd[ix].min())
        year,month,day,hour,minute = utc_time_from_mjd(mjd[ix].min())
        new_hour = np.append(new_hour, hour+minute/60.)

    plt.close()
    fig, ax1 = plt.subplots() ; 
    ax2 = ax1.twinx()
    ax3 = ax1.twiny()

    ax1.scatter(new_x,new_y,c="k", s=15, zorder=10)
    ax1.plot(new_x,new_y,c="b", zorder=9)
    ax1.set_xlabel("slot number")
    ax1.set_ylabel("LIGO probability in slot")
    ax1.set_ylim(0,new_y.max()*1.1)
    ax1.grid(alpha=0.2)

    plotProb, plotSlot,plotN = np.array([]), np.array([]), np.array([])
    for i in np.unique(slotNum) :
        ix = slotNum == i
        if prob[ix].sum() > 0 :
            plotN = np.append(plotN, prob[ix].size)
            plotSlot = np.append(plotSlot,i)
            plotProb = np.append(plotProb,100.*prob[ix].sum())
    ax1.text(0.99,0.02,"total probability: {:5.1f}%".format(prob.sum()*100.),
        transform = ax1.transAxes,   horizontalalignment='right',
        verticalalignment='center')
    avghex = str( np.round(plotN.mean(),1) )
    ax1.text(0.99,0.07,"slot duration: {:5.1f}$^m$".format(duration),
        transform = ax1.transAxes,   horizontalalignment='right',
        verticalalignment='center')
    ax1.text(0.99,0.12,"hexes per slot: {}".format(avghex),
        transform = ax1.transAxes,   horizontalalignment='right',
        verticalalignment='center')

    ix = slotNum == u_slotNum[0]
    year,month,day,hour,minute = utc_time_from_mjd(mjd[ix].min())
    new_cy = np.cumsum(new_y)
    ax2.scatter(new_x,new_cy,color="k",s=15,zorder=10)
    ax2.plot(new_x,new_cy,color="g",zorder=9)
    ax2.tick_params(axis='y', labelcolor="g")
    ax2.set_ylabel("cumulative LIGO probability",color="g")
    ax2.set_ylim(0,new_cy.max())
    fig.tight_layout()
    qtiles = new_cy/new_cy[-1]
    if verbose: print('qtiles:',qtiles, len(qtiles))
    if verbose: print('new_x:',new_x, len(new_x))
    if len (new_x) < 2:
        new_x = np.append(new_x,new_x[0] + 1e-9)
    if len(qtiles) < 2:
        qtiles = np.append(qtiles,qtiles[0] + 1e-9)
    if verbose: print('qtiles after:',qtiles,len(qtiles))
    if verbose: print('new_x after:',new_x,len(new_x))

    interp = interp1d(qtiles, new_x, fill_value="extrapolate")
    for q in [0.25, 0.5, 0.75] :
        ax2.plot([interp(q),interp(q)], [0, new_cy.max()],alpha=0.3,c="r", ls="dotted")
        ax2.text(interp(q), new_cy.max()*0.95, "{:2d}%".format(int(q*100)), 
            verticalalignment='bottom', alpha=0.3, color="r")

    mjd_h = (new_mjd - new_mjd[0])*24.
    interp = interp1d(mjd_h,new_x)
    x_ticks = interp( np.arange(0, mjd_h.max(), 1) )

    mjd_24 = (new_mjd - new_mjd[0])*24. + new_hour[0]
    ix = mjd_24 > 24
    mjd_24[ix] = mjd_24[ix]-24.0
    interp = interp1d(mjd_h, mjd_24)
    hours_labels = interp( np.arange(0, mjd_h.max(), 1) )
    labels = []
    for i in range(hours_labels.size) :
        labels.append( "{:.1f}".format(hours_labels[i]) )
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticks(x_ticks)
    ax3.set_xticklabels(labels)
    ix = slotNum == u_slotNum[0]
    year,month,day,hour,minute = utc_time_from_mjd(mjd[ix].min())
    ix = slotNum == u_slotNum[-1]
    eyear,emonth,eday,ehour,eminute = utc_time_from_mjd(mjd[ix].max())
    ax3.set_xlabel("UT hour, starting {}/{}/{}, ending {}/{}".format(
        year, month, day, emonth, eday))
    print "\t writing {}-slot-probabilities.png".format(data_dir+'/'+trigger_id)
    plt.savefig("{}-slot-probabilities.png".format(data_dir+'/'+trigger_id))



def cumulPlot2(trigger_id, data_dir) :
    from scipy.interpolate import interp1d
    ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)
    ix = np.argsort(prob)[::-1]
    ra,dec,id,prob,mjd,slotNum,dist = \
        ra[ix],dec[ix],id[ix],prob[ix],mjd[ix],slotNum[ix],dist[ix]  
    
    plt.close()
    fig, ax1 = plt.subplots() ; 
    ax2 = ax1.twinx()

    x = range(0,prob.size)
    ax1.scatter(x,prob,c="b",s=15)
    ax1.plot(x,prob,c="k",zorder=100)
    ax1.set_xlabel("hex count,  max to min")
    ax1.set_ylabel("LIGO probability in slot")
    ax1.set_ylim(0,prob.max()*1.1)
    ax1.grid(alpha=0.2)

    new_cy = np.cumsum(prob)
    ax2.scatter(x,new_cy,color="g",s=15, zorder=90)
    ax2.plot(x,new_cy,c="k",zorder=101)
    ax2.tick_params(axis='y', labelcolor="g")
    ax2.set_ylabel("cumulative LIGO probability",color="g")
    ax2.set_ylim(0,new_cy.max())
    fig.tight_layout()

    plt.title("Hexes ordered by probability")
    qtiles = new_cy/new_cy[-1]
    interp = interp1d(qtiles, x, fill_value="extrapolate")
    for q in [0.25, 0.5, 0.75] :
        ax2.plot([interp(q),interp(q)], [0, new_cy.max()],alpha=0.3,c="r", ls="dotted")
        ax2.text(interp(q), new_cy.max()*0.95, "{:2d}%".format(int(q*100)), 
            verticalalignment='bottom', alpha=0.3, color="r")
    print "\t writing {}-hex-probabilities.png".format(data_dir+'/'+trigger_id)
    plt.savefig("{}-hex-probabilities.png".format(data_dir+'/'+trigger_id))

def equalAreaPlot(figure,slot,simNumber,data_dir, title="") :
    import matplotlib.pyplot as plt
    from equalArea import mcplot
    from equalArea import mcbryde
    import insideDesFootprint

    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(data_dir, simNumber, slot)
    # x,y are the mcbryde projection of ra, dec
    # hx,hy are the mcbryde projection of ha, dec
    ra, dec = x, y

    # des footprint
    # plots the DES footprint on the maps
    desra, desdec = insideDesFootprint.getFootprintRaDec()
    desx, desy = mcbryde.mcbryde(desra, desdec)


    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-ligo-eq.png")
    print "making ",name,
    plt.clf();mcplot.plot(ra,dec,ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-maglim-eq.png")
    print name,
    plt.clf();mcplot.plot(ra,dec,maglim,vmin=17);
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-prob-eq.png")
    print name,
    plt.clf();mcplot.plot(ra,dec,prob)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-probXligo-eq.png")
    print name
    plt.clf();mcplot.plot(ra,dec,prob*ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)
    # return the number of plots made
    return 4 

# modify mcbryde to have alpha=center of plot
#   "slot" is roughly hour during the night at which to make plot
def observingPlot(figure, simNumber, slot, data_dir, nslots, n_hexes, camera, 
        allSky=False, gif_resolution=1.0,
        moonRa=-999, moonDec=-999, moonIllumination=0):
    import plotMapAndHex

    if n_hexes > 0 :
        # get the planned observations
        ra,dec,id,prob,mjd,slotNumbers,dist = obsSlots.readObservingRecord(simNumber, data_dir)
        ix = slotNumbers == slot
        if np.any(ix):
            the_mjd = mjd[ix][0]
            time = utcFromMjd(the_mjd, alt=False)
            time = "Shutter open " + time
        else :
            time = "Shutter closed " 
    else :
        time = "Shutter closed " 
        ra = -999; dec = -999
        slotNumbers = -999
    
    if allSky :
        title = "{} Slot {} {}".format(simNumber, slot, time)
    else :
        title = "{} Slot {}".format(simNumber, slot)


# this is useful to debug the plots
    #print "making plotMapAndHex.mapAndHex(figure, ", simNumber, ",", slot, ",", data_dir, ",", nslots, ",ra,dec,", camera, title,"allSky=",allSky,") "

    d=plotMapAndHex.mapAndHex(figure, simNumber, slot, data_dir, nslots, ra, dec, \
        camera, title, slots=slotNumbers, allSky=allSky, scale=gif_resolution,
        moonRa=moonRa, moonDec=moonDec, moonIllumination=moonIllumination)
    return d

def writeObservingRecord(slotsObserving, data_dir, gw_map_trigger, gw_map_control, gw_map_strategy) :
    trigger_id = gw_map_trigger.trigger_id
    just_sort_by_ra = gw_map_control.just_sort_by_ra
    max_number_of_hexes_to_do = gw_map_strategy.max_number_of_hexes_to_do

    name = os.path.join(data_dir, str(trigger_id) + "-ra-dec-id-prob-mjd-slot-dist.txt")
    ra,dec,id,prob,mjd,slotNum,islot = obsSlots.slotsObservingToNpArrays(slotsObserving)

    # let's optionally limit the number of hexes going into the file:
    #    of course, this assumes they are sorted by prob
    if max_number_of_hexes_to_do < 10000 :
        ra = ra[0:max_number_of_hexes_to_do]
        dec = dec[0:max_number_of_hexes_to_do]
        id = id[0:max_number_of_hexes_to_do]
        prob = prob[0:max_number_of_hexes_to_do]
        mjd = mjd[0:max_number_of_hexes_to_do]
        slotNum = slotNum[0:max_number_of_hexes_to_do]
        islot = islot[0:max_number_of_hexes_to_do]
        print "\t ========================================"
        print "\t number of hexes kept for jsons = {}".format( max_number_of_hexes_to_do)
        print "\t ========================================"
    
    if just_sort_by_ra :
        # rearrange inside the slots if desired
        ix = np.argsort(ra)
        ra,dec,id,prob,islot = \
            ra[ix],dec[ix],id[ix],prob[ix],islot [ix]

    map_distance = gw_map_trigger.ligo_dist
    nside = hp.npix2nside(map_distance.size)
    #ang2pix(nside,theta,phi,nest=False)
    pix_num = hp.ang2pix(nside,ra,dec, lonlat=True)
    dist = map_distance[pix_num]
    #fixing the distance
    #dist = 60.
    fd = open(name,'w')
    fd.write("# ra, dec, id, prob, mjd, slotNum, dist\n")
    unique_slots = np.unique(slotNum)
    for slot in unique_slots:
        ix = slot==slotNum
        iy = np.argsort(ra[ix])
        # sort by ra inside the slot
        if (np.any(ra[ix] > 90) and np.any(ra[ix] < -90) ) :
            dummy_ra = ra[ix]
            iy = dummy_ra < 0
            dummy_ra[iy] = dummy_ra[iy] + 360
            iy = np.argsort(dummy_ra)
        else :
            iy = np.argsort(ra[ix])
        for r,d,i,p,m,s,di in zip(
                ra[ix][iy], dec[ix][iy], id[ix][iy],
                prob[ix][iy], mjd[ix][iy], slotNum[ix][iy], dist[ix][iy]):
                #ra[ix], dec[ix], id[ix],
                #prob[ix], mjd[ix], slotNum[ix], dist[ix]):
            fd.write("{:.6f} {:.5f} {:s} {:.7f} {:.4f} {:.1f} {:.2f}\n".format(r,d,i,p,m,s,di))
    fd.close()
    return ra,dec,id,prob,mjd,slotNum,dist

