import numpy as np
import os
import yaml
import getHexObservations
import gw_map_configure
import healpy as hp


#
#   the real recycler will never use 
#       do_make_maps, do_make_hexes, do_make_jsons, do_make_gifs
#       start_slot, do_nslots, mi_map_dir, snarf_mi_maps
#   and the defaults in gw_map_configure.py are correct for it.
#
def recycle (trigger_id, skymap, trigger_type,
        outputDir, propid = -1, do_make_maps = -1, do_make_hexes = -1,
        do_make_jsons = -1, do_make_gifs = -1, 
        start_slot = -1, do_nslots= -1,  mi_map_dir = "./", 
        snarf_mi_maps = False) :

    with open("maininjector.yaml","r") as f:
        config = yaml.safe_load(f); 

    # debug
    debug = config["debug"]    
    # camera
    camera   = config["camera"]

    # control the code running
    if do_make_maps == -1:
        do_make_maps   = config["do_make_maps"]
    if do_make_hexes == -1:
        do_make_hexes  = config["do_make_hexes"]
    if do_make_jsons == -1:
        do_make_jsons  = config["do_make_jsons"]
    if do_make_gifs == -1:
        do_make_gifs   = config["do_make_gifs"]
    allSky = config["allSkyGif"]
    centeredSky = config["centeredSkyGif"]

    # healpix map resolution
    resolution     = config["resolution"]
    # gif pixel resolution
    gif_resolution     = config["gif_resolution"]


    # same day?
    days_since_burst = config["days_since_burst"]

    # strategy
    kasen_fraction         = config["kasen_fraction"]
    hoursAvailable         = config["hoursAvailable"]
    use_teff               = config["use_teff"]
    exposure_length_bright  = config["exposure_length_bright"]
    filter_list_bright     = config["exposure_filter_bright"]
    exposure_tiling_bright = config["exposure_tiling_bright"]
    maxHexesPerSlot_bright = config["maxHexesPerSlot_bright"]
    exposure_length_dark  = config["exposure_length_dark"]
    filter_list_dark      = config["exposure_filter_dark"]
    exposure_tiling_dark  = config["exposure_tiling_dark"]
    maxHexesPerSlot_dark  = config["maxHexesPerSlot_dark"]
    max_number_of_hexes_to_do = config["max_number_of_hexes_to_do"]

    # configure strategy for the event type
    if trigger_type == "bright" :
        exposure_length      = exposure_length_bright
        filter_list          = filter_list_bright
        tiling_list          = exposure_tiling_bright
        maxHexesPerSlot      = maxHexesPerSlot_bright
    elif trigger_type == "dark" :
        exposure_length      = exposure_length_dark
        filter_list          = filter_list_dark 
        tiling_list          = exposure_tiling_dark
        maxHexesPerSlot      = maxHexesPerSlot_dark
    else :
        raise Exception(
            "trigger_type={}  ! Can only compute bright or dark".format(trigger_type))
    exposure_length   = np.array(exposure_length)
    maxHexesPerSlot = np.int(np.round(maxHexesPerSlot/np.size(tiling_list)))

    gw_map_control  = gw_map_configure.control( resolution, outputDir, debug, 
        snarf_mi_maps=snarf_mi_maps, mi_map_dir = mi_map_dir,
        allSky=allSky, centeredSky=centeredSky, 
        gif_resolution = gif_resolution)
    gw_map_trigger  = gw_map_configure.trigger( skymap, trigger_id, trigger_type, 
        resolution, days_since_burst=days_since_burst)
    gw_map_strategy = gw_map_configure.strategy( camera, exposure_length, 
        filter_list, tiling_list, maxHexesPerSlot, hoursAvailable, propid,
        max_number_of_hexes_to_do, kasen_fraction, use_teff)
    gw_map_results = gw_map_configure.results()

    if not os.path.exists(outputDir): os.makedirs(outputDir)
    
    if do_make_maps :
        # make the computationally expensive maps of everything
        getHexObservations.make_maps( 
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)

    if do_make_hexes :
        # compute the best observations
        getHexObservations.make_hexes( 
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results,
            start_slot = start_slot, do_nslots= do_nslots)
        # if start_slot = -1, do_nslots = -1, then do whole night, as if called like:
        #    getHexObservations.make_hexes( 
        #        gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)
        # the web page maker version of main injector should default to -1, -1

    if do_make_jsons :
        # make the jsons 
        getHexObservations.make_jsons( gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)

    if do_make_gifs :
        getHexObservations.makeGifs( gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)

    return 


#
#  so what hexes did we do, given the json files?
def gethexIDfromJson (dir = "/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_night1_json/",
        verbose=True) :
    import json; import glob
    import os
    import hexalate
    cwd = os.getcwd(); print cwd
    os.chdir(dir)
    files = glob.glob("*.json")
    nfiles=len(files);
    ra = np.array([])
    dec = np.array([])
    slot = np.array([])
    for n in range(nfiles) :
        print files[n]
        fd=open(files[n],"r"); data=json.load(fd);fd.close()
        for i in range(len(data)): 
            ra = np.append(ra, data[i]["RA"])
            dec = np.append(dec, data[i]["dec"])
            slot = np.append(slot, files[n][9:11])
            if verbose: print data[i]["RA"],data[i]["dec"],files[n][9:11]

    os.chdir(cwd)
    id = hexalate.getHexId(ra,dec)
    return ra,dec,id, slot
def gethexIDfromDB() :
    import hexalate
    file="/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_exp_table_night1.txt"
    file="/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_exp_table_night1_2.txt" 
    ra,dec = np.genfromtxt(file, unpack=True,usecols=(1,2)); 
    ii=ra>180; ra[ii] = ra[ii]-360
    id = hexalate.getHexId(ra,dec)
    id = np.unique(id)
    return id
#
# shall we run the hexes from a ra-dec-id-prob-mjd-slot-dist.tx file?
def hexID (file = "../Jan4-2017-event/GW170104-ra-dec-id-prob-mjd-slot-dist.txt") :
    ra,dec = np.genfromtxt(file, unpack=True,usecols=(0,1) )
    id = np.genfromtxt(file, unpack=True,usecols=(2), dtype = "str" )
    # for GW170104; doing this cut does the obs right, not doing makes the gif cover better area
    ix = ra < 100; ra = ra[ix]; dec = dec[ix]; id = id[ix]
    return ra,dec,id
