import yaml
import argparse
import datetime
from astropy.time import Time
import os


### create todays date for season
today = datetime.date.today()
start_date = Time(str(today))
start_date.format = 'mjd'
end_date = Time(str(datetime.timedelta(days=30) + today))
end_date.format = 'mjd'
time_of_running_this=datetime.datetime.utcnow()
utc = Time(str(time_of_running_this))
utc.format = 'mjd'


###Function to make yaml for Recycler

def makeYaml(camera, res, propid, start_season=start_date.value, end_season=end_date.value,sendEmail=False, sendTexts=False, debug=False): #camera string name, boolean, boolean
    data=dict(trigger_path  = './OUTPUT/REAL/',
    real_or_sim= 'real', #must switch between 'sim' and 'real'

    # camera
    camera = camera, ##make this an argument with decam as the default

    start_of_season = start_season, ## current date?
    end_of_season = end_season, ## default to 30 days later, some user interaction

    sendEmailsToEveryone = sendEmail, ## argument
    primaryemails = ['djbrout@gmail.com'], ##add others?
    sendtexts = sendTexts, ## should be true?? do we need texts/emails after the first notifications

    ishalfnight=  False,
    isfirsthalf= False,
    issecondhalf= False,

    skipjson= False,

    probHasNS_threshold= .5,

    forceProbHasNS= True,
    probHasNS= 1., ## ??? should this be changable?

    events_observed = ['GW150914' , 'GW151226'], ## fill in from gracedb

    # Static observing params
    area_per_hex = 3.0, # sq deg
    overhead = 30., #sec

    # Optimization params
    resolution= res, #options:64, 128, 256. for map making


    distance = 60.,
    force_distance = False,

    skipAll = False,
    skipPlots= False,

    default_map_name = 'bayestar.fits', ## argument?

    # economics analysis
    time_budget_for_NS = 20., # hours, assuming 5 nights and 10hrs/night, 3/5 for NS
    time_budget_for_BH = 1.5, # hours, assuming 5 nights and 10hrs/night, 2/5 for BH
    hours_lost_to_weather_for_NS = 0, ## changable??
    hours_lost_to_weather_for_BH = 0, ## changable??
    rate_of_bh_in_O2= 20.0,  # n/yr, BH merger triggers per year in observing run 2
    rate_of_ns_in_O2=  2.0, # n/yr, NS merger triggers per year in observing run 2

    # epoch structures. For each epoch the begin date must be set
    # NS strategy
    nepochs_NS = 4,
    epoch1_NS  = 0, # start of epoch 1 in days since burst
    epoch2_NS  = 2, # start of epoch 2 in days since burst
    epoch3_NS  = 4, # start of epoch 3 in days since burst
    enddate_NS = 10, # termination time of observations in days since burst

    # BH strategy
    nepochs_BH = 3,
    epoch1_BH  = 0, # start of epoch 1 in days since burst
    epoch2_BH  = 1, # start of epoch 2 in days since burst
    epoch3_BH  = 14, # start of epoch 3 in days since burst
    enddate_BH = 18, # termination time of observations in days since burst

    # NS strategy
    exposure_length_NS = [ 60., 90. ], #sec
    exposure_filter_NS = [ 'i', 'z' ], ## arguments??
    maxHexesPerSlot_NS = 6,
    exposure_tiling_NS = [ 9, 9 ],
    # BH strategy
    exposure_length_BH = [ 90., ], #sec
    exposure_filter_BH = [ 'i', ],
    exposure_tiling_BH = [ 10, ],
    maxHexesPerSlot_BH = 18,  # related to maxHexesPerSlot_NS by x3 for x3 less images/hex, so same amount of time
    
    propid = propid, #to be used in postprocessing

    ####Only used for running recycler.py standalone!####

    trigger_id = "G297595", ## argument -- where is this generated?
    wrap_all_triggers = False,

    mjd = 57974.3530295, ## what date is this?? one of these is the date of the ligo trigger 
#    force_mjd = False, #if its for making the map could be the start time for making the map

    recycler_mjd = 57979.5, ## what date is this?? # 
#    force_recycler_mjd = True,

    skymap_filename = './OUTPUT/REAL/G297595/bayestar/bayestar.fits.gz',
    
    debug = debug)
    
    ### save output ###
    yamlName='recycler_test.yaml'
    
    with open(yamlName, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
    yamlName.close()
        
    return yamlName


#Change this line to improve your experience.
#yamlName= makeYaml(camera="decam", res="126", propid="x") 

#os.system("sed -i -e 's/false/False/g' "+yamlName)
#os.system("sed -i -e 's/true/True/g' "+yamlName)
