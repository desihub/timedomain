"""
Form the classes for handling the many different
inputs to the map making code.
"""
import fitsio
import numpy as np
import hp2np
import healpy as hp
import yaml

class trigger(object):
    """
    """
    def __init__(self, skymap, trigger_id, trigger_type, resolution, days_since_burst=0.) :
        """
        """
        self.skymap = skymap
        self.trigger_id = trigger_id
        self.trigger_type = trigger_type
        self.days_since_burst = days_since_burst

        hdr = fitsio.read_header(skymap,1)
        burst_mjd = np.float(hdr["mjd-obs"])
        try:
            distance = hdr["distmean"]
            diststd = hdr["diststd"]
        except:
            print('distmean was not in payload... ')
            raise Exception("now")
        self.distance  = distance
        self.diststd  = diststd
        self.burst_mjd = burst_mjd

        # ok, I'm going to declare that we want this routine to start at noon UT on JD of burst
        # As there is automatic sunrise, sunset calculations given MJD, we will pick midnight
        self.start_mjd = np.round(burst_mjd)

        if trigger_type == "bright" :
            named_trigger = "has remnant"
        else:
            named_trigger = "dark"

        print ""
        print "=================================================="
        print "=================  desgw-map ====================="
        print "=================================================="
        print "                                        since 2015"
        print "\ngw_map_trigger: {} map {}, {} at {:.0f} Mpc\n".format(
            trigger_id, skymap, named_trigger, distance)
        if trigger_type == "dark" :
            print('\t strategy=dark,  setting distance, dist_err to 100 Mpc, 30 Mpc')
            self.distance = 100.
            self.diststd = 30.

        ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)
        ligo_dist, ligo_dist_sig, ligo_dist_norm  = \
            distance*np.ones(ra.size), np.zeros(ra.size), np.zeros(ra.size)
        try :
            ligo_dist      = hp.read_map(skymap, field=1, verbose=False)
            ligo_dist_sig  = hp.read_map(skymap, field=2, verbose=False)
            ligo_dist_norm = hp.read_map(skymap, field=3, verbose=False)
            # there are infinite distance pixels in the map. deal with these
            ligo_dist[np.isinf(ligo_dist)] = 10000.

            # change resolution
            ligo_dist      = hp.ud_grade(ligo_dist, resolution)
            ligo_dist_sig  = hp.ud_grade(ligo_dist, resolution)
            ligo_dist_norm = hp.ud_grade(ligo_dist, resolution)
        except:
            print "\n\t !!!!!!!! ------- no distance information in skymap ------ !!!!!!!!\n"
    
        # JTA Hack
        # GW190425
        #print "HACK HACK HACK"
        #ix = ((ra > 150) | (( ra > -180) & (ra < -120)) ) & (dec > -10) & (dec < 40)
        #ix = np.invert(ix)
        #ligo[ix] = 0.0
        #print "HACK HACK HACK"
    
        # GW170217 hack JTA
        #ix = (ra > 0) & ( ra < 180) & (dec >= -30)
        #ix = np.invert(ix)
        #ligo[ix] = 0.0
        # GW170225 hack JTA
        #ix = (dec >= 2)
        #ligo[ix] = 0.0
        # GW170814 hack JTA
        #ix = (ra > -10) & ( ra < 60) & (dec < -20)
        #ix = np.invert(ix)
        #ligo[ix] = 0.0

        # GW170814 hack JTA                                                                             
        ix = (ra > 0) & ( ra < 45) & (dec < 0) & (dec > -40)
        ix = np.invert(ix)                                                                             
        ligo[ix] = 0.0

        self.ligo_ra = ra
        self.ligo_dec = dec
        self.ligo = ligo
        self.ligo_dist = ligo_dist
        self.ligo_dist_sig = ligo_dist_sig
        self.ligo_norm = ligo_dist_norm


# control the observations
class strategy(object) :
    """
    """
    def __init__(self, camera, exposure_list, filter_list, tiling_list, maxHexesPerSlot, 
            hoursAvailable, propid, max_number_of_hexes_to_do, kasen_fraction, use_teff):
        """
        """
        self.camera                    = camera
        self.exposure_list             = np.array(exposure_list)
        self.filter_list               = np.array(filter_list)
        self.tiling_list               = np.array(tiling_list)
        self.maxHexesPerSlot           = maxHexesPerSlot
        self.hoursAvailable            = hoursAvailable
        self.propid                    = propid
        self.max_number_of_hexes_to_do = max_number_of_hexes_to_do
        self.kasen_fraction            = kasen_fraction
        self.use_teff                  = use_teff
        if abs(use_teff - 1.0) > 0.01 :
           print "\t scaling summed exposure time by teff {}".format(use_teff)

        working_filter            = self.filter_list[0]
        ix                        = self.filter_list == working_filter
        sum_exposure_time         = self.exposure_list[ix].sum() * use_teff
        self.summed_exposure_time = sum_exposure_time
        self.working_filter       = working_filter

        if camera == "decam" :
            self.overhead =  30. # seconds
            self.area_per_hex = 3.0 # sq-degrees
        elif camera == "hsc" :
            self.overhead = 20.
            self.area_per_hex = 1.5
        elif camera == "desi" :
            self.overhead = 300.
            self.area_per_hex = 9.8
        else: raise Exception("camera {} not handled".format(camera))



# control the code technically
class control(object):
    """
    """
    def __init__(self, resolution, data_dir, debug=False, allSky=False, centeredSky=True,
            snarf_mi_maps=False, mi_map_dir="/data/des41.a/data/desgw/O3FULL/Main-Injector/OUTPUT/O3REAL/",
            gif_resolution = 1.0 ) :
        """
        """
        # healpix map resolution
        self.resolution = resolution
        self.gif_resolution = gif_resolution
        self.debug = debug
        self.this_tiling = []
        self.reject_hexes= []
        self.allSky = allSky
        self.centeredSky = centeredSky
        self.datadir = data_dir
        # find hexes starting with slot start_slot. Carries on to jsons
        self.start_slot = -1
        # find hexes for the do_nslots starting with slot start_slot. Carries on to jsons
        self.do_nslots = -1
        self.just_sort_by_ra = False # if true, take the calculated hexes and just sort
                                    # the time to be observed as increasing ra
        # should i bypass the mapmaking by copying over from the MI directories?
        self.snarf_mi_maps = snarf_mi_maps
        # what directory to holds the existing maps? Copy over to output_dir
        self.mi_map_dir = mi_map_dir 
 

# save intermediate results
class results(object):
    """
    """
    def __init__(self) :
        """
        """
        self.probability_per_slot = False
        self.time_of_slot = False
        self.isdark = False
        self.made_maps_list = False
        self.slotDuration = False
        self.hoursPerNight = False
        self.n_slots = False
        self.first_slot = False
        self.best_slot = False
        self.slot_numbers = False
        self.n_hexes = False
        self.moonRa = False
        self.moonDec = False

