import numpy as np
import mags
import obsSlots

#
# calc produces information about a single nights observing.
#   the information about what to observe is in the *ra-dec-*txt file that strategy produces
#   the "ra_dec_file_dir" tells the code where that lives
#   tigger_id tells it what the name of the *ra-dec*.txt file is
# The expTime is used only in the limiting magnitude calculation, not in what to observe when
#

def calc(trigger_id, ra_dec_file_dir, mjd, expTime, filter, best=False, camera="decam") :
    # do dark siren stuff
    #ra,dec,prob,filter = np.genfromtxt(trigger_id,unpack=True)
    #filter = np.genfromtxt(trigger_id,unpack=True,dtype="str",usecols=3)
    #filter = "g"

    # get the hexes
    ra,dec,id,prob,obs_mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,ra_dec_file_dir)
    # find the highest prob hex
    best_ix = np.argmax(prob)

    # find night statistics
    night,sunset,sunrise = mags.findNightDuration(mjd, camera)
    night = np.float(night)*24.
    sunset = np.float(sunset)
    sunrise = np.float(sunrise)

    # will work every hour from sunset to sunrise
    mjd_list = np.arange( sunset, sunrise+1./24., 1./24.)
    
    # calculate the limiting magnitudes
    limit_mag = []
    best_limit_mag = []
    moon_sep = []
    moon_phase = []
    obs = mags.observed(ra,dec,prob, sunset, doMaps=False, verbose=False)
    for mjd in mjd_list :
        obs.resetTime(mjd)
        obs.limitMag(filter,exposure=expTime)
        limit_mag.append(obs.maglim)
        best_limit_mag.append(obs.maglim[best_ix])
        moon_sep.append(obs.moonSep[0]*360./2/np.pi)
        moon_phase.append(obs.moonPhase)
    limit_mag = np.vstack(limit_mag)
    ix = limit_mag < 0; limit_mag[ix] = 0
    moon_sep = np.array(moon_sep)
    moon_phase = np.array(moon_phase)
    # moon phase: where 0 = full, 90 equals half, and  180 = new
    # convert to %full
    moon_phase = ((180.-moon_phase)/180.)*100.
    bins = np.arange(21,25.5,0.5)

    # now print the answers
    print("{:15s}  {:10s}".format("MJD".rjust(15),"best".rjust(10)), end=' ')
    print(" {:10s}".format("moon sep".rjust(10)), end=' ')
    print(" {:10s}".format("moon phase".rjust(10)), end=' ')
    print("      hexes w/ limiting mag in ", end=' ')
    print("")
    print("{:15s}  {:10s}".format("days".rjust(15),"hex".rjust(10)), end=' ')
    print(" {:10s}".format("degrees".rjust(10)), end=' ')
    print("  {:10s}".format("% full ".rjust(10)), end=' ')
    for i in range(bins.size-1) : print("{:4.1f}-".format(bins[i]), end=' ')
    print("{:4.1f}".format(bins[i]), end=' ')
    print("")
    for i in range(0,mjd_list.size) :
        print("{:15.4f} ".format(mjd_list[i]), end=' ')
        print("{:10.2f} ".format(best_limit_mag[i]), end=' ')
        print("{:10.2f} ".format(moon_sep[i]), end=' ')
        print("{:10.0f}  ".format(moon_phase[i]), end=' ')
        counts = np.histogram(limit_mag[i], bins)
        #print "\n counts ", np.shape(counts), counts[0], counts[1]
        for j in range(bins.size-1) : 
            if counts[0][j] != 0 :
                print("{:3d}  ".format(counts[0][j]), end=' ')
            else :
                print("     ", end=' ')
        print("")
    return obs


