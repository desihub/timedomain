import numpy as np
import os
import warnings

# ========== do simple calculations on how to divide the night
#
#   one of the questions is how many hours to devote to observations
#       hoursAvailable,  another is the slot duration
#
# if the 1% cut isn't in place in mapsAtTimeT.oneDayOfTotalProbability
# then one expects circa 40-45 maps as there is about 2/hour
#   with the 1% cut, then one expects far fewer. Perhaps zero.

# May 2019
# I added start_slot and do_nslots, in the name of running from a start for a 
#   number of slots. This brought out the mapping from slotNumber and hexData to 
#   the actual slot number, done in the code at def observing line : map_i = i + mapZero
# The issue is that start_slot and do_nslot are in the mapped slot number, so many
#   of the slots don't exist, and much of the infrastructure loops over all slots
#   I had to put in 6 places a "if do_nslots > -1; if i+mapZero < start_slot" technology
#   (marked by JTA)  to deal with this
#   

#===============================================================================
# These calculations used to be spread over hither and yon.
# Bring them together.
#
#    deltaTime = 0.0223  => 32 minute slots  for izzi 90sec exp, 4 hex/slot
#    deltaTime = 0.0446  => 62 minute slots  for izzi 90sec exp, 8 hex/slot
#    deltaTime = 0.0417  => 60 minute slots  for izz 90sec exp, 10 hex/slot
#    deltaTime = 0.0208  => 30 minute slots  for izz 90sec exp,  5 hex/slot
#       The first is flexible to observing conditions,
#       while the second is x2 faster as there are less maps to hexalate
#       The third shows a different tiling/filter complement
#
# So-
#
# Let us redfine this: a slot is the length of time it takes
# to do 6 hexes to completion. That is usually somewhere between 30 minutes
# and one hour, so close to the original defintion, and by force is an
# even number of hexes. Ok. Use n=6 for the forcing definition
#
# This is based on the current "NS" observing strategy, izz at 90s each.
# If instead we move to a BH strategy that is 1xi at 90s, then
# one could still use 6 hexes/slot but that means the number of slots rises by x3.
# Instead, what if we used 18 hexes/slot for the BH case.
#
# Ok, then, code:
#
def slotCalculations(mjd, exposure_lengths, tiling_list, overhead, hexesPerSlot = 6, camera="decam") :
    from getHexObservations import hoursPerNight
    slot_duration = slotDuration(exposure_lengths, tiling_list, overhead, hexesPerSlot) 
    hoursAvailable = hoursPerNight(mjd, camera)
    answers = dict()
    answers["slotDuration"] = slot_duration
    answers["hoursPerNight"] = hoursAvailable
    return answers

def slotDuration(exposure_lengths, tiling_list, overhead, hexesPerSlot = 6) :
    tot_exptime = np.size(tiling_list)*(np.array(overhead)+np.array(exposure_lengths)).sum()
    slot_time = tot_exptime*hexesPerSlot
    slot_duration = slot_time/60. ;# in minutes
    return slot_duration

# find the number of slots per night
def findNSlots(hoursAvailable, slotDuration=32.) :
    verbose = 0
    if verbose:
        print hoursAvailable
        print hoursAvailable*60./slotDuration, round(hoursAvailable*60./slotDuration)
        print int(round(hoursAvailable*60./slotDuration))
    nslots = int(round(hoursAvailable*60./slotDuration))   ;# 32 minutes/slot
    return nslots

# ok, I designed observing for the case
#   where the number of slots in a night
#   was equal to the number of maps made.
# This is unrealistic for two reasons:
#   the night could be longer than the time desgw wishes to allocate
#   the object could be up for less than the time desgw could allocate
# so:
    # possibilities
    # n_maps = n_slots     design case
    # n_maps < n_slots     set n_slots = n_maps
    # n_maps > n_slots     pick highest contiguous set of n_slot maps
# the first two are easy. 
# the third runs into the naming convention of the maps
def findStartMap ( probs, times, n_slots ) :
    n_maps = times.size
    mapNums = np.arange(0,n_maps)
    # sum the probability in each allowed nslot range
    n_mapsToTProb = np.array([])
    n_mapStart = np.array([])
    for map in range(0,n_maps-n_slots) :
        ix = (mapNums >= map) & (mapNums < map+n_slots)
        n_mapsToTProb = np.append(n_mapsToTProb, probs[ix].sum() )
        n_mapStart = np.append(n_mapStart, map)
    minToMax = np.argsort(n_mapsToTProb)
    bestStart = n_mapStart[minToMax[-1]]
    bestStart = int(bestStart)
    return bestStart


# Load all n -hexVals files,
#   pick the highest probability one, 
#   put it into one of the n slots, unless that slots is maxed out
#   remove that hex from all the hexVals lists
#   do it again, untill all n time slots are full.
#       maxHexesPerSlot=4 comes from 32 minute duration slots
#       and 8 minutes/hex (izzi 2 min/image)
#
#   if we do zzi at 2 mins/image then 4 min/hex + 2 min/hex2 = 6 mins
#   call it 60 minute slots  and 10 hexes/slot
def observing(sim, nslots, data_dir, 
        maxHexesPerSlot = 4, mapZero = 0, do_nslots = -1, start_slot=-1, verbose=0) :
    # prep the observing lists
    observingSlots = np.arange(0,nslots)
    slotsObserving = dict()
    slotsObserving["nslots"] = nslots
    slotsObserving["mapZero"] = mapZero
    slotsObserving["do_nslots"] = do_nslots
    slotsObserving["start_slot"] = start_slot
    for i in observingSlots :
        slotsObserving[i] = 0
        slotsObserving[i,"ra"]   = np.array([])
        slotsObserving[i,"dec"]  = np.array([])
        slotsObserving[i,"id"]  = np.array([])
        slotsObserving[i,"prob"] = np.array([])
        slotsObserving[i,"mjd"] = np.array([])
        slotsObserving[i,"slotNum"] = np.array([]) ;#non-zero prob slots
        slotsObserving[i,"islot"] = np.array([]) ;# range(0,nslots)

    # read in the hexelated probability data
    hexData = dict()
    for i in observingSlots :
        map_i = i + mapZero
# JTA 1
        if start_slot > -1 and do_nslots > -1 :
            if map_i < start_slot or map_i >= start_slot+do_nslots : continue
        raHexen, decHexen, idHexen, hexVal, rank, mjd, slotNum = \
           loadHexalatedProbabilities( sim, map_i, data_dir)
        islot = i*np.ones(raHexen.size)
        print "\t", map_i, "map size= {:10d};".format(raHexen.size), 

        impossible = 1e-5
        impossible = 1e-7
        impossible = 1e-8
        ix = np.nonzero(hexVal < impossible)
        raHexen, decHexen, idHexen, hexVal, mjd, slotNum, islot  = \
            np.delete(raHexen, ix), \
            np.delete(decHexen, ix), \
            np.delete(idHexen, ix), \
            np.delete(hexVal, ix), \
            np.delete(mjd, ix) , \
            np.delete(slotNum, ix), \
            np.delete(islot, ix)
        print " n hexes w/ >{} probability=".format(str(impossible)),
        print "{:4d};".format(raHexen.size),
        print "  in slot,all_possible_hexes sum prob= {:7.4f} %".format( 100*hexVal.sum())
        hexData[i] = raHexen, decHexen, idHexen, hexVal, mjd, slotNum, islot
        #print i, np.sort(hexVal[0:10]), hexVal.sum(), 100.*hexVal.sum(),"%"

    # start the search for all max probabilities
    # we'll assume the list is less than 40,000 long, the n-sq-degrees/sky
    for n in range(0,40000) :
        # search for a single max probabilities
        maxRa, maxDec, maxId, maxProb, maxMjd, maxSlotNum, maxIslot  = \
            findMaxProbOfAllHexes(hexData, slotsObserving, observingSlots, n, verbose) 
        maxData = maxRa,maxDec,maxId, maxProb,maxMjd,maxSlotNum, maxIslot

        # we've found the maximum probability on the lists, 
        # so add it to the obs lists # unless not possible. 
        # If the latter, delete it from that slot
        slot = maxIslot
        # if slot is -1, then no max prob found
        if slot > -1 :  
            if slotsObserving[slot] < maxHexesPerSlot : 
                # it is possible to make the observation, 
                # put it onto the observing lists
                slotsObserving = addObsToSlot (slotsObserving, maxData, slot)
                if verbose >= 1: print n, "slot of max:",slot
            else :
                # but if this slot of observing is full, it is not possible 
                # to make the observation,
                # so move on AFTER deleting it from the list
                hexData = deleteHexFromSlot (hexData, slot, maxProb) 
        #if verbose >= 2: 
        #   if n > 7: raise Exception("jack")
    
        # perform the necessary bookkeeping, 
        # eliminating this hex from future observing
        hexData = deleteHexFromAllSlots (
            hexData, slotsObserving, observingSlots, maxRa, maxDec, verbose, n) 

        # do some summary statistics
        sumHexes = 0
        sumObs = 0
        for i in range(0,nslots) :
# JTA 6
            if start_slot > -1 and do_nslots > -1 :
                if i+mapZero < start_slot or i+mapZero >= start_slot+do_nslots : continue
            sumHexes += hexData[i][0].size
            sumObs += slotsObserving[i]

        if verbose >= 2: 
            print "sumHexes =", sumHexes, 
            print "   slots left=", len(observingSlots),
            print "   slots=",observingSlots,
            print "   n_obs=",
            for i in observingSlots:
                print slotsObserving[i],
            print "   sum prob= ",
            for i in observingSlots:
                print " {:8.6f}".format( slotsObserving[i,"prob"].sum()) ,
            print ""

        # eliminate full observing slots
        observingSlots = eliminateFullObservingSlots(\
            hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) 

        # check and see if we are done
        # two conditions: observing is full, or candidates empty
        if (len(observingSlots)==0) | (sumHexes == 0) :
            print "\t======================================== "
            if verbose >= 1: 
                print "n slots =", len(observingSlots)," == 0?"
                print "sumHexes = ", sumHexes, "==? 0"
            print "\tnumber of hexes possible to observe = ", sumObs
            print "\t======================================== "
            return slotsObserving 

        # otherwise, go back and do it again
    
    # we've done everything on the lists, we can observe it all,
    # return this great success that will never be reached.
    return slotsObserving
    
#
# examine the statistics of the observing lists
#
def observingStats( slotsObserving, mapZero=0, do_nslots=-1, start_slot=-1 ) :
    nslots = slotsObserving["nslots"]
    for i in range(0,nslots) :
# JTA 2
        if start_slot > -1 and do_nslots > -1 :
            if i+mapZero < start_slot or i+mapZero >= start_slot+do_nslots : continue
        slot_sum_prob = 100*slotsObserving[i,"prob"].sum() 
        if slot_sum_prob > 1e-7 :
            print "\t {:2d}".format(i+mapZero), 
            #print "slotnum={} ".format( slotsObserving[i,"slotNum"]),
            print "n hexes= {}".format( slotsObserving[i,"ra"].size), 
            print "  sum prob= {:7.4f} %".format( 100*slotsObserving[i,"prob"].sum())
    ra,dec,id,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 

    print "\tobservingStats:  ",
    print "observable prob_tot = {:.1f}%".format(100.*prob.sum()),
    print "\t   (from the prob in each slot, summed)"
    return ra,dec,id,prob,mjd,slotNum,islot

def observingStatsFromRaDecFile( trigger_id, data_dir, slotsObserving, mapZero=0, do_nslots=-1, start_slot=-1 ) :
    nslots = slotsObserving["nslots"]
    ra,dec,id,prob,mjd,slotNum,dist = readObservingRecord(trigger_id,data_dir)
    for i in range(0,nslots) :
# JTA 2
        if start_slot > -1 and do_nslots > -1 :
            if i+mapZero < start_slot or i+mapZero >= start_slot+do_nslots : continue
        ix, = np.where(slotNum == i)
        if not (prob.size == 1 and prob.sum() == 0) and ix.size > 0 :
            slot_sum_prob = 100*prob[ix].sum() 
            if slot_sum_prob > 1e-7 :
                print "\t {:2d}".format(i+mapZero), 
                print "n hexes= {}".format( ix.size ),
                print "  sum prob= {:7.4f} %".format( slot_sum_prob )
    print "\tobservingStats:  ",
    print "observable prob_tot = {:.1f}%".format(100.*prob.sum()),
    print "\t   (from the prob in each slot, summed)"
    return ra,dec,id,prob,mjd,slotNum,dist


#  ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)
def readObservingRecord(trigger_id, data_dir) :
    import os
    name = os.path.join(data_dir, str(trigger_id) + "-ra-dec-id-prob-mjd-slot-dist.txt")
    try :
        ra,dec,prob,mjd,slotNum,dist = \
            np.genfromtxt(name,unpack=True,comments="#",usecols=(0,1,3,4,5,6))
        id = np.genfromtxt(name,unpack=True,comments="#", usecols=(2),dtype="str")
    except :
        ra,dec,id,prob,mjd,slotNum,dist = \
            np.array([0,]),np.array([0,]),np.array(["0",]), \
            np.array([0,]),np.array([0,]),np.array([0,]),np.array([0,])
    return ra,dec,id,prob,mjd,slotNum,dist

def slotsObservingToNpArrays(slotsObserving) :

    do_nslots  = slotsObserving["do_nslots"]
    start_slot = slotsObserving["start_slot"]
    nslots     = slotsObserving["nslots"]
    mapZero    = slotsObserving["mapZero"]

    nslots = slotsObserving["nslots"]
    ra = np.array([])
    dec = np.array([])
    id = np.array([])
    prob = np.array([])
    mjd = np.array([])
    slotNum = np.array([])
    islot = np.array([])
    for i in range(0,nslots) :
# JTA 7
        if start_slot > -1 and do_nslots > -1 :
            if i+mapZero < start_slot or i+mapZero >= start_slot+do_nslots : continue
        ra = np.append(ra, slotsObserving[i,"ra"])
        dec = np.append(dec, slotsObserving[i,"dec"])
        id = np.append(id, slotsObserving[i,"id"])
        prob = np.append(prob, slotsObserving[i,"prob"])
        mjd = np.append(mjd, slotsObserving[i,"mjd"])
        slotNum = np.append(slotNum, slotsObserving[i,"slotNum"])
        islot = np.append(islot, slotsObserving[i,"islot"])
    return ra,dec,id,prob,mjd,slotNum,islot
    

#
# search for the single highest probability hex over all of the possible hexes
# in the hexData slots 
#  slotsObserving is a dictionary, observingSlots is a np array of slot numbers#
def findMaxProbOfAllHexes(hexData, slotsObserving, observingSlots, n="", verbose = 0) :
    do_nslots  = slotsObserving["do_nslots"]
    start_slot = slotsObserving["start_slot"]
    nslots     = slotsObserving["nslots"]
    mapZero    = slotsObserving["mapZero"]

    maxProb = -1
    for i in observingSlots :
# JTA 3
        if start_slot > -1 and do_nslots > -1 :
            if i+mapZero < start_slot or i+mapZero >= start_slot+do_nslots : continue
        data = hexData[i]
        hexRa     = data[0]
        hexDec    = data[1]
        hexId     = data[2]
        hexVal    = data[3]
        hexMjd    = data[4]
        hexMyslot = data[5]
        if hexVal.size == 0: continue
        if verbose >= 2: 
            if i == 2: print n,"====",i, "hexSize =",hexRa.size
        # now check for max prob
        newProb = hexVal.max()
        if verbose >= 4: print n,i, maxProb, ">?", newProb, "     n=",hexVal.size
        if newProb > maxProb :
            if verbose >= 1: print n,"==== new max", i, "       ",newProb , ">", maxProb
            ix = np.where(hexVal == newProb)
            maxRa     = hexRa[ix]
            maxDec    = hexDec[ix]
            maxId     = hexId[ix]
            maxVal    = hexVal[ix]
            maxMjd    = hexMjd[ix]
            maxProb   = newProb
            maxSlot   = hexMyslot[ix]
            islot = i
    #print "observingSlots", observingSlots
    if maxProb == -1 : 
        maxRa, maxDec, maxId, maxVal, maxMjd, maxSlot, islot = \
            -1,-1,-1,-1,-1,-1,-1
        #raise Exception("no max probability found")
    try :
        if len(maxRa) > 1 :
            maxRa, maxDec, maxId, maxVal, maxMjd, maxSlot, islot = \
            maxRa[0], maxDec[0], maxId[0], maxVal[0], maxMjd[0], maxSlot[0], islot
    except :
        pass
    return maxRa, maxDec, maxId, maxVal, maxMjd, maxSlot, islot

# we've found a hex,slot that can be observed so add it the the observing lists
def addObsToSlot (slotsObserving, maxData, slot) :
    maxRa  = maxData[0]
    maxDec = maxData[1]
    maxId  = maxData[2]
    maxVal = maxData[3]
    maxMjd = maxData[4]
    maxSlotNum = maxData[5]
    maxIslot = maxData[6]
    slotsObserving[slot,"ra"]   =  np.append(slotsObserving[slot,"ra"], maxRa)
    slotsObserving[slot,"dec"]  =  np.append(slotsObserving[slot,"dec"], maxDec)
    slotsObserving[slot,"id"]   =  np.append(slotsObserving[slot,"id"], maxId)
    slotsObserving[slot,"prob"] =  np.append(slotsObserving[slot,"prob"], maxVal)
    slotsObserving[slot,"mjd"]  =  np.append(slotsObserving[slot,"mjd"], maxMjd)
    slotsObserving[slot,"slotNum"]   =  np.append(slotsObserving[slot,"slotNum"], maxSlotNum)
    slotsObserving[slot,"islot"] =  np.append(slotsObserving[slot,"islot"], maxIslot)
    slotsObserving[slot] += 1
    return slotsObserving
# there can be no more observing in this slot, so this hex,slotj
# is impossible, delete it from the list.
def deleteHexFromSlot (hexData, slot, maxProb) :
    hexRa, hexDec, hexId, hexVal, hexMjd, hexSlotNum, hexIslot = hexData[slot] 
    ix = np.nonzero(hexVal == maxProb) 
    hexRa  = np.delete(hexRa, ix)
    hexDec = np.delete(hexDec, ix)
    hexId = np.delete(hexId, ix)
    hexVal = np.delete(hexVal, ix)
    hexMjd = np.delete(hexMjd, ix)
    hexSlotNum = np.delete(hexSlotNum, ix)
    hexIslot   = np.delete(hexIslot, ix)
    hexData[slot] = hexRa, hexDec, hexId, hexVal, hexMjd, hexSlotNum, hexIslot
    return hexData
# the hex,slot has made it onto an observing list, so remove the hex
# from all hex lists ( rm hex,*)
def deleteHexFromAllSlots (hexData, slotsObserving, observingSlots, maxRa, maxDec, verbose=0, n="") :
    warnings.filterwarnings("ignore")
    do_nslots  = slotsObserving["do_nslots"]
    start_slot = slotsObserving["start_slot"]
    nslots     = slotsObserving["nslots"]
    mapZero    = slotsObserving["mapZero"]

    for i in observingSlots:
# JTA 4
        if start_slot > -1 and do_nslots > -1 :
            if i+mapZero < start_slot or i+mapZero >= start_slot+do_nslots : continue
        data = hexData[i]
        data = hexData[i]
        hexRa  = data[0]
        hexDec = data[1]
        hexId  = data[2]
        hexVal = data[3]
        hexMjd = data[4]
        hexSlotNum = data[5]
        hexIslot   = data[6]
        #print len(hexRa)
        #print len(maxRa)
        try:
            ix = np.nonzero((hexRa == maxRa) & (hexDec == maxDec))
        except:
            print "len hexRA hexDec "+str(len(hexRa))+" "+str(len(hexDec))
            print "len maxRA maxDec "+str(len(maxRa))+" "+str(len(maxDec))
            maxRa = max(maxRa)
            maxDec = max(maxDec)
            ix = np.nonzero((hexRa == maxRa) & (hexDec == maxDec))
            print i
            print ix              
#        print "ag compare i "+str(i)
#        print ix
        if verbose >=4 : print ix, hexRa, maxRa
        if verbose >=2 : 
            ixs = np.shape(ix)[1]
            print n,"bookkeeping",i,"  nHex=",hexRa.size, "   ix.size", ixs, 
        hexRa  = np.delete(hexRa, ix)
        hexDec = np.delete(hexDec, ix)
        hexId  = np.delete(hexId, ix)
        hexVal = np.delete(hexVal, ix)
        hexMjd = np.delete(hexMjd, ix)
        hexSlotNum = np.delete(hexSlotNum, ix)
        hexIslot   = np.delete(hexIslot, ix)
        hexData[i] = hexRa, hexDec, hexId, hexVal, hexMjd, hexSlotNum, hexIslot

        if verbose >= 2: 
            print "\t after delete nHex=",hexRa.size,
            if ixs > 0 : print "index = ",ix[0][0]
            else: print ""
    return hexData

# it is useful to remove full observing slots from further processing,
# though they survive to be returned in the final lists
def eliminateFullObservingSlots(
        hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) :
    do_nslots = slotsObserving["do_nslots"]
    start_slot = slotsObserving["start_slot"]
    nslots = slotsObserving["nslots"]
    mapZero = slotsObserving["mapZero"]

    full = []
    for i in range(0,len(observingSlots)):
        # either we've observed as much as we can, or there is nothing to see
        slot = observingSlots[i]
# JTA 5
        if start_slot > -1 and do_nslots > -1 :
            if slot+mapZero < start_slot or slot+mapZero >= start_slot+do_nslots : continue
        if (slotsObserving[slot] >= maxHexesPerSlot) | (hexData[slot][0].size ==0): 
            full.append(i)
    if verbose >= 2: 
        print "hiding full observing slot ", 
        for f in full: print observingSlots[f],
        print ""
    observingSlots = np.delete(observingSlots, full)
    return observingSlots


def maxProbabilitySlot(prob,slotNumbers) :
    # find slot with the maximum probability
    maxProb = -1; maxProb_slot = -1
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        sumProb = prob[ix].sum()
        if sumProb > maxProb :
            maxProb = sumProb
            maxProb_slot = slot
    maxProb_slot = np.int(maxProb_slot)
    return maxProb_slot

#===================================================================
#
# Read in all of the hexalated probability files
#   hex inclusion and removal depends on hex id
#
def loadHexalatedProbabilities(sim, slot, data_dir) :
    import hexalate
    nameStem = os.path.join(data_dir, str(sim) + "-{}".format(str(slot)))
    name = nameStem + "-hexVals.txt"
    if os.path.isfile(name) :
        raHexen, decHexen, hexVal, rank, mjd = np.genfromtxt(name, unpack=True, 
            delimiter=",", usecols=(0,1,3,4,5))
        idHexen = np.genfromtxt(name, unpack=True, delimiter=",", usecols=(2), dtype="str")
        slots = np.ones(raHexen.size)*slot
    else :
        raHexen, decHexen, hexVal, rank, mjd, idHexen, slots = \
            7*[np.zeros(0)]

    return raHexen, decHexen, idHexen, hexVal, rank, mjd, slots

