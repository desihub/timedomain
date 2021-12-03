import numpy as np
from scipy.interpolate import interp1d

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
#area_left =  area_per_hex * (time_budget * 3600)/(time_cost_


#how area scale as number of event
def max_scale(N_max,alpha):
    return (2.*(1.-0.5**(1./N_max)))**(alpha/3.)

#solve for (p,N_max)    
def solve_p_N(area_left,rate,T_left, alpha, area_bar, area_bar_p):
    prange=np.arange(0,1,0.01)	#try a range of p
    N_max=np.zeros(np.size(prange))
    N_max_trial=np.arange(1.,50.,0.1)	#try a range of N_max

    # produce a useful version of the area_bar, area_bar_p curve
    area_bar=interp1d(area_bar_p,area_bar)	
    
    # solve the area_left function, 
    # find pairs of (p,N_max) that satisfy the function
    N_comp=max_scale(N_max_trial,alpha)*rate*T_left/N_max_trial	

    for i in range(0,prange.size) :
        pp = prange[i]
        ix=  np.argmin(abs(area_left-area_bar(pp)*N_comp))
        N_max[i]=N_max_trial[ix]
    
    #find the pair (p,N_max) that maximize the p_tot 
    arg=np.argmax(prange*rate*T_left/N_max)
    return prange[arg],N_max[arg]


#
# This is the average event- we evaluate the probabilty for the
# average event, so that we can consider the real event
# this is the decisison we would make for the average event
# (you would use this p on the average map to find the average area)
# 
#   inputs: 
#       area_left area that could be observed this season
#       rate      effective rate of triggers
#       T_left    time left in the season
#       cl, Ap    this is the cumulative  distribution function of area
#                   for a fixed probability p ( p_gw in our notation), 
#                   as measured on the sims
#       area_bar   this pair is the average area of the sims as a 
#       area_bar_p function of probability p (i.e. we vary p_gw, measure area_bar)
#       
def evaluate_average_event(area_left, T_left, rate, cl, Ap, area_bar, area_bar_p) :
    #
    # solve the general case
    #

    #the scaling index of the area, alpha=1.8 for 2015
    alpha = 1.8

    # solve the implicit equation for p and N_max
    p,N_max=solve_p_N( area_left, rate, T_left, alpha, area_bar, area_bar_p )	
    # N_max: how loud is the softest one 
    # p: the probability that maximizes total probability if all events were average

    return p, N_max

