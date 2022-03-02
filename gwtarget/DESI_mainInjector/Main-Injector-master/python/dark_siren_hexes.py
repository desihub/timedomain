import numpy as np
import healpy as hp
import hp2np as hp2np
import matplotlib.pyplot as plt

# make the json files
# Jan 4
# reload(dsh); dsh.get_tiling_ra_dec(raHexen, decHexen, sumProb, nu, ng, nr, ni, nz, slice=24, n=3)
# Feb 2
# reload(dsh); dsh.get_tiling_ra_dec(raHexen, decHexen, sumProb, nu, ng, nr, ni, nz, slice=30, n=4)
# reload(dsh); dsh.get_tiling_ra_dec(raHexen, decHexen, sumProb, nu, ng, nr, ni, nz, slice=62, n=3)
def get_tiling_ra_dec(raHexen, decHexen, sumProb, nu, ng, nr, ni, nz, slice=24, n=3) :

    ix=np.argsort(sumProb)[::-1]
    raHexen, decHexen, sumProb, nu, ng, nr, ni, nz = \
        raHexen[ix], decHexen[ix], sumProb[ix], nu[ix], ng[ix], nr[ix], ni[ix], nz[ix]
    raHexen, decHexen, sumProb, nu, ng, nr, ni, nz = \
        raHexen[0:slice], decHexen[0:slice], sumProb[0:slice], \
        nu[0:slice], ng[0:slice], nr[0:slice], ni[0:slice], nz[0:slice]

    feb2 = True
    if feb2:
        file = "tiling-riz-feb2-ra-dec.txt"
        file = "tiling-rz-feb2-ra-dec.txt"
        fd = open(file,"w")
        for i in range(0,raHexen.size) :
            if nz[i] < 0.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "z", 1))
            #if ni[i] < 0.5 :
            #    fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
            #        raHexen[i], decHexen[i], sumProb[i], "i", 1))
            if nr[i] < 0.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "r", 1))
            if nz[i] < 1.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "z", 2))
            #if ni[i] < 1.5 :
            #    fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
            #        raHexen[i], decHexen[i], sumProb[i], "i", 2))
            if nr[i] < 1.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "r", 2))
            if nz[i] < 2.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "z", 3))
            #if ni[i] < 2.5 :
            #    fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
            #        raHexen[i], decHexen[i], sumProb[i], "i", 3))
            if nr[i] < 2.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "r", 3))
            if nz[i] < 3.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "z", 4))
            #if ni[i] < 3.5 :
            #    fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
            #        raHexen[i], decHexen[i], sumProb[i], "i", 4))
            if nr[i] < 3.5 :
                fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                    raHexen[i], decHexen[i], sumProb[i], "r", 4))
        fd.close()
        sort_on_ra_n_write(file) 
        return
    
    file = "tiling1-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if nz[i] < 0.5 :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "z", 1))
        if ni[i] < 0.5 :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "i", 1))
        if nr[i] < 0.5 :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "r", 1))
        if ng[i] < 0.5 :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "g", 1))
    fd.close()
    sort_on_ra_n_write(file) 
    
    file = "tiling2-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if ( nz[i] < 1.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "z", 2))
        if  ( ni[i] < 1.5):
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "i", 2))
        if  ( nr[i] < 1.5):
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "r", 2))
        if ( ng[i] < 1.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "g", 2))
    fd.close()
    sort_on_ra_n_write(file) 

    file = "tiling3-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if  ( nz[i] < 2.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "z", 3))
        if  ( ni[i] < 2.5):
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "i", 3))
        if  ( nr[i] < 2.5):
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "r", 3))
        if ( ng[i] < 2.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "g", 3))
    fd.close()
    sort_on_ra_n_write(file) 

#    file = "tiling4-ra-dec.txt"
#    fd = open(file,"w")
#    for i in range(0,raHexen.size) :
#        if  ( nz[i] < 3.5) :
#            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
#                raHexen[i], decHexen[i], sumProb[i], "z",4))
#        if  ( ni[i] < 3.5):
#            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
#                raHexen[i], decHexen[i], sumProb[i], "i",4))
#        if  ( nr[i] < 3.5):
#            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
#                raHexen[i], decHexen[i], sumProb[i], "r",4))
#        if ( ng[i] < 3.5) :
#            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
#                raHexen[i], decHexen[i], sumProb[i], "g",4))
#    fd.close()#    sort_on_ra_n_write(file) 


    slice = 20

    raHexen, decHexen, sumProb, nu, ng, nr, ni, nz = \
        raHexen[0:slice], decHexen[0:slice], sumProb[0:slice], \
        nu[0:slice], ng[0:slice], nr[0:slice], ni[0:slice], nz[0:slice]

    file = "tiling1u-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if (nu[i] >= -1.5 ) and ( nu[i] < 0.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "u", 1))
    fd.close()
    sort_on_ra_n_write(file) 

    file = "tiling2u-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if  ( nu[i] < 1.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "u", 2))
    fd.close()
    sort_on_ra_n_write(file) 

    file = "tiling3u-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if  ( nu[i] < 2.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "u", 3))
    fd.close()
    sort_on_ra_n_write(file) 

    file = "tiling4u-ra-dec.txt"
    fd = open(file,"w")
    for i in range(0,raHexen.size) :
        if  ( nu[i] < 3.5) :
            fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
                raHexen[i], decHexen[i], sumProb[i], "u", 4))
    fd.close()
    sort_on_ra_n_write(file) 

def sort_on_ra_n_write(file) :
    a,b,c = np.genfromtxt(file,unpack=True,usecols=(0,1,2))
    d = np.genfromtxt(file,unpack=True,usecols=3, dtype="str")
    e = np.genfromtxt(file,unpack=True,usecols=4, dtype="int")
    ix = np.argsort(a)
    fd = open(file,"w")
    for i in range(0,a.size): 
        fd.write("{:10.5f} {:10.5f} {:10.8f} {:2s} {:2d}\n".format(
            a[ix[i]], b[ix[i]], c[ix[i]], d[ix[i]], e[ix[i]]))
    fd.close()

# for jan 4: see message on Jan 4 to slack channel #gw
# The dark siren run is a Blanco-DECam run, half night on Friday, full night on Saturday. Gray time, last two hours of the night moon free. 
# In the event, we traded the Friday night time for time around Feb 2 with Alfredo. What we're going to do tonight is as follows:
#   a) take the LALinfreance map for GW150914 (our old friend)
#   b) from the telemtry database, extract all u,g,r,i,z images with seeing < 1.5" and t_eff > 0.3. 
#   c) take the hex map from delve and us.
# Then, pixelate the GW map to healpix nside=64, eliminate low prob pixels. Take the decam image centers, filter by filter, 
# assign to the nside=64 healpix map pixel closest, and accumlate exp_time*t_eff. At the end, divide by 90s to turn into fiducial images
# Finally, take these two healpix maps and accumulate them onto the closest hex map pixel, so we obtain n_exp(filter) and total GW prob. 
# Sort on hex probability. Build target lists on the top n hexes.
# We went with n=24 which gets us close to 25% total probability with 3-image coverage in g,r,i,z, leaving about 1.5 hours at the end 
# of the night to get 15 or so u-band 4-image coverage of about half the region.
# That's the plan: the aim is to get photo-z over as much of the GW150914 area as we can. We can use the u-band data to explore the benefits of that
# We still need to delve offsets for tilings after the first- we're looking for the 3 next offset values.
#
# And we finished the complete set of tiles. 24% of spatial probability in 3 tilings in griz, 20% in 4 tilings of u. 
# We'll have to finalize the trade of time with Alfredo- Feb 2, I believe. It will be better placed to do more work on GW150914.

# for Feb 2, we'll used slice=30, n=4. This is 4.1 hours, we'll just tell alfredo to stop when its his time- ours is 3.6 hours

# reload(dsh); dsh.examine(raHexen, decHexen, sumProb, nu, ng, nr, ni, nz, slice=24, n=3)
def examine (raHexen, decHexen, sumProb, nu, ng, nr, ni, nz, slice=30, n=4) :
    ix=np.argsort(sumProb)[::-1]

    print(" i  cumProb   nu  ng  nr  ni  nz     raHexen   decHexen")
    for i in range(0,slice): 
        print("{:3d} {:5.1f}%   {:3.0f} {:3.0f} {:3.0f} {:3.0f} {:3.0f}   {:10.4f} {:10.4f}".format(
            i, 100*sumProb[ix[0:i+1]].sum(), 
            nu[ix[i]], ng[ix[i]], nr[ix[i]], ni[ix[i]],nz[ix[i]], raHexen[ix[i]], decHexen[ix[i]]))

    how (nu[ix],ng[ix],nr[ix],ni[ix],nz[ix], n=n, slice=slice) 


# used in examine
def how (nu,ng,nr,ni,nz, n=2, slice=40) :
    print ("total number of exposures")

    nexp = 0
    sum = 0
    for i in range(0,n+1) :
        ix,=np.where( (ng[0:slice]>i-0.5) & (ng[0:slice]<=i+0.5 ) ); sum = sum + (n-i)*ix.size; 
# feb 2, moon goes down at end, so no g band
    sum = 0
    print("ng", sum)
    nexp = nexp+sum

    sum = 0
    for i in range(0,n+1) :
        ix,=np.where( (nr[0:slice]>i-0.5) & (nr[0:slice]<=i+0.5 ) ); sum = sum + (n-i)*ix.size; 
    print("nr", sum)
    nexp = nexp+sum
    
    sum = 0
    for i in range(0,n+1) :
        ix,=np.where( (ni[0:slice]>i-0.5) & (ni[0:slice]<=i+0.5 ) ); sum = sum + (n-i)*ix.size; 
# feb 2, moon goes down at end, so no g band, and sheer greed- let alex do it
    sum = 0
    print("ni", sum)
    nexp = nexp+sum

    sum = 0
    for i in range(0,n+1) :
        ix,=np.where( (nz[0:slice]>i-0.5) & (nz[0:slice]<=i+0.5 ) ); sum = sum + (n-i)*ix.size; 
    print("nz", sum)
    nexp = nexp+sum
    print("nexp= ", nexp)
    print("time = ", nexp*2, "   hours=", nexp*2/60.)
# feb 2, moon goes down at end, so no u-band
    return

    sum = 0
    slice = int(slice/1.2)
    n = 4
    print("u-band slice = ",slice, " n = ",n)
    for i in range(0,n+1) :
        #if i == 1: print np.where( (nu[0:slice]>i-0.5) & (nu[0:slice]<=i+0.5 ) )
        #if i == 1: print nu[0:slice]
        ix,=np.where( (nu[0:slice]>i-0.5) & (nu[0:slice]<=i+0.5 ) ); sum = sum + (n-i)*ix.size; 
        print(i, ix.size, (n-i)*ix.size)
    print("nu", sum)
    nexp = nexp+sum
    print("nexp= ", nexp)
    print("time = ", nexp*2, "   hours=", nexp*2/60.)



# import dark_siren_hexes as dsh

# reload(dsh); raHexen, decHexen, sumProb = dsh.gw_map_hex(nside=64)

def gw_map_hex (nside=32) :
    ra, dec, vals = gw_map(nside=nside)
    raHexen,decHexen = hexen()
    # gw_hexen_ra, gw_hexen_dec = gw_rahexen_old (ra,dec,vals, raHexen, decHexen)
    raHexen, decHexen, sumProb = gw_sum_prob (ra,dec,vals, raHexen, decHexen, cut_zero=True)

    #return ra,dec,vals,gw_hexen_ra, gw_hexen_dec
    return raHexen, decHexen, sumProb

# reload(dsh); ra,dec,vals,nu,ng,nr,ni,nz = dsh.gw_map_nexp(raHexen, decHexen,  nside=64)

# get the map on raHexen, decHexen of the effective number of exposures
def gw_map_nexp (raHexen, decHexen, file="LALInference_skymap.fits.gz", nside=32) :
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/python/work4/"
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/python/work4b/"
    ra,dec,vals=hp2np.hp2np(dir+file,degrade=nside)
    if nside == 32 : ix,=np.where((vals>2e-3)&(dec<40));
    if nside == 64 : ix,=np.where((vals>5e-4)&(dec<40));
    print(ra[ix].size, vals[ix].sum())
    ra, dec, probs = ra[ix], dec[ix], vals[ix]

# the telemtry db exposures
    data = decam_map()

    # convert to fiducial 90 second exposures
    exp_times = exp_map(data, do_filter="u", nside=nside)
    exp_counts_u= exp_times[ix]/90.
    exp_times = exp_map(data, do_filter="g", nside=nside)
    exp_counts_g= exp_times[ix]/90.
    exp_times = exp_map(data, do_filter="r", nside=nside)
    exp_counts_r= exp_times[ix]/90.
    exp_times = exp_map(data, do_filter="i", nside=nside)
    exp_counts_i= exp_times[ix]/90.
    exp_times = exp_map(data, do_filter="z", nside=nside)
    exp_counts_z= exp_times[ix]/90.

    a,b, exp_counts_u = gw_sum_prob (ra,dec,exp_counts_u, raHexen, decHexen, cut_zero=False)
    a,b, exp_counts_g = gw_sum_prob (ra,dec,exp_counts_g, raHexen, decHexen, cut_zero=False)
    a,b, exp_counts_r = gw_sum_prob (ra,dec,exp_counts_r, raHexen, decHexen, cut_zero=False)
    a,b, exp_counts_i = gw_sum_prob (ra,dec,exp_counts_i, raHexen, decHexen, cut_zero=False)
    a,b, exp_counts_z = gw_sum_prob (ra,dec,exp_counts_z, raHexen, decHexen, cut_zero=False)

    return ra, dec, probs, exp_counts_u, exp_counts_g, exp_counts_r, exp_counts_i, exp_counts_z

def gw_map (file="LALInference_skymap.fits.gz", nside=32) :
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/python/work4/"
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/python/work4b/"
    ra,dec,vals=hp2np.hp2np(dir+file,degrade=nside)
    if nside == 32 : ix,=np.where((vals>2e-3)&(dec<40));
    if nside == 64 : ix,=np.where((vals>5e-4)&(dec<40));
    print(ra[ix].size, vals[ix].sum())

    return ra[ix], dec[ix], vals[ix]

# lmc cut dist= np.sqrt ( ((ra-80)*np.cos(dec*2*np.pi/360.))**2 + (dec- -72)**2 )

def decam_map (file="gw150914.csv") :
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/python/work4/"
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/python/work4b/"
    id,ra,dec,exptime,qc_fwhm,qc_teff,airmass,filter,date,propid,mag_lim_10sigma= \
        np.genfromtxt(file,unpack=True,delimiter=",")
    filter,date,propid= \
        np.genfromtxt(file,unpack=True,delimiter=",",dtype="str",usecols=(7,8,9))

# jan 4
    #ix, = np.where( (qc_fwhm <= 1.5) & ( qc_teff >= 0.3) )
# jan 29
    ix, = np.where( (qc_fwhm <= 1.5) & ( qc_teff >= 0.05) & (id>924002))
    ix, = np.where( (qc_fwhm <= 1.5) & ( qc_teff >= 0.05) )
    id,ra,dec,exptime,qc_fwhm,qc_teff,airmass,filter,date,propid,mag_lim_10sigma= \
        id[ix],ra[ix],dec[ix],exptime[ix],qc_fwhm[ix],qc_teff[ix],airmass[ix],\
        filter[ix],date[ix],propid[ix],mag_lim_10sigma[ix] 
    return id,ra,dec,exptime,qc_fwhm,qc_teff,airmass,filter,date,propid,mag_lim_10sigma


def hexen (file = "all-sky-hexCenters-decam.txt") :
    dir = "/data/des60.a/data/annis/new_mi_desgw/Main-Injector4/data/"
    raHexen,decHexen = np.genfromtxt(dir+file, unpack=True)
    return raHexen,decHexen

# return a list of hexes that contain summed "vals"
#   that could be ligo probability
#   that could be effective number of exposures
def gw_sum_prob (ra,dec,vals, raHexen, decHexen, cut_zero=True) :
    # in order of raHexen
    sumProb = np.zeros(raHexen.size)

    for i in range(0,ra.size) :
        argmin = np.argmin( ((ra[i]-raHexen)*np.cos(dec[i]*2*np.pi/360.))**2 + (dec[i]-decHexen)**2 )
        sumProb[argmin] = sumProb[argmin] + vals[i]

    if cut_zero:
        ix, = np.where(sumProb > 0)
        raHexen, decHexen, sumProb = raHexen[ix], decHexen[ix], sumProb[ix]
    return raHexen, decHexen, sumProb

# return the list of hexen at which ligo map points, in same order as ligo map
# in this routine, the return order is important (gw_hexen_ra etc)
def gw_rahexen_old (ra,dec,vals, raHexen, decHexen) :
    gw_hexen_ra = np.array([])
    gw_hexen_dec = np.array([])

    for i in range(0,ra.size) :
        argmin = np.argmin( ((ra[i]-raHexen)*np.cos(dec[i]*2*np.pi/360.))**2 + (dec[i]-decHexen)**2 )
        rh = raHexen[argmin]
        dh = decHexen[argmin]
        gw_hexen_ra = np.append(gw_hexen_ra, rh)
        gw_hexen_dec = np.append(gw_hexen_dec, dh)
    return gw_hexen_ra, gw_hexen_dec

#
# returned the summed exposure time (corrected by t_eff) for pixels on a healpy map
#   note that this returns an entire healpix map
#
def exp_map (data, do_filter="i", nside=32) :
    id,ra,dec,exptime,qc_fwhm,qc_teff,airmass,filter,date,propid,mag_lim_10sigma = data
    ix, = np.where(do_filter == filter)
    id,ra,dec,exptime,qc_fwhm,qc_teff,airmass,filter,date,propid,mag_lim_10sigma= \
        id[ix],ra[ix],dec[ix],exptime[ix],qc_fwhm[ix],qc_teff[ix],airmass[ix],\
        filter[ix],date[ix],propid[ix],mag_lim_10sigma[ix] 

    vals = np.zeros(hp.nside2npix(nside))

    pix_no = hp.ang2pix(nside, ra, dec, lonlat=True)
    unique_pix = np.unique(pix_no)
    print(do_filter, unique_pix.size)
    for i in range(0,unique_pix.size) :
        up = unique_pix[i]
        ix, = np.where(pix_no == up)

        vals[up] = (exptime[ix]*qc_teff[ix]).sum()

    return vals



import numpy as np

# Support:
#
# How to use JSONs
#   import json
#   from pprint import pprint
#   json_data=open('json_data')
#   data = json.load(json_data)
#   pprint(data)
#   json_data.close()
# 
# What is format of JSON for CTIO?
#   http://dns.ctio.noao.edu/noao/content/Preparing-Observing-Scripts
#   {"count": 3,
#       "expType": "object",
#       "object": "Sgr",
#       "filter": "g",
#       "RA": 281.0,
#       "dec": -29.866,
#       "expTime": 60
#       },

#
# This code wants a list ra,decs and will write a JSON file
# to cause the Blanco to observe them.
#
def jsonall() :
    writeJson("tiling0-ra-dec.txt",90,"tiling-0.json")
    writeJson("tiling1-ra-dec.txt",90,"tiling-1.json")
    writeJson("tiling2-ra-dec.txt",90,"tiling-2.json")
    writeJson("tiling3-ra-dec.txt",90,"tiling-3.json")
    writeJson("tiling1u-ra-dec.txt",90,"tiling-1-u.json")
    writeJson("tiling2u-ra-dec.txt",90,"tiling-2-u.json")
    writeJson("tiling3u-ra-dec.txt",90,"tiling-3-u.json")
    writeJson("tiling4u-ra-dec.txt",90,"tiling-4-u.json")

def writeJson(radecfile, exp, jsonFilename, propid="2019B-0371") :
    ra,dec,prob = np.genfromtxt(radecfile,unpack=True, usecols=(0,1,2))
    filterList = np.genfromtxt(radecfile,unpack=True,usecols=3,dtype="str")
    tilingList = np.genfromtxt(radecfile,unpack=True,usecols=4,dtype="int")
    seqtot, seqid, seqnum = 0,0,0
    nexp = 1

    offsets = tileOffsets()
    fd = open(jsonFilename,"w")
    fd.write("[\n")

    size = ra.size
    seqtot= seqtot*nexp
    print("n targets=",size)
    for i in range(0,size) :
        filter = filterList[i]
        tiling = tilingList[i]
        offsets[tiling]
        #print tiling, offsets[tiling]
        delRa = offsets[tiling][0]
        delDec = offsets[tiling][1]
        tra = ra[i]
        tdec = dec[i]
        tdec = tdec+delDec
        tra = tra + delRa/np.cos(tdec*2*np.pi/360.)
        if tra < 0 : tra = tra+360.
        if tra > 360. : tra = tra-360.
        comment = "dark siren "
        object = comment

        fd.write("{")
        fd.write(" \"expType\" : \"object\",\n")
        fd.write("  \"object\" : \"{}\",\n".format(object))
        fd.write("  \"seqid\" : \"{}\",\n".format(seqid))
        fd.write("  \"seqnum\" : \"{:d}\",\n".format(int(seqnum)))
        fd.write("  \"seqtot\" : \"{:d}\",\n".format(int(seqtot)))
        fd.write("  \"expTime\" : {:d},\n".format(int(exp)))
        fd.write("  \"wait\" : \"False\",\n")
        fd.write("  \"count\" : \"1\",\n")
        fd.write("  \"note\" : \"Added to queue from desgw json file, not obstac\",\n")
        fd.write("  \"filter\" : \"{}\",\n".format(filter))
        fd.write("  \"program\" : \"des gw dark siren\",\n")
        fd.write("  \"RA\" : {:.6f},\n".format(tra))
        fd.write("  \"dec\" : {:.5f},\n".format(tdec))
        fd.write("  \"propid\" : \"{}\",\n".format(propid))
        fd.write("  \"comment\" : \"{}\"\n".format(comment)) 
        # note lack of comma for end
        fd.write("}")
        if (i == size-1) :
                pass
        else :
            fd.write(",")
        fd.write("\n")
    fd.write("]\n")
    fd.close()

# production offsets are in DES docdb 7269
# production offsets are one based, not zero based:
def tileOffsets() :
    offsets = dict()
# blessed offsets 0
    offsets[0] = 0.000, 0.000
# we'll use productions 9,10 (as if we were doing DES year 5)
# production offsets: docdb 7269, offsets_cross.txt
    offsets[1] = 0.000, 0.000
    offsets[2] = -0.76668, 0.473424 
    offsets[3] = -0.543065, -0.828492
    offsets[4] = 0.0479175, 0.777768 
    offsets[5] = 0.06389, 0.287436 
    offsets[6] = -0.4632025, 0.490332 
    offsets[7] = 0.9423775, 0.405792 
    offsets[8] =-0.2395875, -0.135264 
# production offsets 9
    offsets[9] = 0.76668, 0.4227
# production offsets 10
    offsets[10] = -0.0479175, 0.388884
# production offsets 10
    offsets[11] = -0.5257, 0.7222
# production offsets 17,18
    offsets[17] = -1.1388, 0.0166
    offsets[18] =  0.0484, -0.6725
# fake offset 19
    offsets[19] =  0.9423775, -0.405792 
#
    return offsets
