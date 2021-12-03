import numpy as np

def cappedMolly() :
    import rotate
    mra,mdec = molly(do_not_too_close=True)
    sra,sdec = pureSine()
    ix= (mdec>=-65.0) & (mdec < 65.0)
    c_mra = mra[ix]; c_mdec=mdec[ix]
    ix2 = ((sra*np.cos(sdec*2*np.pi/360.))**2+sdec**2) <= 25.5**2
    n_era,n_edec = rotate.rotateAboutEasternPole(0,90,sra[ix2],sdec[ix2]);
    s_era,s_edec = rotate.rotateAboutEasternPole(0,-90,sra[ix2],sdec[ix2]);
    b_ix = n_era > 180; n_era[b_ix]=n_era[b_ix]-360.
    b_ix = s_era > 180; s_era[b_ix]=s_era[b_ix]-360.
    
    n_era,n_edec = not_too_close(n_era, n_edec, c_mra, c_mdec)
    s_era,s_edec = not_too_close(s_era, s_edec, c_mra, c_mdec)

    ra = np.concatenate([n_era, c_mra, s_era])
    dec = np.concatenate([n_edec, c_mdec, s_edec])
    return ra,dec

def not_too_close(ra,dec,master_ra,master_dec) :
    nra, ndec=[],[]
    for i in range(0,ra.size) :
        ix = np.sqrt( ((ra[i]-master_ra)*np.cos(dec[i]*2*np.pi/360.))**2 + 
                      (dec[i]-master_dec)**2 ) >=  1.0
        if np.all(ix) :
            nra.append(ra[i])
            ndec.append(dec[i])
    nra=np.array(nra)
    ndec=np.array(ndec)
    return nra, ndec
    
#
# ra,dec = np.genfromtxt("hexCenters-v20-r13.txt", unpack=True, skiprows=3,usecols=(0,1))
# mra,mdec=np.genfromtxt("hexCenters_sinu-molly14_more82.txt", unpack=True,usecols=(0,1))
#
def molly ( do_not_too_close=False) :
    hexenX, hexenY = euclideanHexen()
    ra, dec = np.array([]),np.array([])

    zone_zero = 6

    seam_list = [[-180,-150], [-150.5,-120], [-120.5,-90], 
                [-90.5,-59.8], [-60.25,-34.7],  [-35.7,-3.0], 
                [-3.6,29.0], [28.3,61.0], [60.0,93.0], 
                [91.9,120.5], [120,150.5], [150,180]]
    zone_mid = []
    for i in range(0,len(seam_list)) :
        zone = seam_list[i]
        zone_a = zone[0]
        zone_b = zone[1]
        zone_mid.append( zone_a+(zone_b-zone_a)/2. )
    zone_mid[4] = -51.1
    zone_mid[5] = -19.15
    zone_mid[zone_zero] = 12.8
    zone_mid[7] = 44.7
    zone_mid[8] = 76.7
    zone_mid[9] = 108.63

        
    zones = range(6,len(seam_list)) + range(5,-1,-1)
    for i in zones :
        zone = seam_list[i]
        zone_a = zone[0]
        zone_b = zone[1]

        zone_ra, zone_dec = invSinProjection (hexenX, hexenY, lon_zero=zone_mid[i] ) 
        ix = (zone_ra > zone_a) & (zone_ra <= zone_b)
        zra  = zone_ra[ix]
        zdec = zone_dec[ix]
        if do_not_too_close :
            zra,zdec = not_too_close(zra, zdec, ra, dec)
        ra = np.concatenate([ra, zra])
        dec = np.concatenate([dec, zdec])
    return ra, dec

def pureSine () :
    hexenX, hexenY = euclideanHexen()
    ra, dec = invSinProjection (hexenX, hexenY)
    return ra, dec


# ix=(mdec<-50);plt.clf();plt.scatter(mra[ix],mdec[ix]);ax = figure.add_subplot(111)
# for i in range(mra[ix].size) :
#   hpath=decam2hp.hexPath(mra[ix][i],mdec[ix][i]);
#   patch = matplotlib.patches.PathPatch(hpath);ax.add_patch(patch)
# plt.show()

#
def euclideanHexen () :
    xgap     = 2.450; 
    ccdScale = 2048*15/1000.
    camHeight = (12*ccdScale + 11*xgap)*0.27*(1/15.)*1000./60./60. ;# degrees
    interrow_spacing = camHeight
    interrow_spacing = interrow_spacing + 0.05
    print "interrow_spacing= ", interrow_spacing

    ygap = 3.096
    ccdScale = 4096*15/1000.
    offset = (5*ccdScale + 4*ygap)*0.27*(1/15.)*1000./60./60. ;# degrees
    intercol_spacing = offset
    intercol_spacing = intercol_spacing/1.0001
    print "intercol_spacing= ", intercol_spacing


    # this is going to be a flat sheet of hexes
    ra_zero = -180-0.1443
    dec_zero = -90-0.6437
    ra_zero = ra_zero+0
    hexenRa,  hexenDec = [], []
    for j in range(0,181,1) :
        for i in range(0,361,1) :
            pseudo_ra   = ra_zero + i*intercol_spacing
            if np.mod(i, 2) == 0 :
                pseudo_dec = dec_zero + j*interrow_spacing
            else :
                pseudo_dec = dec_zero + (j+0.5)*interrow_spacing
            if pseudo_ra  > 179 : continue
            if pseudo_dec >  89 : continue
            if pseudo_ra  < -179 : continue
            if pseudo_dec <  -89 : continue
            hexenRa.append(pseudo_ra)
            hexenDec.append(pseudo_dec)
    hexenRa = np.array(hexenRa)
    hexenDec = np.array(hexenDec)
    return hexenRa, hexenDec

def sinusoidalProjection ( lon, lat, lon_zero=0 ) :
    y = lat
    x = (lon-lon_zero)*np.cos(lat*2*np.pi/360.)
    return x,y

# x,y in degrees
def invSinProjection (x, y, lon_zero=0 ) :
    lat = y
    lon = x/np.cos(lat*2*np.pi/360.) + lon_zero
    ix = (lat > -90) & (lat < 90) & (lon > -180) & (lon < 180)
    return lon[ix], lat[ix]

