from astropy.table import Table, vstack, unique, SortedArray
import glob
import numpy
import time
from datetime import date, timedelta, datetime
import pandas as pd
from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra
from desispec.coaddition import coadd
import matplotlib.pyplot as plt

import os
from desiutil.log import get_logger, DEBUG


def manualDifferencing(dats_group, yyyymmdd):
    start_time = time.time()
    d = dats_group.group_by(['NIGHT'])
    du = unique(d, keys=['NIGHT'])

    s0 = Spectra() #spectra_no_expid() ### all spectra
    s0c = Spectra() #spectra_no_expid() ### only current spectra
    s0r = Spectra() #spectra_no_expid() ### only reference spectra
    s0d = Spectra() #differenced spectra

    f0 = f"/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative/{du[0]['TILEID']}/{yyyymmdd}/spectra-{du[0]['PETAL_LOC']}-{du[0]['TILEID']}-thru{yyyymmdd}.fits"
    #this is being done multiple times: need to find way to implement select(targets = ['TID']) for files already read in
    s0 = read_spectra(f0)
    s0c = s0.select(nights = [yyyymmdd], targets = [du[0]['TARGETID']])
    for t in du:
        if t['NIGHT'] < yyyymmdd: # in future, use getPrevObsdate for days > 30 days ago
            if (s0r.num_spectra() == 0): 
            # first (or only) ref night found
                s0r = s0.select(nights = [t['NIGHT']], targets= [du[0]['TARGETID']])
                print('ref night selected ' + str(t['NIGHT']))
            else:
                s0r.update(s0.select(nights = [t['NIGHT']], targets= [du[0]['TARGETID']]))
                print('ref night updated with ' + str(t['NIGHT']))

    if s0r.num_spectra() != 0:
    # i.e. there is something to difference
        s0d.update(s0c)
        coadd(s0r)
        for i in s0d.flux:
            for j in range(len(s0d.flux[i])):
                s0d.flux[i][j] = s0c.flux[i][j] - s0r.flux[i][j]
        for i in s0d.ivar:
            ok = numpy.logical_and(s0c.ivar[i][0,:] != 0, s0r.ivar[i][0,:] != 0)
            s0d.ivar[i][0,ok] = (1/(1/s0c.ivar[i][0,ok] + 1/s0r.ivar[i][0,ok]))
        coadd(s0d)
        for b in s0d.bands:
            ok = numpy.logical_and(s0d.mask[b][0,:] == 0, s0d.ivar[b][0] != 0)
            plt.plot(s0c.wave[b][ok],s0c.flux[b][0,ok], color = 'red')
            plt.plot(s0r.wave[b][ok],s0r.flux[b][0,ok], color = 'green')
            plt.plot(s0d.wave[b][ok],s0d.flux[b][0,ok], color = 'cyan')
        #plt.xlim(3600, 3650)
        plt.show()
    else:
        print("no reference nights")
    print("--- manual differencing took:  %s seconds ---" % (time.time() - start_time))
    return(s0d)
    