from __future__ import absolute_import, division, print_function
import os, sys, time

import numpy as np

import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg

from astropy.table import Table, Column

import multiprocessing


from desispec.interpolation import resample_flux
from desispec.spectra import Spectra
from desispec.resolution import Resolution
from desispec.fiberbitmasking import get_all_fiberbitmask_with_amp, get_all_nonamp_fiberbitmask_val, get_justamps_fiberbitmask
from desispec.specscore import compute_coadd_scores
from subclassSpectraNoExpid import spectra_no_expid#, coadd_no_expid, coadd_fibermap_no_expid
#from coadd_no_expid import coadd_no_expid, coadd_fibermap_no_expid

from astropy.table import Table, vstack, unique, SortedArray
import glob
import time
from datetime import date, timedelta, datetime
import psycopg2
import sqlite3
import pandas as pd
from functools import reduce
from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra
from desispec.coaddition import coadd

#from updateNoExpid import update_no_expid

import matplotlib.pyplot as plt
import os
from desiutil.log import get_logger
#from astropy import get_logger

def coadd_fibermap_no_expid(fibermap) :
    
    #log = get_logger()
    #log.debug("'coadding' fibermap")

    targets = np.unique(fibermap["TARGETID"])
    ntarget = targets.size

    jj=np.zeros(ntarget,dtype=int)
    for i,tid in enumerate(targets) :
        jj[i]=np.where(fibermap["TARGETID"]==tid)[0][0]
    tfmap=fibermap[jj]

    #- initialize NUMEXP=-1 to check that they all got filled later
    tfmap['COADD_NUMEXP'] = np.zeros(len(tfmap), dtype=np.int16) - 1
    tfmap['COADD_EXPTIME'] = np.zeros(len(tfmap), dtype=np.float32) - 1

    # smarter values for some columns
    mean_cols = [
        'DELTA_X', 'DELTA_Y',
        'FIBER_X', 'FIBER_Y',
        'FIBER_RA', 'FIBER_DEC',
        'FIBERASSIGN_X', 'FIBERASSIGN_Y'
        ]
    rms_cols = ['DELTA_X', 'DELTA_Y']  #- rms_cols must also be in mean_cols
    for k in mean_cols:
        if k in fibermap.colnames :
            if k.endswith('_RA') or k.endswith('_DEC'):
                dtype = np.float64
            else:
                dtype = np.float32
            if k in mean_cols:
                xx = Column(np.zeros(ntarget, dtype=dtype))
                tfmap.add_column(xx,name='MEAN_'+k)
            if k in rms_cols:
                xx = Column(np.zeros(ntarget, dtype=dtype))
                tfmap.add_column(xx,name='RMS_'+k)

            tfmap.remove_column(k)

    first_last_cols = ['NIGHT','EXPID','TILEID','SPECTROID','FIBER','MJD']
    for k in first_last_cols:
        if k in fibermap.colnames :
            if k in ['MJD']:
                dtype = np.float32
            else:
                dtype = np.int32
            if not 'FIRST_'+k in tfmap.dtype.names :
                xx = Column(np.arange(ntarget, dtype=dtype))
                tfmap.add_column(xx,name='FIRST_'+k)
            if not 'LAST_'+k in tfmap.dtype.names :
                xx = Column(np.arange(ntarget, dtype=dtype))
                tfmap.add_column(xx,name='LAST_'+k)
            if not 'NUM_'+k in tfmap.dtype.names :
                xx = Column(np.arange(ntarget, dtype=np.int16))
                tfmap.add_column(xx,name='NUM_'+k)

    for i,tid in enumerate(targets) :
        jj = fibermap["TARGETID"]==tid

        #- coadded FIBERSTATUS = bitwise AND of input FIBERSTATUS
        tfmap['FIBERSTATUS'][i] = np.bitwise_and.reduce(fibermap['FIBERSTATUS'][jj])

        #- Only FIBERSTATUS=0 were included in the coadd
        fiberstatus_nonamp_bits = get_all_nonamp_fiberbitmask_val()
        fiberstatus_amp_bits = get_justamps_fiberbitmask()
        targ_fibstatuses = fibermap['FIBERSTATUS'][jj]
        nonamp_fiberstatus_flagged = ( (targ_fibstatuses & fiberstatus_nonamp_bits) > 0 )
        allamps_flagged = ( (targ_fibstatuses & fiberstatus_amp_bits) == fiberstatus_amp_bits )
        good_coadds = np.bitwise_not( nonamp_fiberstatus_flagged | allamps_flagged )
        tfmap['COADD_NUMEXP'][i] = np.count_nonzero(good_coadds)
        if 'EXPTIME' in fibermap.colnames :
            tfmap['COADD_EXPTIME'][i] = np.sum(fibermap['EXPTIME'][jj][good_coadds])
        for k in mean_cols:
            if k in fibermap.colnames :
                vals=fibermap[k][jj]
                tfmap['MEAN_'+k][i] = np.mean(vals)

        for k in rms_cols:
            if k in fibermap.colnames :
                vals=fibermap[k][jj]
                # RMS includes mean offset, not same as std
                tfmap['RMS_'+k][i] = np.sqrt(np.mean(vals**2))

        for k in first_last_cols:
            if k in fibermap.colnames :
                vals=fibermap[k][jj]
                tfmap['FIRST_'+k][i] = np.min(vals)
                tfmap['LAST_'+k][i] = np.max(vals)
                tfmap['NUM_'+k][i] = np.unique(vals).size

        for k in ['FIBER_RA_IVAR', 'FIBER_DEC_IVAR','DELTA_X_IVAR', 'DELTA_Y_IVAR'] :
            if k in fibermap.colnames :
                tfmap[k][i]=np.sum(fibermap[k][jj])

    #- Remove some columns that apply to individual exp but not coadds
    #- (even coadds of the same tile)
    #for k in ['NIGHT', 'EXPID', 'MJD', 'EXPTIME', 'NUM_ITER']:
    #    if k in tfmap.colnames:
    #        tfmap.remove_column(k)

    return tfmap