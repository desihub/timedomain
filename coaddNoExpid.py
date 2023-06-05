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
from coaddFibermapNoExpid import coadd_fibermap_no_expid

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








def coadd_no_expid(spectra_no_expid, cosmics_nsig=0.) :
    """
    Coaddition the spectra for each target and each camera. The input spectra is modified.

    Args:
       spectra: desispec.spectra.Spectra object

    Options:
       cosmics_nsig: float, nsigma clipping threshold for cosmics rays
    """

    spectra_no_expid.fibermap=coadd_fibermap_no_expid(spectra_no_expid.fibermap)
    return