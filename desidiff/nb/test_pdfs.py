import os
import sys
from glob import glob

module_path = os.path.abspath(os.path.join('../../'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
import desispec
from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

from desiutil.log import get_logger, DEBUG

from desidiff.src.group_tiles import *
from desidiff.src.dates_to_process import *
from desidiff.src.coadd import *
from desidiff.src.scores import *
from desidiff.src.ContinuumFitFilter_desidiff import *
from desidiff.src.dave_addition import *
from desidiff.src.spectra_prep import*

from timedomain.sp_utils import SkyPortal as sp
import requests
import datetime
import heapq
import time
import copy
import numpy
from astropy.time import Time
from astropy.table import Table, vstack, unique, SortedArray
import h5py

#%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from IPython import display


zee_plots_pdf = PdfPages("/project/projectdirs/desi/users/clepart/pdfs/All_plots.pdf")
for i in range(10):
    pdf_per_file = PdfPages("/project/projectdirs/desi/users/clepart/pdfs/"+str(i)+'_'+str(datetime.datetime.now())+".pdf")

pdf_per_file.close()
zee_plots_pdf.close()

