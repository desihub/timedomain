#perl_desidiff_202306.py

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






warnings.filterwarnings('ignore')

# Create a list to store the previous dates in June
#prev_night = [20230601, 20230602, 20230603, 20230604, 20230605, 20230606, 20230607, 20230608]
# Get the current date
#current_date = datetime.date.today()

prev_night = []
# #Iterate over the range of days in June
# for day in range(1, current_date.day):
#     date = datetime.date(current_date.year, 6, day)
#     prev_night.append(date.strftime('%Y%m%d'))

# #prev_night.append(int((datetime.datetime.now() - datetime.timedelta(1)).strftime('%Y%m%d')))

# #night_arr = getUnprocessedDates()

# import datetime

# Get the current date
current_date = datetime.date.today()

# Set the desired month and day
desired_month = 6  # June
desired_day = current_date.day  # Use today's day

# Create a date object with the desired month and day
desired_date = datetime.date(current_date.year, desired_month, desired_day)

# Format the desired date as 'YYYYMMDD'and append to prev_night
prev_night.append(desired_date.strftime('%Y%m%d'))


cumulative = 0

count_per_night = 0
#for night in night_arr:
for night in prev_night:
    start_time = datetime.datetime.now()
    all_plots_pdf = None
    print(night)
    logic_dict_night = {'narrow line':0, 'broad line':0,'Hline':0}
    
    count_per_current_file = 0
    for current_filename in glob(f"/global/cfs/projectdirs/desi/spectro/redux/daily/tiles/cumulative/*/{night}/spectra-*.fits*"):
        tile = current_filename.split('/')[-3]
        petal = current_filename.split('-')[1]
        count_per_current_file +=1 
        logic_dict_ref = {'narrow line':0, 'broad line':0,'Hline':0}
        
        # specify the path to the files you want to process
        ref_path = f"/global/cfs/cdirs/desi/spectro/redux/*/tiles/cumulative/{tile}/*/coadd-{petal}-{tile}*.fits*"
        # use the glob function to get a list of all the file paths matching the pattern
        
        count_per_ref_file = 0
    
        for ref_filename in glob(ref_path):
            # split the file path and check if the 6th last directory is fuji, guadalupe or iron AND not daily
            #if ref_filename.split('/')[-6] in ['fuji', 'guadalupe', 'iron'] and ref_filename.split('/')[-6] != 'daily':
            if ref_filename.split('/')[-6] == 'iron' or (ref_filename.split('/')[-6] == 'daily' and int(ref_filename.split('/')[-2]) < int(night)):
                # read the spectra from the file
                prev_spectra = read_spectra(ref_filename)
                count_per_ref_file += 1  # increment the file counter once a reference is found
                count_passed_logic = 0
            else:
                continue  # skip to the next file if the condition is not met


            ### daily_spectra, the precursor to current spectra, before coadding, to select night, unique target ids, and individual target ids
            ### spectra.select functionality will not work once coadded with coaddition.coadd_cameras
            ### hack to deal with one known case of fibermap['NIGHT' = 20210610] while night = 20210611
            try:
                daily_spectra = ((read_spectra(current_filename)).select(nights = night))
            except:
                continue
            f, table = dave_addition(current_filename)

            

            w = numpy.logical_and(daily_spectra.fibermap['OBJTYPE']=='TGT', daily_spectra.fibermap['FIBERSTATUS']==0)
            unique_targets = numpy.unique(daily_spectra.fibermap['TARGETID'][w])
            relevant_targets = []
            for d in unique_targets:
                if d in prev_spectra.fibermap['TARGETID']:
                    relevant_targets.append(d)
            """
            if len(relevant_targets) > 0:
                daily_spectra = desispec.coaddition.coadd_cameras(daily_spectra)      
                prev_spectra = desispec.coaddition.coadd_cameras(prev_spectra)
            """
            plt_start_time = datetime.datetime.now()
            for t in relevant_targets:
                ### match target ids and coadd_cameras
                coadd_start_time = datetime.datetime.now()
                current_spectra = desispec.coaddition.coadd_cameras(daily_spectra.select(targets = t))            
                ref_spectra = desispec.coaddition.coadd_cameras(prev_spectra.select(targets = t))
                #current_spectra = daily_spectra.select(targets = t)
                #ref_spectra = prev_spectra.select(targets = t)
                ra, dec, zinfo = t_dave_addition(t, table, f)
                coadd_time = datetime.datetime.now() - coadd_start_time
                # print('coadd_time ', coadd_time)
                
                clip_start_time = datetime.datetime.now()
                clip_spectra, dif_spectra, current_spectra, ref_spectra = spectra_prep(current_spectra, ref_spectra)
                clip_time = datetime.datetime.now() - clip_start_time
                # print('clip_time ', clip_time)
                
                # filters 
                # Difference spectrum may have high-frequency signal
                filters_start_time = datetime.datetime.now()
                perres_filter_broad = perconv_SN_redux(clip_spectra, dif_spectra.ivar, dif_spectra.mask, ncon = 100, nsig = 50)
                perres_filter_narrow = perconv_SN_redux(clip_spectra, dif_spectra.ivar, dif_spectra.mask, ncon = 10, nsig = 25)
                # Search for signature lines of TDEs, only interested in Galaxies
                linetable = line_finder(dif_spectra.wave, clip_spectra, dif_spectra.ivar, dif_spectra.mask, zinfo)
                Hline_score = Hline_filter(linetable)
                filters_time = datetime.datetime.now() - filters_start_time 
                # print('filters_time ', filters_time)
                
                
                #logic
                # logic_start_time = datetime.datetime.now()
                narrowlinelogic = perres_filter_narrow >=2
                broadlinelogic = perres_filter_broad >=3
                Hlinelogic = any([Hline_score >= 1])
                logic = [narrowlinelogic,  broadlinelogic, Hlinelogic]
                logic_name = ['narrow line', 'broad line', 'Hline'] #must be in same order as logic!, use as mask
                logic_name = numpy.ma.masked_array(logic_name, mask = [not i for i in logic])

                plt.clf
                #print(logic)
                # logic_time = datetime.datetime.now() - logic_start_time
                # print('logic_time', logic_time)

                if any(logic):
                    # plt_start_time = datetime.datetime.now()
                    count_passed_logic += 1
                    if all_plots_pdf is None:
                        all_plots_pdf = PdfPages("/global/cfs/projectdirs/desi/users/clepart/pdfs/All_plots_"+str(night)+".pdf")
                    ### book-keeping for filter counters
                    for e in logic_name.compressed():
                        if e in logic_dict_night.keys():
                            logic_dict_night[e] += 1
                            logic_dict_ref[e] += 1
                    ### export to database of processed exposures
                    #processed(t, tile, night)
                    plt.figure(figsize = (20,6))
                    for key in (current_spectra.flux).keys():
                        w=numpy.where(dif_spectra.mask[key][0] == 0)

                        plt.plot(dif_spectra.wave[key][w], dif_spectra.flux[key][0][w], color='red', label = 'Difference')
                        plt.plot(current_spectra.wave[key][w], current_spectra.flux[key][0][w], color='blue', alpha=0.5, label = 'New Spectrum')
                        plt.plot(ref_spectra.wave[key][w],ref_spectra.flux[key][0][w],color='green',alpha=0.5,                                  label='Reference Spectrum "' + ref_filename.split('/')[-6] + ' " ' + ref_filename.split('/')[-2])

                        plt.legend()
                        plt.xlabel('Wavelength (Å)')
                        plt.ylabel('Flux') 
                        plt.title(str(t) + "  " + str(night) + "  " + str(petal) + "  " + str(tile)  + "  " + str(logic_name))

                        plt.savefig(all_plots_pdf, format = 'pdf')
                        plt.close()
            #print('plt_time ', datetime.datetime.now() - plt_start_time)
            cumulative += len(relevant_targets)
            count_per_night += count_passed_logic
                    
            print('TIDs passed logic per tile/petal: {}'.format(count_passed_logic))
            print('TIDs passed per filter per tile/petal: {}'.format(logic_dict_ref))
            #closes plots without erroring in case of no plots created that night
        if all_plots_pdf is not None:
           all_plots_pdf.close()      
        
        
    print('Total TIDs passed logic per night: {}'.format(count_per_night))
    print('Cumulative total TIDs passed per filter per night: {}'.format(logic_dict_night))
    print('Tile/petal combinations per current filename: {}'.format(count_per_current_file))
    

    print('Cumulative total TIDs: {}'.format(cumulative))

    end_time = datetime.datetime.now()
    duration = end_time - start_time
    hours = duration.seconds // 3600
    minutes = (duration.seconds // 60) % 60
    seconds = duration.seconds % 60

    print("Duration: {} hours {} minutes {} seconds".format(hours, minutes, seconds))

