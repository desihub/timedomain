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
start_time = datetime.datetime.now()
prev_night = []
#for i in range(4):
#    prev_night.append(int((datetime.datetime.now() - datetime.timedelta(i)).strftime('%Y%m%d')))

prev_night.append(int((datetime.datetime.now() - datetime.timedelta(1)).strftime('%Y%m%d')))

#night_arr = getUnprocessedDates()



cumulative = 0
count_

#for night in night_arr:
for night in prev_night:
    print(night)
    count_night = 0
    logic_dict_night = {'narrow line':0, 'broad line':0,'Hline':0}
    count_filename = 0
    for current_filename in glob(f"/global/cfs/projectdirs/desi/spectro/redux/daily/tiles/cumulative/*/{night}/spectra-*.fits*"):
        tile = current_filename.split('/')[-3]
        petal = current_filename.split('-')[1]
        logic_dict_tp = {'narrow line':0, 'broad line':0,'Hline':0}
        count_tp = 0
   
        # specify the path to the files you want to process
        ref_path = f"/global/cfs/cdirs/desi/spectro/redux/*/tiles/cumulative/{tile}/*/coadd-{petal}-{tile}*.fits*"
        # use the glob function to get a list of all the file paths matching the pattern
        for ref_filename in glob(ref_path):
            # split the file path and check if the 6th last directory is fuji, guadalupe or iron AND not daily
            #if ref_filename.split('/')[-6] in ['fuji', 'guadalupe', 'iron'] and ref_filename.split('/')[-6] != 'daily':
            if ref_filename.split('/')[-6] == 'iron' or (ref_filename.split('/')[-6] == 'daily' and int(ref_filename.split('/')[-2]) < night):
                # read the spectra from the file
                prev_spectra = read_spectra(ref_filename)
                count_filename += 1  # increment the file counter once a reference is found
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
            for t in relevant_targets:
                ### match target ids and coadd_cameras
                current_spectra = desispec.coaddition.coadd_cameras(daily_spectra.select(targets = t))            
                ref_spectra = desispec.coaddition.coadd_cameras(prev_spectra.select(targets = t))

                ra, dec, zinfo = t_dave_addition(t, table, f)

                clip_spectra, dif_spectra, current_spectra, ref_spectra = spectra_prep(current_spectra, ref_spectra)

                # filters 
                # Difference spectrum may have high-frequency signal
                perres_filter_broad = perconv_SN_redux(clip_spectra, dif_spectra.ivar, dif_spectra.mask, ncon = 100, nsig = 32)
                perres_filter_narrow = perconv_SN_redux(clip_spectra, dif_spectra.ivar, dif_spectra.mask, ncon = 5, nsig = 15)
                # Search for signature lines of TDEs, only interested in Galaxies
                linetable = line_finder(dif_spectra.wave, clip_spectra, dif_spectra.ivar, dif_spectra.mask, zinfo)
                Hline_score = Hline_filter(linetable)
                #logic
                narrowlinelogic = perres_filter_narrow >=2
                broadlinelogic = perres_filter_broad >=3
                Hlinelogic = any([Hline_score >= 1])
                logic = [narrowlinelogic,  broadlinelogic, Hlinelogic]
                logic_name = ['narrow line', 'broad line','Hline'] #must be in same order as logic!, use as mask
                logic_name = numpy.ma.masked_array(logic_name, mask = [not i for i in logic])

                plt.clf
                if any(logic):
                    count_tp += 1
                    if count_tp == 1:
                        all_plots_pdf = PdfPages("/global/cfs//projectdirs/desi/users/clepart/pdfs/All_plots_"+str(night)+".pdf")
                    ### book-keeping for filter counters
                    for e in logic_name.compressed():
                        if e in logic_dict_night.keys():
                            logic_dict_night[e] += 1
                            logic_dict_tp[e] += 1
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
                        
            cumul_count += len(relevant_targets)
            count_night += count_tp
            print('TIDs passed logic per tile/petal: {}'.format(count_tp))
            print('TIDs passed per filter per tile/petal: {}'.format(logic_dict_tp))
    print('Total TIDs passed logic per night: {}'.format(count_night))
    print('Cumulative total TIDs passed per filter per night: {}'.format(logic_dict_night))
    print('Tile/petal combinations per current filename: {}'.format(count_filename))
    #closes plots without erroring in case of no plots created that night
    if count_tp >= 1:
        all_plots_pdf.close()
print('Cumulative total TIDs: {}'.format(cumul_count))

end_time = datetime.datetime.now()
duration = end_time - start_time
hours = duration.seconds // 3600
minutes = (duration.seconds // 60) % 60

print("Duration: {} hours {} minutes".format(hours, minutes))

