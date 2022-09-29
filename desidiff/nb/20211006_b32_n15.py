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
import warnings

#SkyPortal token:
secret_file = "/global/cfs/cdirs/desi/science/td/secrets/desidiff_sp.txt"
with open(secret_file, 'r') as file:
    token = file.read().replace('\n', '')
headers = {'Authorization': f'token {token}'}

filter_name = 'DESIDIFF'

warnings.filterwarnings('ignore')





# read in and store in one place all the fibermap information in the spectra files
#night_arr = getUnprocessedDates()
all_plots_pdf = PdfPages("All_plots.pdf")

#for night in night_arr:   
#night_arr = [20210425, 20210428, 20210429, 20210430]
night_arr = [20211006]
#night_arr = [20210428]

a_start_time = time.time()
cumulative_count = 0
for night in night_arr:
    print(night)
    ### counters per night
    count_passed_per_night = 0
    logic_dict_night = {'narrow line':0, 'broad line':0,'Hline':0}
    for current_filename in glob(f"/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative/*/{night}/spectra-*.fits"):
        tile = current_filename.split('/')[-3]
        petal = current_filename.split('-')[1]
        ### counters per tile/petal
        tp_start_time = time.time()
        count_passed_per_tile_petal = 0
        logic_dict_tp = {'narrow line':0, 'broad line':0,'Hline':0}

        ### daily_spectra, the precursor to current spectra, before coadding, to select night, unique target ids, and individual target ids
        ### spectra.select functionality will not work once coadded with coaddition.coadd_cameras
        ### hack to deal with one known case of fibermap['NIGHT' = 20210610] while night = 20210611
        try:
            daily_spectra = ((read_spectra(current_filename)).select(nights = night))
        except:
            continue

        table = Table.read(current_filename, format='fits',hdu=1, memmap=True) 
        ##### DAVE'S ADDITION ##############
        targetcols = [i for i in table.colnames if i[-7:] =='_TARGET']
        nonzerocheck = [True in k for k in [[j != 0 for j in table[targetcols][i]] for i in range(len(table))]]
        #a really ugly line, basically generates a list of bools, 
        #True if there is at least one nonzero element in all columns ending in _TARGET
        table.remove_rows([i for i, val in enumerate(nonzerocheck) if not val])
        #This gets the index of all False values from the previous list and removes those rows
        table = table['TARGETID','TARGET_RA','TARGET_DEC','TILEID','OBJTYPE','PETAL_LOC','FIBERSTATUS','NIGHT']
        ######## END DAVE'S ADDITION ############
        table = table[numpy.logical_and(table['OBJTYPE']=='TGT', table['FIBERSTATUS']==0)]
        
        

        for ref_filename in glob(f"/global/cfs/cdirs/desi/spectro/redux/*/tiles/cumulative/{tile}/*/coadd-{petal}-{tile}*.fits"):
            if ref_filename.split('/')[-6] == 'fuji' or ref_filename.split('/')[-6] == 'guadalupe':
                prev_spectra = read_spectra(ref_filename)
                
                #### get redrock h5 file for z-shift info
                #### once a reference file is found
                try:
                    f = h5py.File((current_filename.replace('spectra','redrock')).replace('fits', 'h5'), "r")
                    print('redrock worked')
                except (FileNotFoundError):
                    f = h5py.File((current_filename.replace('spectra','rrdetails')).replace('fits', 'h5'), "r")
                    print('rrdetails worked')
                

                
                
                num = daily_spectra[numpy.where(numpy.logical_and(daily_spectra.fibermap['OBJTYPE']=='TGT', daily_spectra.fibermap['FIBERSTATUS']==0))].num_targets()
                #print('num: {}'.format(num))
                w = numpy.logical_and(daily_spectra.fibermap['OBJTYPE']=='TGT', 
                                                      daily_spectra.fibermap['FIBERSTATUS']==0)
                unique_targets = numpy.unique(daily_spectra.fibermap['TARGETID'][w])
                relevant_targets = []
                for d in unique_targets:
                    if d in prev_spectra.fibermap['TARGETID']:
                        relevant_targets.append(d)
                print('relevant targets: {}'.format(len(relevant_targets)))
                        
                        
                for t in relevant_targets:
                    
                    ### match target ids and coadd_cameras
                    current_spectra = desispec.coaddition.coadd_cameras(daily_spectra.select(targets = t))
                    ref_spectra = desispec.coaddition.coadd_cameras(prev_spectra.select(targets = t))
                    
                    ### grab ra and dec values for use in SkyPortal functionality later
                    ra, dec = table['TARGETID' == t]['TARGET_RA'], table['TARGETID' == t]['TARGET_DEC']
                    ### grab zinfo for TDE filters later
                    zinfo = f['zfit'][str(t)]['zfit'][0]['z']

                    norm = normalization(current_spectra.flux, current_spectra.mask, ref_spectra.flux, ref_spectra.mask)

                    # need to instantiate a Spectra object for the difference. 
                    ### Using 'dif_spectra = Spectra()' is bugging on dif_spectra.mask type=NoneType, can't assign.
                    #### dif_spectra = current_spectra
                    ### copy.deepcopy() is deprecated as memory expensive
                    dif_spectra = copy.deepcopy(current_spectra)
                    #### any problem with hardcoding in 'brz' for key in the following:
                    for key in (current_spectra.flux).keys():
                        current_spectra.flux[key] = current_spectra.flux[key]/norm
                        current_spectra.ivar[key] = current_spectra.ivar[key]*norm**2 
                        # subtraction of current and reference fluxes
                        dif_spectra.flux[key] = current_spectra.flux[key] - ref_spectra.flux[key]
                        ### couldn't dif_spectra.mask[key] == 2 by summing current.mask and reference.mask
                        # summation of current and reference masks
                        dif_spectra.mask[key] = current_spectra.mask[key] + ref_spectra.mask[key]
                        # inverted summation of current and spectra inverse variance
                        ### still returning RuntimeWarning: divide by zero encountered in true_divide but not in infinite loop for the moment
                        dif_spectra.ivar[key] = 1./(1./current_spectra.ivar[key] + 1./ref_spectra.ivar[key])

                    # mean-subtracted difference
                    dif_flux_clipped = clipmean(dif_spectra.flux, dif_spectra.ivar, dif_spectra.mask)

                    # filters 
                    # Difference spectrum may have broadband signal
                    #perband_filter = perband_SN(dif_flux_clipped, dif_spectra.ivar, dif_spectra.mask)
                    # fractional increase
                    #perband_inc = perband_increase(dif_flux_clipped, dif_spectra.ivar, dif_spectra.mask, ref_spectra.flux, ref_spectra.ivar, ref_spectra.mask)

                    # Difference spectrum may have high-frequency signal
                    perres_filter_broad = perconv_SN(dif_flux_clipped, dif_spectra.ivar, dif_spectra.mask, ncon = 100, nsig = 32)
                    perres_filter_narrow = perconv_SN(dif_flux_clipped, dif_spectra.ivar, dif_spectra.mask, ncon = 5, nsig = 15)

                    # Search for signature lines of TDEs, only interested in Galaxies
                    linetable = line_finder(dif_spectra.wave, dif_flux_clipped, dif_spectra.ivar, dif_spectra.mask, zinfo)
                    Hline_score = Hline_filter(linetable)
                    # deriv_score = deriv_filter(dif_flux_clipped, dif_spectra.ivar,dif_spectra.mask)
                    
                    #broadband
                    #bblogic = (numpy.array(list(perband_filter.values()))>10)#, numpy.array(list(perband_inc.values()))>0.25))
                    narrowlinelogic = perres_filter_narrow >=2
                    broadlinelogic = perres_filter_broad >=3
                    
                    # TDElogic = any([TDE_score == 2, TDE_score == 3, TDE_score == 4, TDE_score == 5])
                    Hlinelogic = any([Hline_score >= 1])
                    # derivlogic = any([deriv_score >= 3])
                    logic = [narrowlinelogic,  broadlinelogic, Hlinelogic]
                    logic_name = ['narrow line', 'broad line','Hline'] #must be in same order as logic!, use as mask
                    logic_name = numpy.ma.masked_array(logic_name, mask = [not i for i in logic])
                    
                    plt.clf()
                    plt.rcParams["figure.figsize"] = (20,6)
                    if any(logic):
                        if count_passed_per_tile_petal == 0:
                            pdf_per_file = PdfPages(str(petal)+'_'+str(tile)+".pdf")
                            
                        ### booking for filter counters
                        for e in logic_name.compressed():
                            if e in logic_dict_night.keys():
                                logic_dict_night[e] += 1
                                logic_dict_tp[e] += 1
                        # export to database of processed exposures
                        #processed(t, tile, night)
                        plt.figure()
                        for key in (current_spectra.flux).keys():
                            w=numpy.where(dif_spectra.mask[key][0] == 0)

                            plt.plot(dif_spectra.wave[key][w], dif_spectra.flux[key][0][w], color='red', label = 'Difference')
                            plt.plot(current_spectra.wave[key][w], current_spectra.flux[key][0][w], color='blue', alpha=0.5, label = 'New Spectrum')
                            plt.plot(ref_spectra.wave[key][w],ref_spectra.flux[key][0][w],color='green',alpha=0.5, label = 'Reference Spectrum')

                            plt.legend()
                            plt.xlabel('Wavelength (Å)')
                            plt.ylabel('Flux') 
                            plt.title(str(t) + "  " + str(night) + "  " + str(tile)  + "  " + str(logic_name))
                            #plt.show()

                            plt.savefig(pdf_per_file, format = 'pdf')
                            plt.savefig(all_plots_pdf, format = 'pdf')
                            plt.close()
                            #print('Time elapsed for 1 TID: {}'.format(datetime.timedelta(seconds=(time.time()-t_start_time))))
        count_passed_per_tile_petal += sum(logic_dict_tp.values())
        if count_passed_per_tile_petal !=0:
            pdf_per_file.close()
                                  
        #cumulative_count += num
        print('Time elapsed for 1 tile/petal: {}'.format(datetime.timedelta(seconds=(time.time()-tp_start_time))))
        count_passed_per_night += count_passed_per_tile_petal
        #print('TIDs passed logic per tile/petal: {}'.format(count_passed_per_tile_petal))
        print('TIDs passed per filter per tile/petal: {}'.format(logic_dict_tp))
all_plots_pdf.close()
#print('Cumulative total: {}'.format(cumulative_count))
#print('Total TIDs passed logic per night: {}'.format(count_passed_per_night))
print('Cumulative total TIDs passed per filter per night: {}'.format(logic_dict_night))

#plt.hist(list_broad, log = True)
#plt.savefig('log_broad.pdf')
#plt.hist(list_narrow, log = True)
#plt.title('log plot perres_filter narrow = orange, broad = blue, {}_{}_{}_{} tids'.format(night, petal, tile, num))
#plt.savefig('log_narrow.pdf')
print(str(datetime.timedelta(seconds=(time.time()-a_start_time))))