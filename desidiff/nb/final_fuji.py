import glob
import sys
import os
import requests
from desispec.io import read_spectra
from astropy.table import Table
import numpy
import matplotlib.pyplot as plt
from desiutil.log import get_logger, DEBUG
from desidiff.src.group_tiles import *
from desidiff.src.dates_to_process import *
from desidiff.src.coadd import *
from desidiff.src.scores import *
from desidiff.src.ContinuumFitFilter_desidiff import *
import desispec

import os
import sys
module_path = os.path.abspath(os.path.join('../../'))
if module_path not in sys.path:
    sys.path.append(module_path)
import numpy
from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

import matplotlib.pyplot as plt
from IPython import display
from astropy.table import Table,vstack
from desiutil.log import get_logger, DEBUG
from desidiff.src.group_tiles import *
from desidiff.src.dates_to_process import *
from desidiff.src.coadd import *
from desidiff.src.scores import *
from desidiff.src.ContinuumFitFilter_desidiff import *
import requests
from timedomain.sp_utils import SkyPortal as sp
import datetime
from astropy.time import Time

from matplotlib.backends.backend_pdf import PdfPages

import heapq

import time

start_time = time.time()

#Set non-default plot size 
plt.rcParams["figure.figsize"] = (20,6)

path = "/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative"


#SkyPortal token:
secret_file = "/global/cfs/cdirs/desi/science/td/secrets/desidiff_sp.txt"
with open(secret_file, 'r') as file:
    token = file.read().replace('\n', '')
headers = {'Authorization': f'token {token}'}

filter_name = 'DESIDIFF'

pdf1 = PdfPages("All_plots.pdf")

lminb=3700.
lminr=5800.
lmaxr=7580.
lmaxz=9100.

count_tile = 0
count_candidates = 0

passed = False

#wrong tiles: 80613, 80693, 80605, 80690, 80650, 80610

for tile in (next(os.walk(path))[1]):
    if (str(tile) != '80610' and (not passed)):
        continue
    if (str(tile) == '80610'):
        passed = True
        continue
    pdf = PdfPages(str(tile)+".pdf")
    path2 = (os.path.join(path, tile))
    for date in (next(os.walk(path2))[1]):
        path3 = (os.path.join(path2, date))
        for file in (next(os.walk(path3))[2]):
            if (file[:7] == "spectra"):
                filename = (os.path.join(path3, file))
                ra_dec = Table.read(filename, format='fits',hdu=1, memmap=True)['TARGETID','TARGET_RA', 'TARGET_DEC']
                #Read in the two files that contain the spectra and redshifts
                spectra = read_spectra(filename)
                rr = Table.read(filename.replace('spectra','redrock'), format='fits',hdu=1, memmap=True)['TARGETID','Z','ZERR','ZWARN','SPECTYPE']
                #If there is only one night then exit
                if len(numpy.unique(spectra.fibermap['NIGHT'])) == 1:
                    continue          

                #Determine the unique targetid's
                unique_targetid = numpy.unique(spectra.fibermap['TARGETID'][numpy.logical_and(spectra.fibermap['OBJTYPE']=='TGT', spectra.fibermap['FIBERSTATUS']==0)])

                #Loop over the unique targetids
                for t in unique_targetid:
                    # a spectra object with the subset of the targetid
                    target_spectra = spectra.select(targets=[t])

                    # the number of unique nights for this target id
                    unique_nights = numpy.unique(target_spectra.fibermap['NIGHT'])

                    # if there is only one night there is no subtraction possible
                    if len(unique_nights) <=1:
                        print("nothing to do here")
                        break

                    count_candidates += 1
                    # the redshift information for this object
                    zinfo = rr[rr['TARGETID']==t]
                    
                    ra_dec_data = ra_dec[ra_dec['TARGETID'] == t]
                    # Do a search for each unique night.  Note that the filters are parity violating.
                    for night in unique_nights:
                        idx = numpy.in1d(unique_nights, night)
                        ref_nights = unique_nights[~idx]
                        ## build reference
                        refSpectra = (spectra.select(nights=ref_nights, targets = t))  

                        ## build new
                        newSpectra = (spectra.select(nights=night, targets = t))

                        ## search

                        newSpectrum = desispec.coaddition.coadd_cameras(newSpectra)
                        refSpectrum = desispec.coaddition.coadd_cameras(refSpectra)

                        newflux = newSpectrum.flux
                        newivar = newSpectrum.ivar
                        newwave = newSpectrum.wave
                        newmask = newSpectrum.mask

                        refflux = refSpectrum.flux
                        refivar = refSpectrum.ivar
                        refwave = refSpectrum.wave
                        refmask = refSpectrum.mask

                        newflux['brz'] = newflux['brz'][0]
                        newivar['brz'] = newivar['brz'][0]
                        newwave['brz'] = newwave['brz']
                        newmask['brz'] = newmask['brz'][0]

                        refflux['brz'] = refflux['brz'][0]
                        refivar['brz'] = refivar['brz'][0]
                        refwave['brz'] = refwave['brz']
                        refmask['brz'] = refmask['brz'][0]


                        # renormalize spectra to match each other
                        # There is a significant background of spectra that have the same shape but different fluxes
                        # This seems to be related to mistaken coordinates of bright sources
                        norm = normalization(newflux, newmask, refflux,refmask)

                        for key in newflux.keys():
                            newflux[key]=newflux[key]/norm
                            newivar[key]=newivar[key]*norm**2        

                        #difflux, difivar, difmask, difwave = dict.fromkeys(["b", "r", "z"]), dict.fromkeys(["b", "r", "z"]), dict.fromkeys(["b", "r", "z"]), dict.fromkeys(["b", "r", "z"])

                        difflux = {key: newflux[key] - refflux[key]
                                       for key in newflux.keys()}
                        difivar = {key: 1./(1./newivar[key] + 1./refivar[key])
                                       for key in newivar.keys()}
                        difmask = {key: newmask[key] + refmask[key]
                                       for key in newmask.keys()}
                        difwave = dict(newwave)

                        # Mask systematically problematic areas of the spectrum
                        # trim red edge 
                        if 'z' in difflux.keys():
                            difmask['z'][difwave['z'] > lmaxz]=1

                        # trim blue edge 
                        if 'b' in difflux.keys():
                            difmask['b'][difwave['b'] < lminb]=1

                        if 'r' in difflux.keys():
                            difmask['r'][difwave['r'] < lminr]=1
                            difmask['r'][difwave['r'] > lmaxr]=1

                        # mean-subtracted difference
                        difflux_clipmean = clipmean(difflux,difivar,difmask)

                        ## filters 

                        #Difference spectrum may have broadband signal
                        perband_filter = perband_SN(difflux_clipmean,difivar,difmask)

                        #fractional increase
                        perband_inc = perband_increase(difflux_clipmean,difivar,difmask, refflux,refivar,refmask)

                        # Difference spectrum may have high-frequency signal
                        perres_filter_broad = perconv_SN(difwave, difflux_clipmean,difivar,difmask, ncon = 100, nsig = 17)
        
                        perres_filter_narrow = perconv_SN(difwave, difflux_clipmean,difivar,difmask, ncon = 5, nsig = 10)

                        # Search for signature lines of TDEs, only interested in Galaxies
                        linetable = line_finder(difwave, difflux_clipmean,difivar,difmask,zinfo['Z'][0])
                        # print(night, linetable)
                        spectype = zinfo['SPECTYPE'][0]
                        Hline_score = Hline_filter(linetable)
                        # deriv_score = deriv_filter(difflux_clipmean, difivar,difmask)

                        #broadband
                        bblogic = any(numpy.logical_and(numpy.array(list(perband_filter.values()))>10, numpy.array(list(perband_inc.values()))>0.25))
                        narrowlinelogic = perres_filter_narrow >=2
                        broadlinelogic = perres_filter_broad >=3
                        # TDElogic = any([TDE_score == 2, TDE_score == 3, TDE_score == 4, TDE_score == 5])
                        Hlinelogic = any([Hline_score >= 1])
                        # derivlogic = any([deriv_score >= 3])

                        logic = [bblogic, narrowlinelogic,  broadlinelogic, Hlinelogic]
                        logic_name = ['Broadband', 'narrow line', 'broad line','Hline'] #must be in same order as logic!, use as mask
                        logic_name = numpy.ma.masked_array(logic_name, mask = [not i for i in logic])
                        plt.clf()
                        if any(logic):
                            #Uncomment next line if you want to print only those TargetIds that get plotted
                            processed(t, tile, night)
                            fig = plt.figure()
                            for b in difflux.keys():
                                w=numpy.where(difmask[b] ==0)[0]
                                if b == list(difflux.keys())[-1]:
                                    plt.plot(difwave[b][w],difflux[b][w],color='red', label = 'Difference')
                                    plt.plot(newwave[b][w],newflux[b][w],color='blue',alpha=0.5, label = 'New Spectrum')
                                    plt.plot(refwave[b][w],refflux[b][w],color='green',alpha=0.5, label = 'Reference Spectrum')
                                else:
                                    plt.plot(difwave[b][w],difflux[b][w],color='red')
                                    plt.plot(newwave[b][w],newflux[b][w],color='blue',alpha=0.5)
                                    plt.plot(refwave[b][w],refflux[b][w],color='green',alpha=0.5)
                            plt.legend()
                            plt.xlim((lminb,lmaxz))
                            plt.xlabel('Wavelength (A)')
                            plt.ylabel('Flux') 
                            plt.title(str(t) + "  " + str(night) + "  " + str(tile)  + "  " + str(logic_name))
                            plt.savefig(pdf, format = 'pdf')
                            plt.savefig(pdf1, format = 'pdf')
                            plt.close()
                            
                            #SkyPortal functionality begins:

                            #SkyPortal's Id for the object
                            objID = 'DESI{}'.format(str(t))

                            #Code to check if this object already exists in SkyPortal. If not, create it

                            response = requests.get("https://desi-skyportal.lbl.gov/api/candidates/{}".format(objID), headers={"Authorization": f"token {token}"})
                            if response.status_code == 400:
                                obj_data = {
                                        "ra": ra_dec_data['TARGET_RA'][0], #RA is required when creating a new object.
                                        "dec": ra_dec_data['TARGET_DEC'][0], #Same for DEC
                                        "id": objID,
                                        "redshift": zinfo['Z'][0],
                                        "filter_ids": [sp.filter_id(filter_name)],
                                        "passed_at": str(datetime.datetime.utcnow()) #UTC time when the object passed the filter
                                        }

                                response = requests.post(
                                        "https://desi-skyportal.lbl.gov/api/candidates",
                                        json= obj_data,
                                        headers={"Authorization": f"token {token}"})

                            #Now, we send the differenced spectrum to SkyPortal
                            #First, check if the same spectrum already exists under the object
                            response = requests.get("https://desi-skyportal.lbl.gov/api/sources/{}/spectra".format(objID),headers={"Authorization": f"token {token}"})
                            post = True
                            # for i in range (len(response.json()['data']['spectra'])):
                            #     if (wavelengths == response.json()['data']['spectra'][i]['wavelengths'] and flux == response.json()['data']['spectra'][i]['fluxes']):
                            #         post = False
                            #Only if the exact same spectrum doesn't already exist, upload it.
                            if post:
                                #Send difference spectrum to SP:
                                diff_spectrum_data = {
                                        "obj_id": objID,
                                        "wavelengths": difwave['brz'].tolist(),
                                        "fluxes": difflux['brz'].tolist(),
                                        "observed_at": str(night)[:4]+'-'+str(night)[4:6]+'-'+str(night)[6:]+' '+'00:00:00.000000', # Date converted into UTC time format
                                        "origin": "DESIDIFF_everest", #Only the difference spectrum gets this tag in order to distinguish it on SkyPortal
                                        "group_ids": [sp.group_id("DESI")],
                                        "instrument_id": sp.instrument_id()
                                        }
                                response = requests.post(
                                            '{}/api/spectrum'.format(sp.url),
                                            json= diff_spectrum_data,
                                            headers={"Authorization": f"token {token}"}) 

                                #Send new spectra to SP:
                                new_spectra_data = {
                                    "obj_id": objID,
                                    "wavelengths":newwave['brz'].tolist(),
                                    "fluxes": newflux['brz'].tolist(),
                                    "observed_at": str(night)[:4]+'-'+str(night)[4:6]+'-'+str(night)[6:]+' '+'00:00:00.000000', # Date converted into UTC time format
                                    "origin": "DESIDIFF_everest",
                                    "group_ids": [sp.group_id("DESI")],
                                    "instrument_id": sp.instrument_id()
                                    }
                                response = requests.post(
                                        '{}/api/spectrum'.format(sp.url),
                                        json= new_spectra_data,
                                        headers={"Authorization": f"token {token}"}) 

                                #Code to find the ref_night that is closest to new_night to be used as 'observed_at' for SP
                                new_night = datetime.datetime(int(str(night)[:4]), int(str(night)[4:6]), int(str(night)[6:]))
                                closest_night = datetime.datetime(int(str(ref_nights[0])[:4]), int(str(ref_nights[0])[4:6]), int(str(ref_nights[0])[6:]))
                                for ref_night in ref_nights[1:]:
                                    d = datetime.datetime(int(str(ref_night)[:4]), int(str(ref_night)[4:6]), int(str(ref_night)[6:]))
                                    time_delta = new_night - d
                                    best = new_night - closest_night 
                                    if (time_delta.total_seconds() < best.total_seconds()):
                                        closest_night = d
                                    

                                ref_spectra_data = {
                                    "obj_id": objID,
                                    "wavelengths": refwave['brz'].tolist(),
                                    "fluxes": refflux['brz'].tolist(),
                                    "observed_at": str(closest_night),
                                    "origin": "DESIDIFF_everest",
                                    "group_ids": [sp.group_id("DESI")],
                                    "instrument_id": sp.instrument_id()
                                    }
                                response = requests.post(
                                        '{}/api/spectrum'.format(sp.url),
                                        json= ref_spectra_data,
                                        headers={"Authorization": f"token {token}"})    

    pdf.close()
pdf1.close()
