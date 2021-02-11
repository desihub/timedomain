#!/usr/bin/env python
# coding: utf-8

# # Apply CNN Classifier to DESI Spectra
# 
# ## SV0 and mini-SV2
# 
# Mini-SV2 tiles from February-March 2020:
# - https://desi.lbl.gov/trac/wiki/TargetSelectionWG/miniSV2
# 
# See also the DESI tile picker with (limited) SV0 tiles from March 2020:
# - https://desi.lbl.gov/svn/data/tiles/trunk/
# - https://desi.lbl.gov/svn/data/tiles/trunk/SV0.html
# 
# ## SV1
# 
# SV1 tiles from December-January 2020 are being updated daily on the wiki:
# https://desi.lbl.gov/trac/wiki/SurveyValidation/SV1
# 
# See also the DESI e-log for a less curated list:
# http://desi-www.kpno.noao.edu:8090/ECL/desi/E/search?text=sv1&search=Search

import os
import sys
sys.path.append('/global/homes/d/divij18/timedomain/desitrip/py/') #Note:change this path as needed!

from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra
from desispec.coaddition import coadd_cameras
from desitarget.sv1.sv1_targetmask import bgs_mask

from desitrip.preproc import rebin_flux, rescale_flux

from astropy.io import fits
from astropy.table import Table, vstack, hstack, join
import fitsio

from glob import glob
from datetime import date

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from tensorflow import keras

mpl.rc('font', size=14)

from argparse import ArgumentParser

import sqlite3

#DESITRIP_daily
from timedomain.sp_utils import *
from timedomain.filters import *
from timedomain.iterators import *

#sqlite commands
db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'
conn = sqlite3.connect(db_filename)
c = conn.cursor()


if __name__ == '__main__':
    parser = ArgumentParser(description='Run CNN classifier on selected tile and day')

    parser.add_argument('-d', '--obsdate', type=str,default='20201221',
                        help='Night to be processed')
    parser.add_argument('-t', '--tilenum', nargs='+', type=str,default='80611',
                        help='Tile Numbers to be processed')
    parser.add_argument('-r', '--redux', default='daily',
                        help='Spectroscopic reduction: daily, andes, blanc, ...')

    args = parser.parse_args()
    
    #DESITRIP_daily
#     logic = getattr(sys.modules[__name__], args.logic)

    tile_numbers = args.tilenum
    obsdate = args.obsdate

    base_path='/global/u2/p/palmese/desi/timedomain/cronjobs/'
    td_path='/global/cfs/cdirs/desi/science/td/daily-search/desitrip/'
    plot_path=td_path+'plots/'
    out_path=td_path+'out/'
    # Set up BGS target bit selection.
    sv1_bgs_bits = '|'.join([_ for _ in bgs_mask.names() if 'BGS' in _])


    # ## Load the Keras Model
    # 
    # Load a model trained on real or simulated data using the native Keras output format. In the future this could be updated to just load the Keras weights.


    tfmodel = '/global/homes/l/lehsani/timedomain/desitrip/docs/nb/models_9label_first/6_b65_e200_9label/b65_e200_9label_model'
    # tfmodel = '/'.join([os.environ['HOME'], 'desi/timedomain/desitrip/docs/nb', '6label_cnn_restframe'])
    if os.path.exists(tfmodel):
        classifier = keras.models.load_model(tfmodel)
    else:
        classifier = None
        print('Sorry, could not find {}'.format(tfmodel))

    if classifier is not None:
        classifier.summary()

    label_names = ['Galaxy',
              'SN Ia',
              'SN Ib',
              'SN Ib/c',
              'SN Ic',
              'SN IIn',
              'SN IIL/P',
              'SN IIP',
              'KN']
    label_names_arr=np.array(label_names)

    # ## Loop Through Spectra and Classify

    def get_petal_id(filename):
        return int(filename.split('-')[1])

    def match_files(cafiles, zbfiles):  
        matched = [] 
        for _l1 in cafiles: 
            p1 = get_petal_id(_l1)  
            for _l2 in zbfiles: 
                p2 = get_petal_id(_l2) 
                if p1 == p2: 
                    matched.append([_l1, _l2]) 
                    break
        return matched 

    #Prepare variables to save for transients (tr)
    tr_targetid=None
    tr_z=None
    tr_label=None
    tr_tile=None
    tr_date=None
    tr_spectrum=None

    # Access redux folder.
    redux = '/'.join([os.environ['DESI_SPECTRO_REDUX'], args.redux, 'tiles'])
    prefix_in = '/'.join([redux, tile_number, obsdate])
    if not os.path.isdir(prefix_in):
        print('{} does not exist.'.format(prefix_in))
    else:
        cframefiles = sorted(glob('{}/cframe*.fits'.format(prefix_in)))
        cframe_fibermap = fitsio.read_header(cframefiles[0],'FIBERMAP')
        if ('BGS' in cframe_fibermap['PROGRAM']):
            # Access the zbest and coadd files.
            # Files are organized by petal number.
            zbfiles = sorted(glob('{}/zbest*.fits'.format(prefix_in)))
            cafiles = sorted(glob('{}/coadd*.fits'.format(prefix_in)))

            # Loop through zbest and coadd files for each petal.
            # Extract the fibermaps, ZBEST tables, and spectra.
            # Keep only BGS targets passing basic event selection.
            allzbest = None
            allfmap = None
            allwave = None
            allflux = None
            allivar = None
            allmask = None
            allres  = None

            print('Tile: {} - {}'.format(tile_number, obsdate))

            for cafile, zbfile in match_files(cafiles, zbfiles):
                print('  - Petal {}'.format(get_petal_id(cafile)))

                # Access data per petal.
                zbest = Table.read(zbfile, 'ZBEST')
                pspectra = read_spectra(cafile)
                cspectra = coadd_cameras(pspectra)
                fibermap = cspectra.fibermap

                # Apply standard event selection.
                isTGT = fibermap['OBJTYPE'] == 'TGT'
                isGAL = zbest['SPECTYPE'] == 'GALAXY'
                isBGS = fibermap['SV1_BGS_TARGET'] & bgs_mask.mask(sv1_bgs_bits) != 0
                isGoodFiber = fibermap['FIBERSTATUS'] == 0
                isGoodZbest = (zbest['DELTACHI2'] > 25.) & (zbest['ZWARN'] == 0)
                select = isTGT & isGAL & isBGS & isGoodFiber & isGoodZbest

                print('     + selected: {}'.format(np.sum(select)))

                # Accumulate spectrum data.
                if (np.sum(select) > 0):
                    if allzbest is None:
                        allzbest = zbest[select]
                        allfmap = fibermap[select]
                        allwave = cspectra.wave['brz']
                        allflux = cspectra.flux['brz'][select]
                        allivar = cspectra.ivar['brz'][select]
                        allmask = cspectra.mask['brz'][select]
                        allres  = cspectra.resolution_data['brz'][select]
                    else:
                        allzbest = vstack([allzbest, zbest[select]])
                        allfmap = vstack([allfmap, fibermap[select]])
                        allflux = np.vstack([allflux, cspectra.flux['brz'][select]])
                        allivar = np.vstack([allivar, cspectra.ivar['brz'][select]])
                        allmask = np.vstack([allmask, cspectra.mask['brz'][select]])
                        allres  = np.vstack([allres, cspectra.resolution_data['brz'][select]])

        #                 # Break loop over petals - for now this is the easy way of checking if an exposure is 
        #                 break
                #fitsio.read_header(cafile)

            if allzbest is None:
                print("No useful BGS observations in this tile")
            else:
                # Apply the DESITRIP preprocessing to selected spectra.
                rewave, reflux, reivar = rebin_flux(allwave, allflux, allivar, allzbest['Z'],
                                                    minwave=2500., maxwave=9500., nbins=150,
                                                    log=True, clip=True)
                rsflux = rescale_flux(reflux)

                # Run the classifier on the spectra.
                # The output layer uses softmax activation to produce an array of label probabilities.
                # The classification is based on argmax(pred).
                pred = classifier.predict(rsflux)
                ymax = np.max(pred, axis=1)
                    
                ### Selection on Classifier Output
                # To be conservative we can select only spectra where the classifier is very confident in its output, e.g., ymax > 0.99. See the [CNN training notebook](https://github.com/desihub/timedomain/blob/master/desitrip/docs/nb/cnn_multilabel-restframe.ipynb) for the motivation behind this cut.
                idx = np.argwhere(ymax > 0.99).flatten()
                labels = np.argmax(pred, axis=1)
                ntr = len(idx)
                print(ntr,' transients found')
#                     print(idx)
                
                # Save data to file - we need to add ra dec, and some id
                
                if ntr > 0:
                #    if tr_z is None:
                #        tr_z = allzbest[idx]['Z']
                #        tr_label = label_names_arr[labels[idx]]
                #        tr_tile = np.full(ntr,tile_number)
                #        tr_date = np.full(ntr,obsdate)
                #        tr_spectrum = rsflux[idx]
                #    else:
                #        tr_z = np.vstack((tr_z,allzbest[idx]['Z']))
                #        tr_label = np.vstack((tr_label,label_names_arr[labels[idx]]))
                #        tr_tile = np.vstack((tr_tile,np.full(ntr,tile_number)))
                #        tr_date = np.vstack((tr_date,np.full(ntr,obsdate)))
                #        tr_spectrum = np.vstack((tr_spectrum,rsflux[idx]))

                    # Save classification info to a table.
                    classification = Table()
                    classification['TARGETID'] = allfmap[idx]['TARGETID']
                    classification['CNNPRED'] = pred[idx]
                    classification['CNNLABEL'] = label_names_arr[labels[idx]]

                    # Merge the classification and redrock fit to the fibermap.
                    fmap = join(allfmap[idx], allzbest[idx], keys='TARGETID')
                    fmap = join(fmap, classification, keys='TARGETID')

                    # Pack data into Spectra and write to FITS.
                    cand_spectra = Spectra(bands=['brz'],
                                           wave={'brz' : allwave},
                                           flux={'brz' : allflux[idx]},
                                           ivar={'brz' : allivar[idx]},
                                           resolution_data={'brz' : allres[idx]},
                                           fibermap=fmap
                                       )

                    
                    
                    outfits = '{}/transient_candidate_spectra_{}_{}.fits'.format(out_path, obsdate, tile_number)
                    write_spectra(outfits, cand_spectra)
                    print('Output file saved in {}'.format(outfits))
                    


                    #DESITRIP_daily - Send the selected fluxes to SkyPortal - Divij Sharma
                    for i in range(len(fmap['TARGETID'])):
                        SkyPortal.postCandidate(i, fmap)
                        SkyPortal.postSpectra(fmap['TARGETID'][i].astype('str'), cand_spectra)
            

                    # Make a plot of up to 16 transients

                    selection = sorted(np.random.choice(idx.flatten(), size=np.minimum(len(idx), 16), replace=False))

                    fig, axes = plt.subplots(4,4, figsize=(15,10), sharex=True, sharey=True,
                                             gridspec_kw={'wspace':0, 'hspace':0})

                    for j, ax in zip(selection, axes.flatten()):
                        ax.plot(rewave, rsflux[j], alpha=0.7, label='label: '+label_names[labels[j]]+'\nz={:.2f}'.format(allzbest[j]['Z']))
                        ax.legend(fontsize=10)

                    fig.tight_layout()
                    outplot = '{}/transient_candidates_{}_{}.png'.format(plot_path, obsdate, tile_number)
                    fig.savefig(outplot, dpi=200)
                    print('Figure saved in {}', outplot)
                    
                #Now add this tile info to the sqlite db
                #Maybe expid is important for spectraldiff? Here we coadd
                #expids = set([int(cframefile.split('-')[-1][:-5]) for cframefile in cframefiles])
                prog=cframe_fibermap['PROGRAM']
                sql = """ INSERT OR IGNORE INTO desitrip_exposures(tileid,program,obsdate)
              VALUES(?,?,?) """
                c.execute(sql,(tile_number, prog, obsdate))

            else:
                print('Not a BGS tile')

    if tr_z is None:
        print("No transients found on night ",obsdate)
   else:
       c1 = fits.Column(name='Z', array=tr_z, format='F')
       c2 = fits.Column(name='LABEL', array=tr_label, format='6A')
       c3 = fits.Column(name='TILE', array=tr_tile, format='10A')
       c4 = fits.Column(name='DATE', array=tr_date, format='10A')
       #Will think of a more clever way to save the spectra
       #c5 = fits.Column(name='SPECTRUM', array=tr_spectrum[0], format='F')
       t = fits.BinTableHDU.from_columns([c1, c2, c3, c4]) #, c5])
       t.writeto(out_path+'transients_'+obsdate+'.fits', overwrite=True)
       print('Output file saved in ', out_path)

#Close connection to sqlite
conn.commit()
conn.close()  