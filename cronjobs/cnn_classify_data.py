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
sys.path.append('/global/homes/p/palmese/desi/timedomain/desitrip/py/') #Note:change this path as needed!
sys.path.append('/global/homes/p/palmese/desi/timedomain/timedomain/')

from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra
from desispec.coaddition import coadd_cameras
from desitarget.sv1.sv1_targetmask import bgs_mask as bgs_mask_sv1
from desitarget.sv2.sv2_targetmask import bgs_mask as bgs_mask_sv2
from desitarget.sv3.sv3_targetmask import bgs_mask as bgs_mask_sv3

from desitrip.preproc import rebin_flux, rescale_flux
from desitrip.deltamag import delta_mag

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

#DESITRIP_daily
from timedomain.sp_utils import *
from timedomain.filters import *
#from timedomain.iterators import *

import sqlite3
import exposure_db

# Add redrock templates to look at residuals 

import redrock.templates

templates = dict()
for filename in redrock.templates.find_templates():
    t = redrock.templates.Template(filename)
    templates[(t.template_type, t.sub_type)] = t


#sqlite db path
db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'


if __name__ == '__main__':
    parser = ArgumentParser(description='Run CNN classifier on selected tile and day')

    parser.add_argument('-d', '--obsdate', type=str,default='20201221',
                        help='Night to be processed')
    parser.add_argument('-t', '--tilenum', nargs='+', type=str,default='80611',
                        help='Tile Numbers to be processed')
    parser.add_argument('-r', '--redux', default='daily',
                        help='Spectroscopic reduction: daily, andes, blanc, ...')
    parser.add_argument('-g', '--gradcam', default=True,
                        help='Set to true if you want to apply GRADCam to the plotting -takes more time')
    
    #If there are more than 1 obsdate, provide a 2D array
    parser.add_argument('-o', '--obsdates_tilenumbers', nargs='+', type=str,default=None,
                        help='str array with columns obsdate, tilenumber, separated by |')

    args = parser.parse_args()
    gradcam = args.gradcam
    
    if args.obsdates_tilenumbers!=None:
        obsdates_tilenumbers_str = args.obsdates_tilenumbers
        obsdates_tilenumbers = np.chararray((len(obsdates_tilenumbers_str),2),itemsize=10,unicode=True)
        for i in range(len(obsdates_tilenumbers_str)):
            obsdates_tilenumbers[i,:]=obsdates_tilenumbers_str[i].split('|')
        print(obsdates_tilenumbers_str,obsdates_tilenumbers)
    else:    
        obsdates = args.obsdate
        tilenums = args.tilenum
        if type(obsdates) is str: nobs=1
        else: nobs=len(obsdates)
        if type(tilenums) is str: ntiles=1
        else: ntiles=len(tilenums)
        obsdates_tilenumbers = np.chararray((ntiles,2),itemsize=10,unicode=True)
        obsdates_tilenumbers[:,0] =obsdates
        obsdates_tilenumbers[:,1] = tilenums#np.chararray((nobs,ntiles),itemsize=10)
        #obsdates_tilenumbers[:,0]=obsdates
        #obsdates_tilenumbers[:,1]=tilenums
        
    base_path='/global/u2/p/palmese/desi/timedomain/cronjobs/'
    td_path='/global/cfs/cdirs/desi/science/td/daily-search/desitrip/'
    plot_path=td_path+'plots/'
    out_path=td_path+'out/'
    # Set up BGS target bit selection.
    sv1_bgs_bits = '|'.join([_ for _ in bgs_mask_sv1.names() if 'BGS' in _])
    sv2_bgs_bits = '|'.join([_ for _ in bgs_mask_sv2.names() if 'BGS' in _])
    sv3_bgs_bits = '|'.join([_ for _ in bgs_mask_sv3.names() if 'BGS' in _])

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
    
    #Only import tf and define stuff if you want to applt gradcam - otherwise save time
    #This function needs to be changed if the keras model architecture is changed
    
    if gradcam:
        import tensorflow as tf
        last_conv_layer_name = "conv1d_23"
        classifier_layer_names = [
        "batch_normalization_23",
        "activation_23",
        "max_pooling1d_23",
        "flatten_5",
        "dense_5",
        "dropout_5",
        "Output_Classes"
        ]
        
        def make_gradcam_heatmap(
            img_array, model, last_conv_layer_name, classifier_layer_names
        ):
            # First, we create a model that maps the input image to the activations
            # of the last conv layer
            last_conv_layer = model.get_layer(last_conv_layer_name)
            last_conv_layer_model = keras.Model(model.inputs, last_conv_layer.output)

            # Second, we create a model that maps the activations of the last conv
            # layer to the final class predictions
            classifier_input = keras.Input(shape=last_conv_layer.output.shape[1:])
            x = classifier_input
            for layer_name in classifier_layer_names:
                #print(layer_name,x.shape)
                x = model.get_layer(layer_name)(x)
            classifier_model = keras.Model(classifier_input, x)

            # Then, we compute the gradient of the top predicted class for our input image
            # with respect to the activations of the last conv layer
            with tf.GradientTape() as tape:
                # Compute activations of the last conv layer and make the tape watch it
                last_conv_layer_output = last_conv_layer_model(img_array)
                tape.watch(last_conv_layer_output)
                # Compute class predictions
                preds = classifier_model(last_conv_layer_output)
                top_pred_index = tf.argmax(preds[0])
                top_class_channel = preds[:, top_pred_index]

            # This is the gradient of the top predicted class with regard to
            # the output feature map of the last conv layer
            grads = tape.gradient(top_class_channel, last_conv_layer_output)
            # This is a vector where each entry is the mean intensity of the gradient
            # over a specific feature map channel
            pooled_grads = tf.reduce_mean(grads, axis=(0, 1))
            #print(grads.shape,pooled_grads.shape)

            # We multiply each channel in the feature map array
            # by "how important this channel is" with regard to the top predicted class
            last_conv_layer_output = last_conv_layer_output.numpy()[0]
            pooled_grads = pooled_grads.numpy()
            for i in range(pooled_grads.shape[-1]):
                last_conv_layer_output[:, i] *= pooled_grads[i]

            # The channel-wise mean of the resulting feature map
            # is our heatmap of class activation
            heatmap = np.mean(last_conv_layer_output, axis=-1)

            #We apply ReLU here and select only elements>0
            # For visualization purpose, we will also normalize the heatmap between 0 & 1
            heatmap = np.maximum(heatmap, 0) / np.max(heatmap)
            return heatmap
        
        preprocess_input = keras.applications.xception.preprocess_input
        decode_predictions = keras.applications.xception.decode_predictions
        
        
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
    for obsdate,tile_number in obsdates_tilenumbers:
        #20210503 is when the directory structure changed, no more daily reductions in the tiles directory
        redux = '/'.join([os.environ['DESI_SPECTRO_REDUX'], args.redux])
        if int(obsdate)<20210503:
            prefix_in = '/'.join([redux, 'tiles', tile_number, obsdate])
        else:
            prefix_in = '/'.join([redux, 'tiles/cumulative', tile_number])
            
            #Keeping these lines for later - they are useful if we want to use the cframe files in the future
            #We need to search in the db the expid for that night      
            #db_filename = '/global/cfs/cdirs/desi/science/td/daily-search/transients_search.db'
            #db = exposure_db.ReduxDB(db_filename)
            #conn = sqlite3.connect(db_filename)
            #c = conn.cursor()    
            #check_expids_query="select expid from exposures where tileid="+tile_number+" and obsdate="+obsdate+";"
            #cursor = c.execute(check_expids_query)
            #expids_list = cursor.fetchall()
            #Only pick first one if multiple expids
            #this_expid  = '000'+str(expids_list[0][0])
            #conn.close()          
            #prefix_in = '/'.join([redux, 'exposures', obsdate,this_expid])
            
        if not os.path.isdir(prefix_in):
            print('{} does not exist.'.format(prefix_in))
        else:
            #Since cframe files do not live there anymore, commenting this out
            # We have already checked in cron_db_clever that this is a bgs tile
            
            #cframefiles = sorted(glob('{}/cframe*.fits'.format(prefix_in)))
            #cframe_fibermap = fitsio.read_header(cframefiles[0],'FIBERMAP')
                            #Check which survey it is - sv1 sv2 etc
            #if (cframe_fibermap['FA_SURV']=='sv1'): 
            #    bgs_mask=bgs_mask_sv1
            #    bgs_bits=sv1_bgs_bits
            #if (cframe_fibermap['FA_SURV']=='sv2'): 
            #    bgs_mask=bgs_mask_sv2
            #    bgs_bits=sv2_bgs_bits
            #if (cframe_fibermap['FA_SURV']=='sv3'): 
            #    bgs_mask=bgs_mask_sv3
            #    bgs_bits=sv3_bgs_bits
            
            if int(obsdate)>20210404:
                bgs_mask=bgs_mask_sv3
                bgs_bits=sv3_bgs_bits
            else:
                bgs_mask=bgs_mask_sv2
                bgs_bits=sv2_bgs_bits

            # Access the zbest and coadd files.
            # Files are organized by petal number.
            if int(obsdate)<20210503:
                zbfiles = sorted(glob('{}/zbest*.fits'.format(prefix_in)))
                cafiles = sorted(glob('{}/coadd*.fits'.format(prefix_in)))
            else:
                #Only grabbing stuff from the most recent cumulative directory
                zb_ca_dir = sorted(glob('{}/*/'.format(prefix_in)))[-1]
                zbfiles = sorted(glob('{}/zbest*.fits'.format(zb_ca_dir)))
                cafiles = sorted(glob('{}/spectra*.fits'.format(zb_ca_dir)))
                #cafiles = sorted(glob('{}/cframe*.fits'.format(prefix_in)))

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
                if int(obsdate)>20210503:
                    select_nite = pspectra.fibermap['NIGHT'] == int(obsdate)
                    pspectra = pspectra[select_nite]

                #Should we move the coaddition to after the spectra selection?
                #It would require to only select a single night first though.
                cspectra = coadd_cameras(pspectra)
                fibermap = cspectra.fibermap

                # Apply standard event selection.
                isTGT = fibermap['OBJTYPE'] == 'TGT'
                isGAL = zbest['SPECTYPE'] == 'GALAXY'

                #This is old selection, does not work anymore! BGS target is always 0
                #isBGS = fibermap['SV1_BGS_TARGET'] & bgs_mask.mask(sv1_bgs_bits) != 0
                isBGS = bgs_mask.mask(bgs_bits) != 0
                isGoodFiber = fibermap['FIBERSTATUS'] == 0
                isGoodZbest = (zbest['DELTACHI2'] > 25.) & (zbest['ZWARN'] == 0)
                select = isTGT & isGAL & isBGS & isGoodFiber & isGoodZbest
                
                fibermap = delta_mag(cspectra, fibermap, select, nsigma=3)

                print('     + selected: {}'.format(np.sum(select)))

                # Accumulate spectrum data.
                arm_avail = [key for key in cspectra.wave.keys()][0]
                if (np.sum(select) > 0):
                    if allzbest is None:
                        allzbest = zbest[select]
                        allfmap = fibermap[select]
                        allwave = cspectra.wave[arm_avail]
                        allflux = cspectra.flux[arm_avail][select]
                        allivar = cspectra.ivar[arm_avail][select]
                        allmask = cspectra.mask[arm_avail][select]
                        allres  = cspectra.resolution_data[arm_avail][select]
                    else:
                        allzbest = vstack([allzbest, zbest[select]])
                        allfmap = vstack([allfmap, fibermap[select]])
                        allflux = np.vstack([allflux, cspectra.flux[arm_avail][select]])
                        allivar = np.vstack([allivar, cspectra.ivar[arm_avail][select]])
                        allmask = np.vstack([allmask, cspectra.mask[arm_avail][select]])
                        allres  = np.vstack([allres, cspectra.resolution_data[arm_avail][select]])

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
                    cand_spectra = Spectra(bands=[arm_avail],
                                           wave={arm_avail : allwave},
                                           flux={arm_avail : allflux[idx]},
                                           ivar={arm_avail : allivar[idx]},
                                           resolution_data={arm_avail : allres[idx]},
                                           fibermap=fmap
                                       )

                    outfits = '{}/transient_candidate_spectra_{}_{}.fits'.format(out_path, obsdate, tile_number)
                    write_spectra(outfits, cand_spectra)
                    print('Output file saved in {}'.format(outfits))

                    #DESITRIP_daily - Send the selected fluxes to SkyPortal - Divij Sharma
                    #Modified by AP in Apr 2021 to fix new error
                    for i in range(len(fmap['TARGETID'])):
                        # Collect the extra data specific to DESITRIP to be saved
                        data = dict()
                        data['redshift'] = fmap['Z'].data[i]
                        altdata_dict = dict()
                        altdata_dict['classifier']={'CNNLABEL' : fmap['CNNLABEL'].data[i]}
                        data['altdata'] = altdata_dict
                            
                        try:
                            SkyPortal.postCandidate(i, fmap, 'DESITRIP', data_override=data)
                        except UniqueViolationError as err:
                            log.info(err.response.json())
                        except Exception as err:
                            print(err)
#                             sys.exit(1)

                        try:
                            SkyPortal.postSpectra(fmap['TARGETID'][i], cand_spectra)
                            #Also post the residuals and redrock to skyportal!
                
                        except Exception as err:
                            print(err)

                    # Make a plot of up to 16 transients

                    selection = sorted(np.random.choice(idx.flatten(), size=np.minimum(len(idx), 16), replace=False))

                    fig, axes = plt.subplots(4,4, figsize=(15,10), sharex=True, sharey=True,
                                             gridspec_kw={'wspace':0, 'hspace':0})

                    #these lines are to add to the output plot the wavelengths between the arms
                    br_band=[5600,6000]
                    rz_band=[7400,7800]

                    for j, ax in zip(selection, axes.flatten()):

                        delta_fibermag_g=None
                        delta_fibermag_r=None
                        delta_fibermag_z=None
                        if 'b' in arm_avail:
                            delta_fibermag_g=allfmap['DELTAMAG_G'][j]
                        if 'r' in arm_avail:
                            delta_fibermag_r=allfmap['DELTAMAG_R'][j]
                        if 'z' in arm_avail:
                            delta_fibermag_z=allfmap['DELTAMAG_Z'][j]

                        if gradcam:
                            this_flux=rsflux[j,:].reshape((1,150)) 

                            # Generate class activation heatmap
                            heatmap = make_gradcam_heatmap(
                                this_flux, classifier, last_conv_layer_name, classifier_layer_names
                            )

                            color='blue'
                            
                            #If a candidate has a deltamag>2 sigma away from the mean deltamag distribution
                            #it's flagged as a deltamag candidate and we plot it in red
                            if allfmap['DELTAMAG_CANDIDATE'][j]:
                                print("FOUND A DELTAMAG and DESITRIP CANDIDATE!")
                                color='red'
                            rewave_nbin_inblock=rewave.shape[0]/float(heatmap.shape[0])
                            first_bin=0
                            for i in range(1,heatmap.shape[0]+1):
                                alpha=np.min([1,heatmap[i-1]+0.2])
                                last_bin=int(i*rewave_nbin_inblock)
                                if (i==1):
                                    ax.plot(rewave[first_bin:last_bin+1], this_flux[0,first_bin:last_bin+1],c=color,alpha=alpha,\
                                            label=str(allfmap['TARGETID'][j])+'\n'+label_names[labels[j]]+'\nz={:.2f}'.format(allzbest[j]['Z'])\
                                           +'\n dg={:.2f}'.format(delta_fibermag_g)\
                                            +'\n dr={:.2f}'.format(delta_fibermag_r)\
                                           +'\n dz={:.2f}'.format(delta_fibermag_z)\
                                           )
                                else:
                                    ax.plot(rewave[first_bin:last_bin+1], this_flux[0,first_bin:last_bin+1],c=color,alpha=alpha)
                                first_bin=last_bin

                        else:
                            ax.plot(rewave, rsflux[j], alpha=0.7, label='label: '+label_names[labels[j]]+'\nz={:.2f}'.format(allzbest[j]['Z']))
                        this_br_band=br_band/(1.+(allzbest[j]['Z']))
                        this_rz_band=rz_band/(1.+(allzbest[j]['Z']))  
                        ax.fill_between(this_br_band,[0,0],[1,1],alpha=0.1,color='k')
                        ax.fill_between(this_rz_band,[0,0],[1,1],alpha=0.1,color='k')                      
                        ax.legend(fontsize=10)

                    fig.tight_layout()
                    outplot = '{}/transient_candidates_{}_{}.png'.format(plot_path, obsdate, tile_number)
                    fig.savefig(outplot, dpi=200)
                    print('Figure saved in {}', outplot)
                    
                    #Here plot redrock residuals
                    #Turning off residuals plots until it works
                    plot_residuals=False
                    
                    if plot_residuals:
                        fig, axes = plt.subplots(4,4, figsize=(15,10), sharex=True, sharey=True,
                                        gridspec_kw={'wspace':0, 'hspace':0})

                        for j, ax in zip(selection, axes.flatten()):

                            spectype = allzbest['SPECTYPE'][j].strip()
                            subtype = allzbest['SUBTYPE'][j].strip()
                            fulltype = (spectype, subtype)
                            ncoeff = templates[fulltype].flux.shape[0]
                            coeff = allzbest['COEFF'][j][0:ncoeff]
                            tflux = templates[fulltype].flux.T.dot(coeff)
                            twave = templates[fulltype].wave * (1+allzbest[j]['Z'])

                            plt.plot(allwave[j], allflux[j]-tflux, 'k-', alpha=0.5)
                            #maxflux = max(maxflux, np.max(spectra.flux[band][ispec]))

                        fig.tight_layout()
                        outplot = '{}/residual_candidates_{}_{}.png'.format(plot_path, obsdate, tile_number)
                        fig.savefig(outplot, dpi=200)
                        print('Figure saved in {}', outplot)

                #Now add this tile info to the sqlite db
                #Maybe expid is important for spectraldiff? Here we coadd
                #expids = set([int(cframefile.split('-')[-1][:-5]) for cframefile in cframefiles])
                #Open db only here
                conn = sqlite3.connect(db_filename)
                c = conn.cursor()
                #prog=cframe_fibermap['PROGRAM']
                prog='BGS'
                sql = """ INSERT OR IGNORE INTO desitrip_exposures(tileid,program,obsdate) VALUES(?,?,?) """
                c.execute(sql,(tile_number, prog, obsdate))
                #Close connection to sqlite
                #We open and close in the loop to minimize the amount of time
                #the db is open to minimize clashes with outer jobs opening the
                #same db
                conn.commit()
                conn.close()  

        #else:
        #    print('Not a BGS tile, program: ',cframe_fibermap['PROGRAM'])


