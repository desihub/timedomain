"""
Docs » Module code » desispec.coaddition
Source code for desispec.coaddition

Coadd spectra
"""

from __future__ import absolute_import, division, print_function
import os, sys, time

import numpy as np

import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg

from astropy.table import Table, Column

import multiprocessing

from desiutil.log import get_logger

from desispec.interpolation import resample_flux
from desispec.spectra import Spectra
from desispec.resolution import Resolution
from desispec.fiberbitmasking import get_all_fiberbitmask_with_amp, get_all_nonamp_fiberbitmask_val, get_justamps_fiberbitmask
from desispec.specscore import compute_coadd_scores

class spectra_no_expid(Spectra):
    def show(self):
        self.fibermap
        return
    

    def coadd_fibermap_no_expid(fibermap) :

        log = get_logger()
        log.debug("'coadding' fibermap")

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
    
    
    
    def coadd_no_expid(spectra_no_expid, cosmics_nsig=0.) :
        """
        Coaddition the spectra for each target and each camera. The input spectra is modified.

        Args:
           spectra: desispec.spectra.Spectra object

        Options:
           cosmics_nsig: float, nsigma clipping threshold for cosmics rays
        """
        log = get_logger()
        targets = np.unique(spectra_no_expid.fibermap["TARGETID"])
        ntarget=targets.size
        log.debug("number of targets= {}".format(ntarget))
        for b in spectra_no_expid.bands :
            log.debug("coadding band '{}'".format(b))
            nwave=spectra_no_expid.wave[b].size
            tflux=np.zeros((ntarget,nwave),dtype=spectra_no_expid.flux[b].dtype)
            tivar=np.zeros((ntarget,nwave),dtype=spectra_no_expid.ivar[b].dtype)
            if spectra_no_expid.mask is not None :
                tmask=np.zeros((ntarget,nwave),dtype=spectra_no_expid.mask[b].dtype)
            else :
                tmask=None
            trdata=np.zeros((ntarget,spectra_no_expid.resolution_data[b].shape[1],nwave),dtype=spectra_no_expid.resolution_data[b].dtype)

            fiberstatus_bits = get_all_fiberbitmask_with_amp(b)
            good_fiberstatus = ( (spectra_no_expid.fibermap["FIBERSTATUS"] & fiberstatus_bits) == 0 )
            for i,tid in enumerate(targets) :
                jj=np.where( (spectra_no_expid.fibermap["TARGETID"]==tid) & good_fiberstatus )[0]

                #- if all spectra were flagged as bad (FIBERSTATUS != 0), contine
                #- to next target, leaving tflux and tivar=0 for this target
                if len(jj) == 0:
                    continue

                if cosmics_nsig is not None and cosmics_nsig > 0  and len(jj)>2 :
                    # interpolate over bad measurements
                    # to be able to compute gradient next
                    # to a bad pixel and identify outlier
                    # many cosmics residuals are on edge
                    # of cosmic ray trace, and so can be
                    # next to a masked flux bin
                    grad=[]
                    gradvar=[]
                    for j in jj :
                        if spectra_no_expid.mask is not None :
                            ttivar = spectra_no_expid.ivar[b][j]*(spectra_no_expid.mask[b][j]==0)
                        else :
                            ttivar = spectra_no_expid.ivar[b][j]
                        good = (ttivar>0)
                        bad  = ~good
                        if np.sum(good)==0 :
                            continue
                        nbad = np.sum(bad)
                        ttflux = spectra_no_expid.flux[b][j].copy()
                        if nbad>0 :
                            ttflux[bad] = np.interp(spectra_no_expid.wave[b][bad],spectra_no_expid.wave[b][good],ttflux[good])
                        ttivar = spectra_no_expid.ivar[b][j].copy()
                        if nbad>0 :
                            ttivar[bad] = np.interp(spectra_no_expid.wave[b][bad],spectra_no_expid.wave[b][good],ttivar[good])
                        ttvar = 1./(ttivar+(ttivar==0))
                        ttflux[1:] = ttflux[1:]-ttflux[:-1]
                        ttvar[1:]  = ttvar[1:]+ttvar[:-1]
                        ttflux[0]  = 0
                        grad.append(ttflux)
                        gradvar.append(ttvar)

                tivar_unmasked= np.sum(spectra_no_expid.ivar[b][jj],axis=0)
                if spectra_no_expid.mask is not None :
                    ivarjj=spectra_no_expid.ivar[b][jj]*(spectra_no_expid.mask[b][jj]==0)
                else :
                    ivarjj=spectra_no_expid.ivar[b][jj]
                if cosmics_nsig is not None and cosmics_nsig > 0 and len(jj)>2  :
                    grad=np.array(grad)
                    gradvar=np.array(gradvar)
                    gradivar=(gradvar>0)/np.array(gradvar+(gradvar==0))
                    nspec=grad.shape[0]
                    sgradivar=np.sum(gradivar)
                    if sgradivar>0 :
                        meangrad=np.sum(gradivar*grad,axis=0)/sgradivar
                        deltagrad=grad-meangrad
                        chi2=np.sum(gradivar*deltagrad**2,axis=0)/(nspec-1)

                        bad  = (chi2>cosmics_nsig**2)
                        nbad = np.sum(bad)
                        if nbad>0 :
                            log.info("masking {} values for targetid={}".format(nbad,tid))
                            badindex=np.where(bad)[0]
                            for bi in badindex  :
                                k=np.argmax(gradivar[:,bi]*deltagrad[:,bi]**2)
                                ivarjj[k,bi]=0.
                                log.debug("masking spec {} wave={}".format(k,spectra_no_expid.wave[b][bi]))

                tivar[i]=np.sum(ivarjj,axis=0)
                tflux[i]=np.sum(ivarjj*spectra_no_expid.flux[b][jj],axis=0)
                for r in range(spectra_no_expid.resolution_data[b].shape[1]) :
                    trdata[i,r]=np.sum((spectra_no_expid.ivar[b][jj]*spectra_no_expid.resolution_data[b][jj,r]),axis=0) # not sure applying mask is wise here
                bad=(tivar[i]==0)
                if np.sum(bad)>0 :
                    tivar[i][bad] = np.sum(spectra_no_expid.ivar[b][jj][:,bad],axis=0) # if all masked, keep original ivar
                    tflux[i][bad] = np.sum(spectra_no_expid.ivar[b][jj][:,bad]*spectra_no_expid.flux[b][jj][:,bad],axis=0)
                ok=(tivar[i]>0)
                if np.sum(ok)>0 :
                    tflux[i][ok] /= tivar[i][ok]
                ok=(tivar_unmasked>0)
                if np.sum(ok)>0 :
                    trdata[i][:,ok] /= tivar_unmasked[ok]
                if spectra_no_expid.mask is not None :
                    tmask[i]      = np.bitwise_and.reduce(spectra_no_expid.mask[b][jj],axis=0)
            spectra_no_expid.flux[b] = tflux
            spectra_no_expid.ivar[b] = tivar
            if spectra_no_expid.mask is not None :
                spectra_no_expid.mask[b] = tmask
            spectra_no_expid.resolution_data[b] = trdata

        if spectra_no_expid.scores is not None:
            orig_scores = Table(spectra_no_expid.scores.copy())
            orig_scores['TARGETID'] = spectra_no_expid.fibermap['TARGETID']
        else:
            orig_scores = None

        spectra_no_expid.fibermap=coadd_fibermap_no_expid(spectra_no_expid.fibermap)
        spectra_no_expid.scores=None
        compute_coadd_scores(spectra_no_expid, orig_scores, update_coadd=True)
        return