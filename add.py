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
from desispec.coaddition import coadd_fibermap








def add(spectra, cosmics_nsig=0.) :
    """
    Coaddition the spectra for each target and each camera. The input spectra is modified.

    Args:
       spectra: desispec.spectra.Spectra object

    Options:
       cosmics_nsig: float, nsigma clipping threshold for cosmics rays
    """
    log = get_logger()
    targets = np.unique(spectra.fibermap["TARGETID"])
    ntarget=targets.size
    log.debug("number of targets= {}".format(ntarget))
    for b in spectra.bands :
        log.debug("coadding band '{}'".format(b))
        nwave=spectra.wave[b].size
        tflux=np.zeros((ntarget,nwave),dtype=spectra.flux[b].dtype)
        tivar=np.zeros((ntarget,nwave),dtype=spectra.ivar[b].dtype)
        if spectra.mask is not None :
            tmask=np.zeros((ntarget,nwave),dtype=spectra.mask[b].dtype)
        else :
            tmask=None
        trdata=np.zeros((ntarget,spectra.resolution_data[b].shape[1],nwave),dtype=spectra.resolution_data[b].dtype)

        fiberstatus_bits = get_all_fiberbitmask_with_amp(b)
        good_fiberstatus = ( (spectra.fibermap["FIBERSTATUS"] & fiberstatus_bits) == 0 )
        for i,tid in enumerate(targets) :
            jj=np.where( (spectra.fibermap["TARGETID"]==tid) & good_fiberstatus )[0]

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
                    if spectra.mask is not None :
                        ttivar = spectra.ivar[b][j]*(spectra.mask[b][j]==0)
                    else :
                        ttivar = spectra.ivar[b][j]
                    good = (ttivar>0)
                    bad  = ~good
                    if np.sum(good)==0 :
                        continue
                    nbad = np.sum(bad)
                    ttflux = spectra.flux[b][j].copy()
                    if nbad>0 :
                        ttflux[bad] = np.interp(spectra.wave[b][bad],spectra.wave[b][good],ttflux[good])
                    ttivar = spectra.ivar[b][j].copy()
                    if nbad>0 :
                        ttivar[bad] = np.interp(spectra.wave[b][bad],spectra.wave[b][good],ttivar[good])
                    ttvar = 1./(ttivar+(ttivar==0))
                    ttflux[1:] = ttflux[1:]-ttflux[:-1]
                    ttvar[1:]  = ttvar[1:]+ttvar[:-1]
                    ttflux[0]  = 0
                    grad.append(ttflux)
                    gradvar.append(ttvar)

            #tivar_unmasked= np.sum(spectra.ivar[b][jj],axis=0)
            tivar_unmasked = 1 / np.sum(1/spectra.ivar[b][jj],axis=0)
            if spectra.mask is not None :
                ivarjj=spectra.ivar[b][jj]*(spectra.mask[b][jj]==0)
            else :
                ivarjj=spectra.ivar[b][jj]
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
                            log.debug("masking spec {} wave={}".format(k,spectra.wave[b][bi]))

            #tivar[i]=np.sum(ivarjj,axis=0)
            tivar[i]= 1 / np.sum(1/ivarjj,axis=0)
            tflux[i]=np.sum(spectra.flux[b][jj],axis=0)
            for r in range(spectra.resolution_data[b].shape[1]) :
                trdata[i,r]=np.sum((spectra.resolution_data[b][jj,r]),axis=0) # not sure applying mask is wise here
            bad=(tivar[i]==0)
            if np.sum(bad)>0 :
                tivar[i][bad] = 1 / np.sum(1/spectra.ivar[b][jj][:,bad],axis=0) # if all masked, keep original ivar
                tflux[i][bad] = np.sum(spectra.flux[b][jj][:,bad],axis=0)
            ok=(tivar[i]>0)
            #if np.sum(ok)>0 :
                #tflux[i][ok] /= tivar[i][ok]
            ok=(tivar_unmasked>0)
            if np.sum(ok)>0 :
                trdata[i][:,ok] /= tivar_unmasked[ok]
            if spectra.mask is not None :
                tmask[i]      = np.bitwise_or.reduce(spectra.mask[b][jj],axis=0)
        spectra.flux[b] = tflux
        spectra.ivar[b] = tivar
        if spectra.mask is not None :
            spectra.mask[b] = tmask
        spectra.resolution_data[b] = trdata

    if spectra.scores is not None:
        orig_scores = Table(spectra.scores.copy())
        orig_scores['TARGETID'] = spectra.fibermap['TARGETID']
    else:
        orig_scores = None

    spectra.fibermap=coadd_fibermap(spectra.fibermap)
    spectra.scores=None
    compute_coadd_scores(spectra, orig_scores, update_coadd=True)
