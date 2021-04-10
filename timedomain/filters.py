import numpy as np
import numpy.ma as ma

# from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

import matplotlib.pyplot as plt
from desiutil.log import get_logger
log = get_logger()
from . import plot_utils

# difference between two spectra
# for now the two spectra are assumed to be single night coadd of matching tile

def maskskylines(s):
    skylam = [4340.,4368]
    R=1000.
    mask = dict()
    for dindex in s.bands:
        ans = np.zeros(len(s.wave[dindex]))
        for wa in skylam:
            wmin = wa * np.exp(-1/R)
            wmax = wa * np.exp(1/R)
            ans = np.logical_or(ans, np.logical_and(s.wave[dindex] >= wmin, s.wave[dindex] < wmax))
        mask[dindex]=ans

    return mask
        

    # s0-s1
def difference(s0, s1):
    diff = dict()
    ivar = dict()
    mask = dict()
    wave = dict()
    common = list(set(s1.bands).intersection(s0.bands))
    for dindex in common:

#         diff[dindex] = ma.array(data=s1.flux[dindex],mask=s1.mask[dindex])
#         diff[dindex] = diff[dindex] - ma.array(data=s0.flux[dindex],mask=s0.mask[dindex])
        
        diff[dindex] = s0.flux[dindex]-s1.flux[dindex]
        
#         ivar0 = ma.array(data=s0.ivar[dindex],mask=s0.mask[dindex])
#         ivar1 = ma.array(data=s1.ivar[dindex],mask=s1.mask[dindex])
#         ivar[dindex] = 1/(1/ivar0 + 1/ivar1)

        ivar[dindex] = 1/(1/s0.ivar[dindex] + 1/s1.ivar[dindex])

        mask[dindex] = 1-(1-s0.mask[dindex])*(1-s1.mask[dindex]).astype('int')
        
        wave[dindex] = s1.wave[dindex]
        
    ans = Spectra(bands=common, wave=wave, flux=diff,ivar=ivar,mask=mask,fibermap = s0.fibermap, resolution_data=s0.resolution_data)
    return ans

# renormalize two spectra
# for now the two spectra are assumed to be single night coadd of matching tile

def renorm(s0, s1):

    common = list(set(s1.bands).intersection(s0.bands))
    for dindex in common:
        norm = ma.array(data=s1.flux[dindex],mask=s1.mask[dindex])/ ma.array(data=s0.flux[dindex],mask=s0.mask[dindex])
        norm.filled(np.nan)
        norm = np.nanpercentile(norm,(50),axis=1)
        
        s1.flux[dindex] = s1.flux[dindex]/norm[:,None]
        s1.ivar[dindex] = s1.ivar[dindex]*((norm*norm)[:,None])

    return s0,s1


# Sometimes the spectrum has very low signal, perhaps due to mispointed fiber

class HasSignal:
    
    @staticmethod
    def filter(s0,s1,ston_cut=20.,anyband=True):
        ans = np.full(s0.num_spectra(),False)
        common = list(set(s1.bands).intersection(s0.bands))
        
        for i, dindex in enumerate(common):
            s1flux = ma.array(data=s1.flux[dindex],mask=s1.mask[dindex])
            s0flux = ma.array(data=s0.flux[dindex],mask=s0.mask[dindex])

            s1ivar = ma.array(data=s1.ivar[dindex],mask=s1.mask[dindex])
            s0ivar = ma.array(data=s0.ivar[dindex],mask=s0.mask[dindex])

            ston1 = s1flux.sum(axis=1) / ma.sqrt((1/s1ivar).sum(axis=1))
            ston0 = s0flux.sum(axis=1) / ma.sqrt((1/s0ivar).sum(axis=1))

            if i==0:
                ans = np.logical_and(ston0 > ston_cut, ston1 > ston_cut)
            else:
                if anyband:
                    ans = np.logical_or(ans, np.logical_and(ston0 > ston_cut, ston1 > ston_cut))
                else:
                    ans = np.logical_and(ans, np.logical_and(ston0 > ston_cut, ston1 > ston_cut))

        return ans
    

# Cataclysmic Variable

class CVLogic:
    target_wave = (6562.79, 4861.35, 4340.472, 4101.734, 3970.075)
    R=1000.
    plotter = plot_utils.diffplot_CV
    
    @staticmethod
    def filter(pspectra0, pspectra1, zbest, norm=True, ston_cut=7., frac_inc_cut= .20):
        fibermap = pspectra0.fibermap #Table.read(datafile0, 'FIBERMAP')
        isTGT = fibermap['OBJTYPE'] == 'TGT'
        okFibers = np.logical_and(pspectra0.fibermap['FIBERSTATUS'] == 0, pspectra1.fibermap['FIBERSTATUS'] == 0)
        isStar = zbest['SPECTYPE']=='STAR'
        
        # the interesting guy is wedge 6 index 331
    #     print(fibermap['TARGETID'].data[0], np.where(fibermap['TARGETID'].data== 35191288252861933)[0])
    #     pspectra0 = read_spectra(datafile0)
    #     pspectra1 = read_spectra(datafile1)
    
        hasSignal = HasSignal.filter(pspectra0,pspectra1)
        if norm:
            pspectra0, pspectra1 = renorm(pspectra0,pspectra1)
        diff = difference(pspectra0,pspectra1)
        
        skymask = maskskylines(diff)
        
        nspec =diff.num_spectra()
        
        signal=np.zeros(nspec)
        var=np.zeros(nspec)
        ref_signal=np.zeros(nspec)

        for dindex in diff.bands:
            
            #mask containing lines of interest
            lmask = np.zeros(len(diff.wave[dindex]))
            
            for wa in CVLogic.target_wave:
                wmin = wa * np.exp(-1/CVLogic.R/2.)
                wmax = wa * np.exp(1/CVLogic.R/2.)
                lmask = np.logical_or(lmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))
                
            #remove lines that are by the sky lines
            lmask = np.logical_and(lmask, skymask[dindex] ==0)

            for sindex in range(nspec):
                # only include unmasked
                nmask = np.logical_and(lmask, diff.mask[dindex][sindex,:]==0)
                signal[sindex] += diff.flux[dindex][sindex,nmask].sum()
                var[sindex] += (1/diff.ivar[dindex][sindex,nmask]).sum()
                ref_signal[sindex] += pspectra1.flux[dindex][sindex,nmask].sum()

        brighter = np.logical_or(np.abs(signal/ref_signal) >= frac_inc_cut, ref_signal <=0)
        significant = (np.abs(signal)/ma.sqrt(var) >= ston_cut)
        triggered = np.logical_and.reduce((significant, isTGT, hasSignal, okFibers, isStar, brighter))
        return triggered, diff
    
# TDE

class TDELogic:
    plotter = plot_utils.diffplot_CV
    
    @staticmethod
    def filter(pspectra0, pspectra1, zbest, norm=False, ston_cut=7., frac_inc_cut= .50):

        fibermap = pspectra0.fibermap #Table.read(datafile0, 'FIBERMAP')
        isTGT = fibermap['OBJTYPE'] == 'TGT'
        okFibers = np.logical_and(pspectra0.fibermap['FIBERSTATUS'] == 0, pspectra1.fibermap['FIBERSTATUS'] == 0)
        isGalaxy = zbest['SPECTYPE']=='GALAXY'
        
        hasSignal = HasSignal.filter(pspectra0,pspectra1)
        if norm:
            pspectra0, pspectra1 = renorm(pspectra0,pspectra1)
        diff = difference(pspectra0,pspectra1)
        
        skymask = maskskylines(diff)
        
        nspec = diff.flux[diff.bands[0]].shape[0]
        
        signal=np.zeros(nspec)
        var=np.zeros(nspec)
        ref_signal=np.zeros(nspec)
        
        if 'b' in diff.bands:           
            for sindex in range(nspec):
                nmask = diff.mask['b'][sindex,:]==0
                signal[sindex] = diff.flux['b'][sindex,nmask].sum()
                var[sindex] = (1/diff.ivar['b'][sindex,nmask]).sum()
                ref_signal[sindex] = pspectra1.flux['b'][sindex,nmask].sum()
        
        # nan's should fail here
        brighter = np.logical_or(np.abs(signal/ref_signal) >= frac_inc_cut, ref_signal <=0)
        significant = (np.abs(signal)/ma.sqrt(var) >= ston_cut)
        triggered = np.logical_and.reduce((significant, isTGT, hasSignal, okFibers, isGalaxy, brighter))
        return triggered, diff
    
# Cataclysmic Variable

class HydrogenLogic:
    target_wave = np.array((6562.79, 4861.35, 4340.472, 4101.734, 3970.075))
    R=200.
#     RH=1000     # cutting a notch out right at the line
    plotter = plot_utils.diffplot_CV
    maxwave = 7400.
    zmin=0.001
    
    @staticmethod
    def filter(pspectra0, pspectra1, zbest, norm=True, ston_cut=7., frac_inc_cut= .25):

        fibermap = pspectra0.fibermap #Table.read(datafile0, 'FIBERMAP')
        isTGT = fibermap['OBJTYPE'] == 'TGT'
        okFibers = np.logical_and(pspectra0.fibermap['FIBERSTATUS'] == 0, pspectra1.fibermap['FIBERSTATUS'] == 0)
        isGalaxy = zbest['SPECTYPE']=='GALAXY'
    
        hasSignal = HasSignal.filter(pspectra0,pspectra1)
        if norm:
            pspectra0, pspectra1 = renorm(pspectra0,pspectra1)
        diff = difference(pspectra0,pspectra1)
        
        skymask = maskskylines(diff)
        
        nspec = diff.flux[diff.bands[0]].shape[0]
        
        signal=np.zeros(nspec)
        signal[:]=np.zeros(nspec)
        var=np.zeros(nspec)
#         ref_signal=np.zeros(nspec)
#         w=np.where(diff.fibermap['TARGETID']==39627700813438689)[0]
        for dindex in diff.bands:
            wavecut = diff.wave[dindex]<HydrogenLogic.maxwave
            for sindex in range(nspec):
                z = zbest['Z'][sindex]
                if z>=zmin:
                    z_target_wave = (1+z)*HydrogenLogic.target_wave

                    #mask containing lines of interest
                    lmask = np.zeros(len(diff.wave[dindex]))
                    bmask = np.zeros(len(diff.wave[dindex]))
                    for wa in z_target_wave:
    #                     wmin = wa * np.exp(-1/HydrogenLogic.R/2.)
    #                     wmax = wa * np.exp(-1/HydrogenLogic.RH/2.)
    # #                     if sindex == w:
    # #                         print( wmin,wmax)
    #                     lmask = np.logical_or(lmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))
    #                     wmin = wa * np.exp(1/HydrogenLogic.RH/2.)
    #                     wmax = wa * np.exp(1/HydrogenLogic.R/2.)
    # #                     if sindex == w:
    # #                         print( wmin,wmax)
    #                     lmask = np.logical_or(lmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))
                        wmin = wa * np.exp(-1/HydrogenLogic.R/2.)
    #                     wmax = wa * np.exp(-1/HydrogenLogic.R/5.)
    #                     lmask = np.logical_or(lmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))
    #                     wmin = wa * np.exp(1/HydrogenLogic.R/5.)
                        wmax = wa * np.exp(1/HydrogenLogic.R/2.)
                        lmask = np.logical_or(lmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))

                        wmin = wa * np.exp(-3/HydrogenLogic.R)
                        wmax = wa * np.exp(-1/HydrogenLogic.R)
    #                     if sindex == w:
    #                         print( wmin,wmax)
                        bmask = np.logical_or(bmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))
                        wmin = wa * np.exp(1/HydrogenLogic.R)
                        wmax = wa * np.exp(3/HydrogenLogic.R)
    #                     if sindex == w:
    #                         print( wmin,wmax)
                        bmask = np.logical_or(bmask, np.logical_and.reduce((diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)))

                    #remove lines that are by the sky lines
                    lmask = np.logical_and.reduce((lmask, skymask[dindex] ==0,diff.mask[dindex][sindex,:]==0,wavecut))
                    bmask = np.logical_and.reduce((bmask, skymask[dindex] ==0,diff.mask[dindex][sindex,:]==0,wavecut))   

                    if (lmask.sum()>0 and bmask.sum()>0):
                        background = np.nanpercentile(diff.flux[dindex][sindex,bmask],50)
                        # only include unmasked
        #                 if sindex == w:
        #                     print(nmask.sum())
        #                     print( diff.flux[dindex][sindex,nmask])
    #                     print(background,lmask.sum(),diff.flux[dindex][sindex,lmask].sum())
                        signal[sindex] += diff.flux[dindex][sindex,lmask].sum() - lmask.sum()*background
                        var[sindex] += (1/diff.ivar[dindex][sindex,lmask]).sum()
        #                 ref_signal[sindex] += pspectra1.flux[dindex][sindex,nmask].sum()
    #             print(signal[sindex],var[sindex])
                else:
                    signal[sindex]=float('NaN')  # give a low redshift a bad signal
        ok = np.logical_and.reduce((isTGT, okFibers, isGalaxy))
        signal = signal - np.nanpercentile(signal[ok],50)
        onesig = np.nanpercentile(signal[ok],(15.865,84.135))

#         print(signal[w],signal[w]/np.sqrt(var[w]),5*onesig)
        large = np.logical_or(signal > 7*onesig[1], signal < 7* onesig[0])
        significant = (np.abs(signal)/ma.sqrt(var) >= ston_cut)
        triggered = np.logical_and.reduce((significant, isTGT, hasSignal, okFibers, isGalaxy, large))
        log.info("Number with signal {}; Number significant {}; Number large {}".format(hasSignal.sum(), significant.sum(), large.sum()))
        log.info("Number triggered {}".format(triggered.sum()))
        return triggered, diff

    # Color

class ColorLogic:
    plotter = plot_utils.diffplot_CV
    
    @staticmethod
    def filter(pspectra0, pspectra1, zbest, norm=False, ston_cut=7., frac_inc_cut= 1.):
        nspec = pspectra0.num_spectra()

        if not ('b' in pspectra1.bands and 'z' in pspectra1.bands and 'b' in pspectra0.bands and 'z' in pspectra0.bands):
            return np.zeros(nspec), None
        
        fibermap = pspectra0.fibermap #Table.read(datafile0, 'FIBERMAP')
        isTGT = fibermap['OBJTYPE'] == 'TGT'
        okFibers = np.logical_and(pspectra0.fibermap['FIBERSTATUS'] == 0, pspectra1.fibermap['FIBERSTATUS'] == 0)
        isGalaxy = zbest['SPECTYPE']=='GALAXY'
        
        hasSignal = HasSignal.filter(pspectra0,pspectra1,anyband=False)
                
        if norm:
            pspectra0, pspectra1 = renorm(pspectra0,pspectra1)
        
#         diff = difference(pspectra0,pspectra1)
#         w = np.where(fibermap['TARGETID']==39627712863667987)[0]
#         if len(w) !=0:
#             w=w[0]
#         else:
#             w=None

        signal=np.zeros(nspec)
        var=np.zeros(nspec)
        bignumber = 1e10
        for sindex in range(nspec):        
            nmask = np.logical_and(pspectra0.mask['b'][sindex,:]==0, pspectra1.mask['b'][sindex,:]==0)
            b0 = pspectra0.flux['b'][sindex,nmask].sum(); b1 = pspectra1.flux['b'][sindex,nmask].sum()
            b0var = (1/pspectra0.ivar['b'][sindex,nmask]).sum(); b1var = (1/pspectra1.ivar['b'][sindex,nmask]).sum()
            nmask = np.logical_and(pspectra0.mask['z'][sindex,:]==0, pspectra1.mask['z'][sindex,:]==0)
            z0 = pspectra0.flux['z'][sindex,nmask].sum(); z1 = pspectra1.flux['z'][sindex,nmask].sum()
            z0var = (1/pspectra0.ivar['z'][sindex,nmask]).sum(); z1var = (1/pspectra1.ivar['z'][sindex,nmask]).sum()
#             if sindex==w:
#                 log.info("{} {} {} {} {} {} {} {}".format(b0,np.sqrt(b0var),b1,np.sqrt(b1var),z0,np.sqrt(z0var),z1, np.sqrt(z1var)))
            
            if b0 <=0:
                logb0=-bignumber
            else:
                logb0= np.log(b0)
            if b1 <=0:
                logb1=-bignumber
            else:
                logb1= np.log(b1)                
            if z0 <=0:
                logz0=-bignumber
            else:
                logz0= np.log(z0)
            if z1 <=0:
                logz1=-bignumber
            else:
                logz1= np.log(z1)
                
            signal[sindex]  = logb0-logb1 - logz0 + logz1
            var[sindex] = b0var/b0**2 + b1var/b1**2 + z0var/z0**2 + z1var/z1**2
            
        # there is systematic calibration offset
        ok = np.logical_and.reduce((isTGT, okFibers,np.abs(signal)<bignumber/2))
        signal = signal - np.nanpercentile(signal[ok],50)
        onesig = np.nanpercentile(signal[ok],(15.865,84.135))
#         print(onesig)
#         fin=np.isfinite(signal)
#         fin=np.where(fin)[0]
#         print(np.nanpercentile(signal[ok],50))
#         print(signal[fin].mean(),signal[fin].std())
        
#         if w is not None:        
#             print(signal[w], var[w],hasSignal[w])
#             wef

        significant = (np.abs(signal)/ma.sqrt(var) >= ston_cut)
        large = np.logical_or(signal > 5*onesig[1], signal < 5* onesig[0])
        triggered = np.logical_and.reduce((significant, isTGT, hasSignal, isGalaxy, okFibers, large))
        log.info("Number with signal {}; Number significant {}; Number large {}".format(hasSignal.sum(), significant.sum(), large.sum()))
        log.info("Number triggered {}".format(triggered.sum()))
        return triggered, None
    
# single element logic
class SingleElementLogic:
    ston_cut=7.
    
    @staticmethod
    def filter(diff):    
        for i, dindex in enumerate(diff.bands):
            if i==0:
                ans = np.any(np.abs(diff.flux[dindex])*np.sqrt(diff.ivar[dindex]) > SingleElementLogic.ston_cut, axis=1)
            else:
                ans = np.logical_or(ans,np.any(np.abs(diff.flux[dindex])*np.sqrt(diff.ivar[dindex]) > SingleElementLogic.ston_cut, axis=1))
        return ans