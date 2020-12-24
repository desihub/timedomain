import numpy as np
import numpy.ma as ma

from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

from . import plot_utils

# difference between two spectra
# for now the two spectra are assumed to be single night coadd of matching tile

def difference(s0, s1):
    diff = dict()
    ivar = dict()
    mask = dict()
    for dindex in s1.bands:

        diff[dindex] = ma.array(data=s1.flux[dindex],mask=s1.mask[dindex])
        diff[dindex] = diff[dindex] - ma.array(data=s0.flux[dindex],mask=s0.mask[dindex])

        ivar0 = ma.array(data=s0.ivar[dindex],mask=s0.mask[dindex])
        ivar1 = ma.array(data=s1.ivar[dindex],mask=s1.mask[dindex])
        ivar[dindex] = 1/(1/ivar0 + 1/ivar1)
        mask[dindex] = diff[dindex].mask.astype('int')
        
    return Spectra(bands=s1.bands.copy(), wave=dict(s1.wave), flux=diff,ivar=ivar,mask=mask,fibermap = s0.fibermap)

# renormalize two spectra
# for now the two spectra are assumed to be single night coadd of matching tile

def renorm(s0, s1):

    for dindex in s1.bands:
        norm = ma.array(data=s1.flux[dindex],mask=s1.mask[dindex])/ ma.array(data=s0.flux[dindex],mask=s0.mask[dindex])
        norm = np.percentile(norm,(50),axis=1)
        s0.flux[dindex] = s0.flux[dindex]*norm[:,None]
        s0.ivar[dindex] = s0.ivar[dindex]/((norm*norm)[:,None])

    return s0,s1


# Sometimes the spectrum has very low signal, perhaps due to mispointed fiber

class HasSignal:
    ston_cut=20.
    
    @staticmethod
    def filter(s0,s1):    
        for i, dindex in enumerate(s0.bands):
            ston1 = s1.flux[dindex].sum(axis=1) / np.sqrt((1/s1.ivar[dindex]).sum(axis=1)) 
            ston0 = s0.flux[dindex].sum(axis=1) / np.sqrt((1/s0.ivar[dindex]).sum(axis=1)) 

            if i==0:
                ans = np.logical_and(np.greater(ston0,HasSignal.ston_cut), np.greater(ston1,HasSignal.ston_cut))
            else:
                ans = np.logical_or(ans, np.logical_and(np.greater(ston0,HasSignal.ston_cut), np.greater(ston1,HasSignal.ston_cut)))
        return ans
    

# Cataclysmic Variable

class CVLogic:
    target_wave = (6562.79, 4861.35, 4340.472, 4101.734, 3970.075)
    R=200.
    ston_cut=5.
    plotter = plot_utils.diffplot_CV
    @staticmethod
    def filter(pspectra0, pspectra1, renorm=False):
        
        
        
        fibermap = pspectra0.fibermap #Table.read(datafile0, 'FIBERMAP')
        isTGT = fibermap['OBJTYPE'] == 'TGT'
        # the interesting guy is wedge 6 index 331
    #     print(fibermap['TARGETID'].data[0], np.where(fibermap['TARGETID'].data== 35191288252861933)[0])
    #     pspectra0 = read_spectra(datafile0)
    #     pspectra1 = read_spectra(datafile1)
        hasSignal = HasSignal.filter(pspectra0,pspectra1)

        if renorm:
            pspectra0, pspectra1 = renorm(pspectra0,pspectra1)

        diff = difference(pspectra0,pspectra1)
        
        signal=0.
        var=0.
        for dindex in diff.bands:
            for wa in CVLogic.target_wave:
                wmin = wa * np.exp(-1/CVLogic.R)
                wmax = wa * np.exp(1/CVLogic.R)
                w = np.logical_and(diff.wave[dindex] >= wmin, diff.wave[dindex] < wmax)
                signal +=  diff.flux[dindex][:,w].sum(axis=1)
                var += (diff.ivar[dindex][:,w]).sum(axis=1)
                
        significant = (np.abs(signal)/np.sqrt(var) >= CVLogic.ston_cut)
        triggered = np.logical_and.reduce((significant, isTGT, hasSignal))
        return triggered, diff
    
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