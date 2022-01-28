import numpy
import numpy.ma as ma
import copy
def clipmean_one(y,ivar,mask,nsig=3):

    w=numpy.where(mask==0)[0]
    ansivar = ivar[w].sum()
    
    ansmean = numpy.sum(y[w]*ivar[w])/ansivar
    newy = y-ansmean
    w=numpy.where(numpy.logical_and.reduce([mask==0, numpy.abs(newy*numpy.sqrt(ivar)) < nsig]))[0]
    ansivar = ivar[w].sum()
    ansmean = numpy.sum(y[w]*ivar[w])/ansivar
    return y-ansmean

def clipmean(y,ivar,mask,nsig=3):
    ans = copy.deepcopy(y)
    for k in ans.keys():
        ans[k]=clipmean_one(y[k],ivar[k],mask[k],nsig=nsig)
    return ans

def perband_SN(y,ivar,mask,nsig=10):
    ans=dict()
    for k in y.keys():
        w=numpy.where(mask[k]==0)[0]
        ansivar = ivar[k][w].sum()
        ansmean = numpy.sum(y[k][w]*ivar[k][w])/ansivar
        ston=numpy.abs(ansmean)*numpy.sqrt(ansivar)
        ans[k]=ston
    return ans

def perband_increase(y,ivar,mask,refy, refivar, refmask):
    ans=dict()
    for k in y.keys():
        w=numpy.where(mask[k]==0)[0]
        ansivar = ivar[k][w].sum()
        ansmean = numpy.sum(y[k][w]*ivar[k][w])/ansivar
        ans[k]=ansmean
        w=numpy.where(refmask[k]==0)[0]
        ansivar = refivar[k][w].sum()
        ansmean = numpy.sum(refy[k][w]*refivar[k][w])/ansivar
        ans[k]=ans[k]/ansmean
    return ans

def perres_SN(y,ivar,mask,nsig=10):
    ans=dict()
    # use observed dispersion rather than statistical
    for k in y.keys():
        w=numpy.where(mask[k]==0)[0]
        std=y[k][w].std()
        ston=numpy.abs(y[k][w])/std
#         ston=numpy.abs(y[k][w])*numpy.sqrt(ivar[k][w])
        ans[k]=(ston>10).sum()
    return ans

def perconv_SN(wave, y, ivar, mask, ncon=5, nsig=10):
    newmask = dict(mask)
    newston = dict(y)
    newy = dict(y)
    conk = numpy.zeros(ncon) + 1.

    nTrue = 0
    for b in newston.keys():
        # mask lines that are smaller than the ncov

        f_logic = y[b]*numpy.sqrt(ivar[b]) > nsig/numpy.sqrt(ncon)
        index=0
        for _, g in groupby(f_logic):
            ng  = sum(1 for _ in g)
            if _ and ng<ncon//2:
                newy[b][index:index+ng] = 0
            index=index+ng

        newston[b]=numpy.convolve(newy[b], conk, mode='valid') * numpy.sqrt(numpy.convolve(1/ivar[b], conk, mode='valid'))
        newmask[b] = numpy.convolve(mask[b], conk, mode='valid')

        w = numpy.where(newmask[b] == 0)[0]
        for gb0,_ in groupby(newston[b][w] > nsig):
            if gb0:
                nTrue=nTrue+1
    return nTrue

def Hlines(wave, y,ivar,mask, z):
    target_wave = (6562.79, 4861.35, 4340.472, 4101.734, 3970.075)
    R=500.
    
    target_wave = numpy.array(target_wave)*(1+z)
    
    signal=0.
    var=0.
    
    for dindex in y.keys():        
        for wa in target_wave:
            wmin = wa * numpy.exp(-1/R/2.)
            wmax = wa * numpy.exp(1/R/2.)
            lmask = numpy.logical_and.reduce((wave[dindex] >= wmin, wave[dindex] < wmax))
            lmask = numpy.logical_and(lmask, mask[dindex]==0)
            if lmask.sum() >1:
                
                #Window for background
                wmin = wa * numpy.exp(-5/(R/5)/2.)
                wmax = wa * numpy.exp(-3/(R/5)/2.)
                bmask = numpy.logical_and.reduce((wave[dindex] >= wmin, wave[dindex] < wmax))

                wmin = wa * numpy.exp(3/(R/5)/2.)
                wmax = wa * numpy.exp(5/(R/5)/2.)
                bmask = numpy.logical_or(bmask, numpy.logical_and.reduce((wave[dindex] >= wmin, wave[dindex] < wmax)))
                bmask = numpy.logical_and(bmask, mask[dindex]==0)
                background=(y[dindex][bmask]*ivar[dindex][bmask]).sum()/ivar[dindex][bmask].sum()
 
                # only include unmasked
                signal += (y[dindex][lmask].sum()-background*lmask.sum())
                var += (1/ivar[dindex][lmask]).sum()

    ston = numpy.abs(signal)/numpy.sqrt(var)
    return ston


# renormalize two spectra
# for now the two spectra are assumed to be single night coadd of matching tile

def normalization(newflux, newmask, refflux, refmask):

    common = list(set(newflux.keys()).intersection(refflux.keys()))
    norms=[]
    for dindex in common:
        norm = ma.array(data=newflux[dindex],mask=newmask[dindex])/ ma.array(data=refflux[dindex],mask=refmask[dindex])
        norm.filled(numpy.nan)
        norms.append(norm)
        
    norms=numpy.concatenate(norms)  
    norm = numpy.nanpercentile(norms,(50))

    return norm