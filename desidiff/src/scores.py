import numpy

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

def perconv_SN(y,ivar,mask,ncon=3,nsig=10):
    newy=dict(y)
    newivar=dict(ivar)
    newmask=dict(mask)
    ncon = numpy.zeros(ncon)+1.
    for b in newy.keys():
        newivar=numpy.convolve(ivar[b],ncon,mode='valid')
        newy = numpy.convolve(y[b]*ivar[b],ncon,mode='valid')
        newy = newy/newivar
        newmask=numpy.convolve(mask[b],ncon,mode='valid')
        
    return perres_SN(y,ivar,mask,nsig=nsig)

def Hlines(wave, y,ivar,mask, z):
    target_wave = (6562.79, 4861.35, 4340.472, 4101.734, 3970.075)
    R=1000.
    
    target_wave = numpy.array(target_wave)*(1+z)
    
    signal=0.
    var=0.
    
    for dindex in diff.bands:
        #mask containing lines of interest
        lmask = numpy.zeros(len(wave[dindex]))

        for wa in target_wave:
            wmin = wa * np.exp(-1/R/2.)
            wmax = wa * np.exp(1/R/2.)
            lmask = numpy.logical_or(lmask, np.logical_and.reduce((wave[dindex] >= wmin, wave[dindex] < wmax)))

        # only include unmasked
        lmask = np.logical_and(lmask, mask[dindex]==0)
        signal += y[dindex][lmask].sum()
        var += (1/ivar[dindex][lmask]).sum()

    ston = numpy.abs(signal)/ma.sqrt(var)