from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from . import filters

# savedir  = '/global/cscratch1/sd/akim/project/timedomain/output/'
savedir = '/global/cfs/cdirs/desi/science/td/daily-search/'
def diffplot_CV(sig,pspectra0, pspectra1, diff,savepdf=False):

    lims = []
    fig, axes = plt.subplots(3,1, figsize=(15,10), sharex=True, sharey=False, gridspec_kw={'wspace':0, 'hspace':0})
    for dindex in pspectra0.bands:
        wave = pspectra0.wave[dindex]
#             diff = ma.array(data=pspectra1.flux[dindex],mask=pspectra1.mask[dindex])
#             diff = diff - ma.array(data=pspectra0.flux[dindex],mask=pspectra0.mask[dindex])                        
        axes[0].plot(wave,pspectra0.flux[dindex][sig,:],alpha=0.7)
        axes[1].plot(wave,pspectra1.flux[dindex][sig,:],alpha=0.7)
        axes[2].plot(wave,diff.flux[dindex][sig,:],alpha=0.7)
#             axes[2].plot(wave,smooth(diff[sig,:]),color='red',alpha=0.5)
#                        axes[2].plot(wave,diff[whereSky[0][0],:],color='yellow',alpha=0.5)
        lims.append(np.percentile(diff.flux[dindex][sig,:],(0.01,99.99)))
    lims=np.array(lims)
    lims = [lims.min()*1.05,lims.max()*1.05]
    axes[2].set_ylim(lims)
    for wav in filters.CVLogic.target_wave:
        axes[2].axvline(wav,ls=':',color='black')
    axes[0].set_title("{} RA: {} DEC: {}".format(diff.fibermap['TARGETID'].data[sig],diff.fibermap['TARGET_RA'].data[sig],diff.fibermap['TARGET_DEC'].data[sig]))
    fig.tight_layout()
    fig.show()
    
    if savepdf:
        tile = pspectra0.fibermap["TILEID"].data[sig]
        expid0 = pspectra0.fibermap["EXPID"].data[sig]
        expid1 = pspectra1.fibermap["EXPID"].data[sig]
        outdir = os.path.join(savedir,*savepdf)
        Path(outdir).mkdir(parents=True, exist_ok=True)        

        pdf = PdfPages(os.path.join(outdir,'{}_{}_{}.pdf'.format(tile,expid0,expid1)))
        pdf.savefig()
        pdf.close()
        
        
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError( "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[int(window_len/2-1):-int(window_len/2)-1]