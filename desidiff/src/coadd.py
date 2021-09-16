import numpy

def coadd(newSpectra):
    sflux,  sivar, swave, smask = dict(), dict(), dict(), dict()
        
    first = True
    for s in newSpectra: # why/how would there be multiple spectra in s? 
        for i in range(s.num_spectra()):
            if first: 
                for b in s.bands:
                    sflux[b]=s.flux[b][i,:]*s.ivar[b][i,:]
                    sivar[b]=s.ivar[b][i,:]
                    smask[b]=numpy.abs(s.mask[b][i,:])
                    swave[b]=s.wave[b][:]
                first=False
            else:
                for b in s.bands:
                    sflux[b] += s.flux[b][i,:]*s.ivar[b][i,:]
                    sivar[b] += s.ivar[b][i,:]
                    smask[b] += numpy.abs(s.mask[b][i,:])
                    
    for b in s.bands:
        # Mask pixels with no signal-to-noise
        dum = sivar[b]==0
        w=numpy.where(dum)[0]
        smask[b][w]=1
        w=numpy.where(~dum)[0]    
        sflux[b][w] = sflux[b][w]/sivar[b][w]
    
    return(sflux, sivar, swave, smask)