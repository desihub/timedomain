import numpy

def coadd(newSpectra):
    sflux,  sivar, swave, smask = dict(), dict(), dict(), dict()
        
    first = True
    for s in newSpectra: # why/how would there be multiple spectra in s? 
        if first: 
            for b in s.bands:
                sflux[b]=s.flux[b]*s.ivar[b]
                sivar[b]=s.ivar[b]
                smask[b]=s.mask[b]
                swave[b]=s.wave[b]
            first=False
        else:
            for b in s.bands:
                sflux[b] += s.flux[b]*s.ivar[b]
                sivar[b] += s.ivar[b]
                smask[b] = smask[b]+s.mask[b]
    for b in s.bands:
        sflux[b] = sflux[b]/sivar[b] 
    
    return(sflux, sivar, swave, smask)
