import numpy

def clipmean(y,ivar,mask,nsig=3):
    w=numpy.where(mask==0)[0]
    ansivar = ivar[w].sum()
    ansmean = numpy.sum(y[w]*ivar[w])/ansivar
    
    newy = y-ansmean
    w=numpy.where(numpy.logical_and(mask==0,numpy.abs(newy*sqrt(ivar)) < nsig))[0]
    
    ansivar = ivar[w].sum()
    ansmean = numpy.sum(y[w]*ivar[w])/ansivar
    
    return y-ansmean

