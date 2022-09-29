import healpy as hp
import numpy as np
from scipy.stats import norm

#Following https://iopscience.iop.org/article/10.3847/0067-0049/226/1/10/pdf
def dist_from_map(filename):
    prob, distmu, distsigma, distnorm =hp.read_map(filename,field=[0,1,2,3])
    npix = len(prob)
    nside = hp.npix2nside(npix)
    r = np.linspace(0,5000,1000)
    pix_of_maxprob = np.argmax(prob)
    dp_dr = r**2 * distnorm[pix_of_maxprob] * norm(distmu[pix_of_maxprob],\
                                              distsigma[pix_of_maxprob]).pdf(r)
    return np.average(r,weights=dp_dr)

#dist_from_map('../OUTPUT/O3REAL/S190408an/bayestar/bayestar.fits.gz')
