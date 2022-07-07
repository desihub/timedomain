import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
from scipy.ndimage import gaussian_filter1d

### Different Line Profiles
def continuum(x,m,b):
    y = m*x + b
    return y

def lorentzian(x,x0,gamma,A):
    y = (A*0.25*gamma**2)/((x - x0)**2 + (gamma/2)**2) 
    return y
def gaussian(x,mu,sigma,a):
    y = a*np.exp(-(x-mu)**2/(2*sigma)**2)
    return y

#defining the wavelength range to be imilar to desi, constructing spectrum, adding noise, plotting
wave = np.linspace(3000,10000,8000)
spectrum = continuum(wave, -0.0025,50) + gaussian(wave,6562,120,10)
spectrum = [rd.normal(i,3) for i in spectrum]


plt.plot(wave,spectrum)
plt.show()

spectrum = np.array(spectrum)

#stolen from a stackexchange page, just shifts the entire array by n spaces
def shift(xs, n):
    if n >= 0:
        return np.concatenate((np.full(n, np.nan), xs[:-n]))
    else:
        return np.concatenate((xs[-n:], np.full(-n, np.nan)))

#shift by 1, subtract, huge gaussian filter 
specplus = shift(spectrum,1)

firstderiv = specplus -spectrum
firstderiv = gaussian_filter1d(firstderiv, 30)
plt.plot(wave,firstderiv)
plt.show()

#do it again for teh second derivative
secondderiv = shift(firstderiv,1) - firstderiv
secondderiv = gaussian_filter1d(secondderiv, 30)
plt.plot(wave,secondderiv)
plt.show()