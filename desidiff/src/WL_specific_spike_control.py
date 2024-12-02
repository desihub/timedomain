### implemention looks like:
### difmask = WL_specific_spike_control(difflux, difwave, difmask, 10, 10, 3)
import copy

#takes in a mask, spectra (i.e. wavelength, flux), and an applicable threshold for spikes in each band, only a mask is returned
def WL_specific_spike_control(flux, wave, mask, b_threshold, r_threshold, z_threshold):
    # deepcopy of flux and wave does not affect spectra data
    yc = copy.deepcopy(flux)
    newmask = dict(mask)
    xc = copy.deepcopy(wave)
    
    # description of applied wavelengths 
    b_window = [2464, 2485] #5571.2 A -> 5588 A
    r_window = [2301, 2322] #7600.8 A -> 7616.8 A
    z_window = [0, 1974] #7520 A -> 9099.2 A

    for b in yc.keys():
        if b == 'b':
            yb = yc[b][2464:2485]
            dely = abs(numpy.diff(yb))
            for i in range(len(dely)):
                if dely[i] > b_threshold:
                    newmask[b][i+2464] = 1
                    
        if b == 'r':
            yb = yc[b][2301:2322]
            dely = abs(numpy.diff(yb))
            for i in range(len(dely)):
                if dely[i] > r_threshold:
                    newmask[b][i+2301] = 1

        if b == 'z':
            yb = yc[b][:]
            dely = abs(numpy.diff(yb))
            for i in range(len(dely)):
                if dely[i] > z_threshold:
                    newmask[b][i+1] = 1

        return(newmask)
