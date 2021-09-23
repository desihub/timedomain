                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                d by the continuum to remove its slope.
        """
        blue_min, blue_max = self.cont_blue
        blue_avg = 0.5*(blue_min + blue_max)
        select_blue = (roi_wave > blue_min) & (roi_wave < blue_max)
        flux_blue = roi_flux[select_blue]
        
        red_min, red_max = self.cont_red
        red_avg = 0.5*(red_min + red_max)
        select_red = (roi_wave > red_min) & (roi_wave < red_max)
        flux_red = roi_flux[select_red]
        
        # Perform a linear fit to the mean flux vs mean wavelength in the blue and red sidebands of the line.
        # a = intercept, b = slope.
        b = (np.nanmean(flux_red) - np.nanmean(flux_blue)) / (red_avg - blue_avg)
        a = np.nanmean(flux_blue) - (b * blue_avg)

        # Calculate the continuum using the linear fit.
        continuum = a + b*roi_wave
        
        # Scale the flux by the linear continuum.
        scaled_roi_flux = roi_flux / continuum
        return a, b, continuum, scaled_roi_flux
    
    def estimate_peqw(self, roi_wave, scaled_roi_flux):
        """Compute the pseudo-equivalent width of the spectral line.
        
        Parameters
        ----------
        roi_wave : ndarray
            Wavelength in the fit region of interest of the line.
        scaled_roi_flux : ndarray
            Flux in the ROI scaled by the continuum to remove its slope.
            
        Returns
        -------
        """
        continuum = np.ones_like(scaled_roi_flux)
        _sum = np.sum((scaled_roi_flux - continuum) / continuum)
        _delta = np.diff(roi_wave)[0]
        return _sum * _delta

    def get_linewidth(self, wave, flux):
        """Return the line width (pseudo-equivalent width).
        
        Parameters
        ----------
        wave : ndarray
            Input wavelength from a spectrum.
        flux : ndarray
            Input flux from a spectrum.
        
        Returns
        -------
        peqw : float
            Pseudo-equivalent width of the spectral line.
        """
        vs, roi_wave, roi_flux = self.get_roi(wave, flux,self.line - np.nanmean(self.cont_blue),\
                                              np.nanmean(self.cont_red)-self.line)
        a, b, continuum, scaled_roi_flux = self.estimate_continuum(roi_wave, roi_flux)
        peqw = self.estimate_peqw(roi_wave, scaled_roi_flux)
        
        # Store the details of the continuum fit for debugging.
        self.vshift = vs
        self.roi_wave = roi_wave
        self.roi_flux = roi_flux
        self.roi_flux_scaled = scaled_roi_flux
        self.roi_continuum = continuum
        
        return peqw

    
    def __str__(self):
        s = '{}: {}\nBlue continuum fit: {} - {}\nRed continuum fit:  {} - {}'.format(
            self.name, self.line,
            self.cont_blue[0], self.cont_blue[1],
            self.cont_red[0], self.cont_red[1])
        return s
    
def is_TDE(Ha,HB,HeII,OIII,NIII, flux):
    '''
    Ha,HB, HeII, and OIII should be pEW for each of those lines respectively.
    '''
    TDE_score = 0 
    blue_end = np.nanmean(flux[1000:3000])
    red_end = np.nanmean(flux[5000:7000])
    mean_flux = np.nanmean(flux)
    try:
        if OIII/Ha > 0.1:
            TDE_score -= 1.1 #if OIII is present, will be marked by an index with a noninteger part
    except ZeroDivisionError:
        pass
    if Ha > 4: #Just looks for generally wide Ha lines
        TDE_score += 1
    if HB > 4:
        TDE_score += 1
    try:
        if HeII/Ha >0.75: #looks for comparable He and Ha lines
            TDE_score += 1
    except ZeroDivisionError:
        pass
    if NIII > 4: #looks for NIII sometimes seen in TDE-bowen 
        TDE_score += 1
    if blue_end/red_end >= 1.25:
        TDE_score += 1
    if mean_flux < 1:
        TDE_score -= 1.1
    
    return TDE_score
def TDE_Check(wave,flux, redshift, scale = 3.,cont_width = 10):
    '''
    Creates a scaled flux (assumed to be roughly gaussian) for each line of interest. Then fits this scaled flux to a
    Gaussian, pulling the $\sigma$. Recalculates pEW using an ROI of the line $\pm$ scale*sigma. Appends this pEW to
    a Dictionary (lines), each line should have a list in format [wavelength, pEW]. pEW used in logic to determine whether 
    or not it's a TDE.
    Parameters:
    wave, flux: wavelength and flux of spectrum
    scale: # of standard deviations to fit continuum
    cont_width: width of region to sample cont_blue and cont_red, default 5 grabs area 2.5 angstroms on either side of 
        bound determined by scale*sigma
    '''
    wave = wave/(1+redshift)
    flux = gaussian_filter1d(flux, sigma = 3)
    half_cont = cont_width/2
    lines = {'Halpha':[6562.79],'Hbeta':[4861.4], 'HeII4686':[4686],'OIII':[5007],'NIII':[4100]}
    for line in lines:
        pass1 = SpectralLine_pEW(line,lines[line][0])
        try:
            peqw_pass1 = pass1.get_linewidth(wave, flux)
            #plt.plot(pass1.roi_wave, pass1.roi_flux)
            #plt.show()
            #print(line)
            #print(peqw_pass1)

            try:
                popt = curve_fit(gaus, pass1.roi_wave, pass1.roi_flux_scaled, p0 = [10.,pass1.line,max(pass1.roi_flux_scaled),1],check_finite = False)[0]
                std = popt[0]
                #print(std)
                if scale*std <= half_cont: #Fit failed, line probably not present enough for measurement. 
                    peqw_specline = 0
                    #print('Failed')

                else:
                    cont_blue=[pass1.line - scale*std-half_cont, pass1.line - scale*std+half_cont]
                    cont_red=[pass1.line + scale*std-half_cont, pass1.line + scale*std+half_cont]
                    specline = SpectralLine_pEW(line, lines[line][0], cont_blue=cont_blue, \
                                cont_red=cont_red)
                    try:
                        peqw_specline = specline.get_linewidth(wave, flux)
                        #print(cont_blue,cont_red)
                        #print(peqw_specline)
                        #plt.plot(specline.roi_wave, specline.roi_flux)
                        #plt.show()
                    except IndexError:
                        peqw_specline = 0
            #print('###########')
            except RuntimeError:
                peqw_specline = 0
        except:
            peqw_specline = 0
        lines[line].append(peqw_specline)
    return is_TDE(lines['Halpha'][-1],lines['Hbeta'][-1],lines['HeII4686'][-1],lines['OIII'][-1],\
                  lines['NIII'][-1], flux)
