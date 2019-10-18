import gzip
import numpy as np
from glob import glob
from astropy import units as u


class K17Sim:
    
    def __init__(self, filename):
        """Generate a simulation record.
        
        Parameters
        ----------
        filename : str
            Name of gzipped ASCII file with K17-KNe simulations.
        """
        self.time = None
        self.wave = None
        self.flux = None
        
        m, v, x = self.get_model_pars(filename)
        self.m_ej = m
        self.v_ej = v
        self.xlan = x
        
        for i, (day, wave, flux) in enumerate(self.get_records(filename)):
            flux = np.asarray(flux)
            if self.time is None:
                self.time = [day]
                self.wave = wave
                self.flux = [flux]
            else:
                self.time.append(day)
                self.flux.append(flux)
        
        self.time = np.asarray(self.time) * u.day
        self.wave = np.asarray(self.wave) * u.angstrom
        
    def get_model_pars(self, filename):
        """Extract K17 simulation model parameters from filename.

        Parameters
        ----------
        filename : str
            Name of a gzipped ASCII file.

        Returns
        -------
        meject : float
            Ejected mass.
        veject : float
            Ejection velocity.
        xlanth : float
            Lanthanide fraction.
        """
        tokens = filename.strip().split('_')
        meject = float(tokens[3][1:])
        veject = float(tokens[4][2:])
        xlanth = float(tokens[6][4:-9])
        
        return meject, veject, xlanth
        
    def get_records(self, filename):
        """Read one of the K17 kilnova models and extract one time record.

        Parameters
        ----------
        filename : str
            Name of a gzipped ASCII file.

        Returns
        -------
        day : float
            Date of explosion w.r.t. t0.
        wave : list
            List of wavelength values, in Angstroms.
        flux : list
            List of fluxes, in ?
        """
        day, wave, flux = None, [], []

        with gzip.open(filename, 'r') as f:
            for line in f:
                if line.startswith(b'#'):
                    continue

                t, wl, fl = [float(x) for x in line.strip().split()]
                if t != day:
                    if len(wave) > 0:
                        yield day, wave, flux
                        wave, flux = [], []
                    day = t

                wave.append(wl)
                flux.append(fl)

            if len(wave) > 0:
                yield day, wave, flux
                
    def get_flux(self, time):
        """Get the flux at (or near) a given time.
        
        Parameters
        ----------
        time : float
            Time of the flux sample.
            
        Returns
        -------
        time : float
            Simulated time nearest the user-specified time.
        flux : ndarray
            Flux at the simulated time nearest the user-specified time.
        """
        k = (np.abs(self.time-time)).argmin()
        return self.time[k], self.flux[k]


class K17Model:
    
    def __init__(self, a, b, m_ej, v_ej, xlan, time, wave, flux):
        """Initialize a kilonova model as the sum of two kilonova simulations.
        
        Parameters
        ----------
        a : float
            Weight (flux scale factor) of the first kilonova simulation.
        b : float
            Weight (flux scale factor) of the second kilonova simulation.
        m_ej : list
            Ejecta masses of the first and second simulations.
        v_ej : list
            Ejecta velocities of the first and second simulations.
        xlan : list
            Lanthanide fractions of the first and second simulations.
        time : ndarray
            List of available times w.r.t. t0.
        wave : ndarray
            List of wavelengths simulated.
        flux : list
            2D list of summed fluxes as a function of time.
        """
        self.a = a
        self.b = b
        self.m_ej = m_ej
        self.v_ej = v_ej
        self.xlan = xlan
        self.time = time
        self.wave = wave
        self.flux = flux
        
    def get_flux(self, time):
        """Get the flux at (or near) a given time.
        
        Parameters
        ----------
        time : float
            Time of the flux sample.
            
        Returns
        -------
        time : float
            Simulated time nearest the user-specified time.
        flux : ndarray
            Flux at the simulated time nearest the user-specified time.
        """
        k = (np.abs(self.time-time)).argmin()
        return self.time[k], self.flux[k]


class K17Generator:
    """Randomly pick two kilonova simulations and add the fluxes, picking
    a random scale factor between 0 and 1 for each flux."""
    
    def __init__(self, prefix=None):
        self.m_ej = []      # Ejecta mass (proxy for luminosity).
        self.v_ej = []      # Ejecta velocity.
        self.xlan = []      # Lanthanide fraction.
        self.idx_lo = []    # Low-lanthanide simulations (file index).
        self.idx_hi = []    # High-lanthanide simulations (file index).
        
        if prefix is None:
            self.ksfiles = sorted(glob('knova_d1_n10_m0.*txt.gz'))
        else:
            self.ksfiles = sorted(glob(prefix + '/knova_d1_n10_m0.*txt.gz'))
            
        for idx, ksfile in enumerate(self.ksfiles):
            m, v, x = self.get_model_pars(ksfile)
            self.m_ej.append(m)
            self.v_ej.append(v)
            self.xlan.append(x)
            if np.log10(x) < -2:
                self.idx_lo.append(idx)
            else:
                self.idx_hi.append(idx)
        
    def get_model_pars(self, filename):
        """Extract K17 simulation model parameters from filename.

        Parameters
        ----------
        filename : str
            Name of a gzipped ASCII file.

        Returns
        -------
        meject : float
            Ejected mass.
        veject : float
            Ejection velocity.
        xlanth : float
            Lanthanide fraction.
        """
        tokens = filename.strip().split('_')
        meject = float(tokens[3][1:])
        veject = float(tokens[4][2:])
        xlanth = float(tokens[6][4:-9])
        
        return meject, veject, xlanth
    
    def generate_model(self):
        # Grab two models (get indices without replacement.)
        i = np.random.choice(self.idx_lo)
        j = np.random.choice(self.idx_hi)
        #a, b = np.random.uniform(low=0., high=1., size=2)
        a = np.random.uniform()
        b = 1. - a
        
        ks_i = K17Sim(self.ksfiles[i])
        ks_j = K17Sim(self.ksfiles[j])
        
        flux = []
        for k in range(len(ks_i.time)):
            fl = a*ks_i.flux[k] + b*ks_j.flux[k]
            flux.append(fl)
        
        return K17Model(a, b, [ks_i.m_ej, ks_j.m_ej],
                              [ks_i.v_ej, ks_j.v_ej],
                              [ks_i.xlan, ks_j.xlan],
                              ks_i.time, ks_i.wave, flux)
