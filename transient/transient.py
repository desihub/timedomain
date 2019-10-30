"""Module for defining interface to Transient simulators, plus subclasses.
"""

from abc import ABC, abstractmethod
import numpy as np
import sncosmo
from astropy import units as u


class Transient(ABC):

    def __init__(self, ratio=1., phase=0*u.day):
        """Initialize transient.

        Parameters
        ----------
        ratio : float
            Transient to host flux ratio.
        phase : float
            Phase/epoch of the transient.
        """
        # Ratio of transient to host flux (default: r-band).
        self.ratio = ratio

        # Phase/epoch within the transient.
        self.phase = phase

    @abstractmethod
    def sed(self, date):
        """Return the spectral energy distribution of the transient.
        """
        pass


class Supernova(Transient):
    
    def __init__(self, ratio=1., phase=0*u.day, modeltype=None, modelname=None):
        """A supernova model based on the sncosmo package.

        Parameters
        ----------
        ratio : float
            Transient to host flux ratio.
        phase : float
            Phase/epoch of the transient.
        modeltype : str or None
            Supernova type (SN Ia, SN IIP, ...).
        modelname : str or None
            Supernova sncosmo model name.
        """
        super().__init__(ratio, phase)

        # Built-in sncosmo supernova models.
        self.models = {
            'SN Ia'     : ['nugent-sn1a', 'nugent-sn91t', 'nugent-sn91bg',
                           'salt2-h17',   'hsiao',        'salt2-extended-h17',
                           'snf-2011fe'],
            
            'SN Ib'     : ['s11-2005hl',   's11-2005hm',   's11-2006jo',
                           'snana-2004gv', 'snana-2006ep', 'snana-2007y',
                           'snana-2004ib', 'snana-2005hm', 'snana-2006jo',
                           'snana-2007nc'],

            'SN Ib/c'   : ['nugent-sn1bc'],

            'SN Ic'     : ['s11-2006fo',   'snana-2004fe',     'snana-sdss004012',
                           'snana-2006fo', 'snana-sdss014475', 'snana-2006lc',
                           'snana-04d1la', 'snana-04d4jv'],

            'SN IIn'    : ['nugent-sn2n', 'snana-2006ez', 'snana-2006ix'],

            'SN IIP'    : ['nugent-sn2p',  's11-2005lc',   's11-2005gi',
                           's11-2006jl',   'snana-2004hx', 'snana-2005gi',
                           'snana-2006gq', 'snana-2006kn', 'snana-2006jl',
                           'snana-2006iw', 'snana-2006kv', 'snana-2006ns',
                           'snana-2007iz', 'snana-2007nr', 'snana-2007kw',
                           'snana-2007ky', 'snana-2007lj', 'snana-2007lb',
                           'snana-2007ll', 'snana-2007nw', 'snana-2007ld',
                           'snana-2007md', 'snana-2007lz', 'snana-2007lx',
                           'snana-2007og', 'snana-2007nv', 'snana-2007pg'],
            
            'SN IIL'    : ['nugent-sn2l'],

            'SN IIL/P'  : ['s11-2004hx'],

            'SN II-pec' : ['snana-2007ms']
        }

        self.modelname = None
        self.modeltype = None

        # Either use an existing model or randomly select one. Add necessary
        # checks for invalid type names and model names.
        if modeltype is None:
            if modelname is None:
                self.modeltype = np.random.choice(self.models.keys())
                self.modelname = np.random.choice(self.models[self.modeltype])
            else:
                for t, n in self.models.items():
                    if modelname in n:
                        self.modelname = modelname
                        self.modeltype = t
        else:
            if modeltype in self.models.keys():
                self.modeltype = modeltype
                if modelname is None:
                    self.modelname = np.random.choice(self.models[modeltype])
                else:
                    if modelname in self.models[self.modeltype]:
                        self.modelname = modelname
        
        if self.modelname is None or self.modeltype is None:
            raise ValueError('Invalid type {} or model name {}'.format(modelname, modeltype))

        # Initialize the model.
        self.model = sncosmo.Model(self.modelname)
        self.model.set(z=0., t0=0.)
        self.wave = np.arange(np.ceil(self.model.minwave()),
                              np.floor(self.model.maxwave()), 1.) * u.Angstrom
        self.time = np.arange(np.ceil(self.model.mintime()),
                              np.floor(self.model.maxtime()), 1.) * u.day

    def sed(self, date=None):
        """Return the spectral energy distribution of the supernova.
        """
        if date is None:
            date = np.random.choice(self.time)
        return self.model.flux(date, self.wave)

    def get_types(self):
        """Get a list of available SN types.
        """
        return list(self.models.keys())

    def get_models(self):
        """Get a list of available SN model names.
        """
        models = []
        for m in self.models.values():
            models += m
        return models

