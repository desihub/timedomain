"""
Docs » Module code » desispec.coaddition
Source code for desispec.coaddition

Coadd spectra
"""

from __future__ import absolute_import, division, print_function
import os, sys, time

import numpy as np

import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg

from astropy.table import Table, Column

import multiprocessing

from desiutil.log import get_logger

from desispec.interpolation import resample_flux
from desispec.spectra import Spectra
from desispec.resolution import Resolution
from desispec.fiberbitmasking import get_all_fiberbitmask_with_amp, get_all_nonamp_fiberbitmask_val, get_justamps_fiberbitmask
from desispec.specscore import compute_coadd_scores

class spectra_no_expid(Spectra):
    def show(self):
        self.fibermap
        print(type(self))
        return
    

    
    
    def coadd_no_expid(spectra_no_expid, cosmics_nsig=0.) :
        """
        Coaddition the spectra for each target and each camera. The input spectra is modified.

        Args:
           spectra: desispec.spectra.Spectra object

        Options:
           cosmics_nsig: float, nsigma clipping threshold for cosmics rays
        """

        spectra_no_expid.fibermap=coadd_fibermap_no_expid(spectra_no_expid.fibermap)
        return

    def select_no_expid(self, nights=None, exposures=None, bands=None, targets=None, fibers=None, invert=False, return_index=False):
        """
        Select a subset of the data.

        This filters the data based on a logical AND of the different
        criteria, optionally inverting that selection.

        Args:
            nights (list): optional list of nights to select.
            exposures (list): optional list of exposures to select.
            bands (list): optional list of bands to select.
            targets (list): optional list of target IDs to select.
            fibers (list): list/array of fiber indices to select.
            invert (bool): after combining all criteria, invert selection.
            return_index (bool): if True, also return the indices of selected spectra.

        Returns:
            spectra: a new Spectra object containing the selected data.
            indices (list, optional): indices of selected spectra. Only provided if return_index is True.
        """
        if bands is None:
            keep_bands = self.bands
        else:
            keep_bands = [ x for x in self.bands if x in bands ]
        if len(keep_bands) == 0:
            raise RuntimeError("no valid bands were selected!")

        keep_rows = np.ones(len(self.fibermap), bool)
        for fm_select,fm_var in zip([nights, exposures, targets, fibers],
                                    ['NIGHT', 'EXPID', 'TARGETID', 'FIBER']):
            if fm_select is not None:
                keep_selection = np.isin(self.fibermap[fm_var], fm_select)
                if sum(keep_selection) == 0:
                    raise RuntimeError("no valid "+fm_var+" were selected!")
                keep_rows = keep_rows & keep_selection

        if invert:
            keep_rows = np.invert(keep_rows)

        keep, = np.where(keep_rows)
        if len(keep) == 0:
            raise RuntimeError("selection has no spectra")

        sp = self._get_slice(keep_rows, bands=keep_bands)

        if return_index:
            return (sp, keep)

        return sp
    
    
    def coadd_fibermap_no_expid(fibermap) :
        #log = get_logger()
        #log.debug("'coadding' fibermap")

        targets = np.unique(fibermap["TARGETID"])
        ntarget = targets.size

        jj=np.zeros(ntarget,dtype=int)
        for i,tid in enumerate(targets) :
            jj[i]=np.where(fibermap["TARGETID"]==tid)[0][0]
        tfmap=fibermap[jj]

        #- initialize NUMEXP=-1 to check that they all got filled later
        tfmap['COADD_NUMEXP'] = np.zeros(len(tfmap), dtype=np.int16) - 1
        tfmap['COADD_EXPTIME'] = np.zeros(len(tfmap), dtype=np.float32) - 1

        # smarter values for some columns
        mean_cols = [
            'DELTA_X', 'DELTA_Y',
            'FIBER_X', 'FIBER_Y',
            'FIBER_RA', 'FIBER_DEC',
            'FIBERASSIGN_X', 'FIBERASSIGN_Y'
            ]
        rms_cols = ['DELTA_X', 'DELTA_Y']  #- rms_cols must also be in mean_cols
        for k in mean_cols:
            if k in fibermap.colnames :
                if k.endswith('_RA') or k.endswith('_DEC'):
                    dtype = np.float64
                else:
                    dtype = np.float32
                if k in mean_cols:
                    xx = Column(np.zeros(ntarget, dtype=dtype))
                    tfmap.add_column(xx,name='MEAN_'+k)
                if k in rms_cols:
                    xx = Column(np.zeros(ntarget, dtype=dtype))
                    tfmap.add_column(xx,name='RMS_'+k)

                tfmap.remove_column(k)

        first_last_cols = ['NIGHT','EXPID','TILEID','SPECTROID','FIBER','MJD']
        for k in first_last_cols:
            if k in fibermap.colnames :
                if k in ['MJD']:
                    dtype = np.float32
                else:
                    dtype = np.int32
                if not 'FIRST_'+k in tfmap.dtype.names :
                    xx = Column(np.arange(ntarget, dtype=dtype))
                    tfmap.add_column(xx,name='FIRST_'+k)
                if not 'LAST_'+k in tfmap.dtype.names :
                    xx = Column(np.arange(ntarget, dtype=dtype))
                    tfmap.add_column(xx,name='LAST_'+k)
                if not 'NUM_'+k in tfmap.dtype.names :
                    xx = Column(np.arange(ntarget, dtype=np.int16))
                    tfmap.add_column(xx,name='NUM_'+k)

        for i,tid in enumerate(targets) :
            jj = fibermap["TARGETID"]==tid

            #- coadded FIBERSTATUS = bitwise AND of input FIBERSTATUS
            tfmap['FIBERSTATUS'][i] = np.bitwise_and.reduce(fibermap['FIBERSTATUS'][jj])

            #- Only FIBERSTATUS=0 were included in the coadd
            fiberstatus_nonamp_bits = get_all_nonamp_fiberbitmask_val()
            fiberstatus_amp_bits = get_justamps_fiberbitmask()
            targ_fibstatuses = fibermap['FIBERSTATUS'][jj]
            nonamp_fiberstatus_flagged = ( (targ_fibstatuses & fiberstatus_nonamp_bits) > 0 )
            allamps_flagged = ( (targ_fibstatuses & fiberstatus_amp_bits) == fiberstatus_amp_bits )
            good_coadds = np.bitwise_not( nonamp_fiberstatus_flagged | allamps_flagged )
            tfmap['COADD_NUMEXP'][i] = np.count_nonzero(good_coadds)
            if 'EXPTIME' in fibermap.colnames :
                tfmap['COADD_EXPTIME'][i] = np.sum(fibermap['EXPTIME'][jj][good_coadds])
            for k in mean_cols:
                if k in fibermap.colnames :
                    vals=fibermap[k][jj]
                    tfmap['MEAN_'+k][i] = np.mean(vals)

            for k in rms_cols:
                if k in fibermap.colnames :
                    vals=fibermap[k][jj]
                    # RMS includes mean offset, not same as std
                    tfmap['RMS_'+k][i] = np.sqrt(np.mean(vals**2))

            for k in first_last_cols:
                if k in fibermap.colnames :
                    vals=fibermap[k][jj]
                    tfmap['FIRST_'+k][i] = np.min(vals)
                    tfmap['LAST_'+k][i] = np.max(vals)
                    tfmap['NUM_'+k][i] = np.unique(vals).size

            for k in ['FIBER_RA_IVAR', 'FIBER_DEC_IVAR','DELTA_X_IVAR', 'DELTA_Y_IVAR'] :
                if k in fibermap.colnames :
                    tfmap[k][i]=np.sum(fibermap[k][jj])

        #- Remove some columns that apply to individual exp but not coadds
        #- (even coadds of the same tile)
        #for k in ['NIGHT', 'EXPID', 'MJD', 'EXPTIME', 'NUM_ITER']:
        #    if k in tfmap.colnames:
        #        tfmap.remove_column(k)

        return tfmap
    
    def read_spectra_no_expid(infile, single=False):
        """
        Read Spectra object from FITS file.

        This reads data written by the write_spectra function.  A new Spectra
        object is instantiated and returned.

        Args:
            infile (str): path to read
            single (bool): if True, keep spectra as single precision in memory.

        Returns (Spectra):
            The object containing the data read from disk.

        """
        log = get_logger()
        ftype = np.float64
        if single:
            ftype = np.float32

        infile = os.path.abspath(infile)
        if not os.path.isfile(infile):
            raise IOError("{} is not a file".format(infile))

        t0 = time.time()
        hdus = fits.open(infile, mode="readonly")
        nhdu = len(hdus)

        # load the metadata.

        meta = dict(hdus[0].header)

        # initialize data objects

        bands = []
        fmap = None
        expfmap = None
        wave = None
        flux = None
        ivar = None
        mask = None
        res = None
        extra = None
        extra_catalog = None
        scores = None

        # For efficiency, go through the HDUs in disk-order.  Use the
        # extension name to determine where to put the data.  We don't
        # explicitly copy the data, since that will be done when constructing
        # the Spectra object.

        for h in range(1, nhdu):
            name = hdus[h].header["EXTNAME"]
            if name == "FIBERMAP":
                fmap = encode_table(Table(hdus[h].data, copy=True).as_array())
            elif name == "EXP_FIBERMAP":
                expfmap = encode_table(Table(hdus[h].data, copy=True).as_array())
            elif name == "SCORES":
                scores = encode_table(Table(hdus[h].data, copy=True).as_array())
            elif name == 'EXTRA_CATALOG':
                extra_catalog = encode_table(Table(hdus[h].data, copy=True).as_array())
            else:
                # Find the band based on the name
                mat = re.match(r"(.*)_(.*)", name)
                if mat is None:
                    raise RuntimeError("FITS extension name {} does not contain the band".format(name))
                band = mat.group(1).lower()
                type = mat.group(2)
                if band not in bands:
                    bands.append(band)
                if type == "WAVELENGTH":
                    if wave is None:
                        wave = {}
                    wave[band] = native_endian(hdus[h].data.astype(ftype))
                elif type == "FLUX":
                    if flux is None:
                        flux = {}
                    flux[band] = native_endian(hdus[h].data.astype(ftype))
                elif type == "IVAR":
                    if ivar is None:
                        ivar = {}
                    ivar[band] = native_endian(hdus[h].data.astype(ftype))
                elif type == "MASK":
                    if mask is None:
                        mask = {}
                    mask[band] = native_endian(hdus[h].data.astype(np.uint32))
                elif type == "RESOLUTION":
                    if res is None:
                        res = {}
                    res[band] = native_endian(hdus[h].data.astype(ftype))
                else:
                    # this must be an "extra" HDU
                    if extra is None:
                        extra = {}
                    if band not in extra:
                        extra[band] = {}
                    extra[band][type] = native_endian(hdus[h].data.astype(ftype))

        hdus.close()
        duration = time.time() - t0
        log.info(iotime.format('read', infile, duration))

        # Construct the Spectra object from the data.  If there are any
        # inconsistencies in the sizes of the arrays read from the file,
        # they will be caught by the constructor.

        spec = spectra_no_expid(bands, wave, flux, ivar, mask=mask, resolution_data=res,
            fibermap=fmap, exp_fibermap=expfmap,
            meta=meta, extra=extra, extra_catalog=extra_catalog,
            single=single, scores=scores)

        return spec
    
    def update_no_expid(self, other): 
        """
        Overwrite or append new data.

        Given another Spectra object, compare the fibermap information with
        the existing one.  For spectra that already exist, overwrite existing 
        data with the new values.  For spectra that do not exist, append that 
        data to the end of the spectral data.

        Args:
            other (Spectra): the new data to add.

        Returns:
            nothing (object updated in place).

        Note: if fibermap, scores and extra_catalog exist in the new data, they
        are appended to the existing tables. If those new tables have different columns,
        only columns with identical names will be appended. Spectra.meta is unchanged.
        """

        # Does the other Spectra object have any data?

        if other.num_spectra() == 0:
            return

        # Do we have new bands to add?

        newbands = []
        for b in other.bands:
            if b not in self.bands:
                newbands.append(b)
            else:
                if not np.allclose(self.wave[b], other.wave[b]):
                    raise RuntimeError("band {} has an incompatible wavelength grid".format(b))

        bands = list(self.bands)
        bands.extend(newbands)

        # Are we adding mask data in this update?

        add_mask = False
        if other.mask is None:
            if self.mask is not None:
                raise RuntimeError("existing spectra has a mask, cannot "
                    "update it to a spectra with no mask")
        else:
            if self.mask is None:
                add_mask = True

        # Are we adding resolution data in this update?

        ndiag = {}

        add_res = False
        if other.resolution_data is None:
            if self.resolution_data is not None:
                raise RuntimeError("existing spectra has resolution data, cannot "
                    "update it to a spectra with none")
        else:
            if self.resolution_data is not None:
                for b in self.bands:
                    ndiag[b] = self.resolution_data[b].shape[1]
                for b in other.bands:
                    odiag = other.resolution_data[b].shape[1]
                    if b not in self.bands:
                        ndiag[b] = odiag
                    else:
                        if odiag != ndiag[b]:
                            raise RuntimeError("Resolution matrices for a"
                                " given band must have the same dimensoins")
            else:
                add_res = True
                for b in other.bands:
                    ndiag[b] = other.resolution_data[b].shape[1]

        # Are we adding extra data in this update?

        add_extra = False
        if other.extra is None:
            if self.extra is not None:
                raise RuntimeError("existing spectra has extra data, cannot "
                    "update it to a spectra with none")
        else:
            if self.extra is None:
                add_extra = True

        # Compute which targets / exposures are new

        nother = len(other.fibermap)
        exists = np.zeros(nother, dtype=np.int)

        indx_original = []

        if ( (self.fibermap is not None) and
            all([x in fm.keys() for x in ['EXPID', 'FIBER']
                                for fm in [self.fibermap, other.fibermap]]) ):
            for r in range(nother):
                expid = other.fibermap[r]["EXPID"]
                fiber = other.fibermap[r]["FIBER"]
                for i, row in enumerate(self.fibermap):
                    if (expid == row["EXPID"]) and (fiber == row["FIBER"]):
                        indx_original.append(i)
                        exists[r] += 1

        if len(np.where(exists > 1)[0]) > 0:
            #raise RuntimeError("found duplicate spectra (same FIBER) in the fibermap")
            raise RuntimeError("found duplicate spectra (same EXPID and FIBER) in the fibermap")

        indx_exists = np.where(exists == 1)[0]
        indx_new = np.where(exists == 0)[0]

        # Make new data arrays of the correct size to hold both the old and 
        # new data

        nupdate = len(indx_exists)
        nnew = len(indx_new)

        if self.fibermap is None:
            nold = 0
            newfmap = other.fibermap.copy()
        else:
            nold = len(self.fibermap)
            newfmap = encode_table(np.zeros( (nold + nnew, ),
                                   dtype=self.fibermap.dtype))
        
        newscores = None
        if self.scores is not None:
            newscores = encode_table(np.zeros( (nold + nnew, ),
                                   dtype=self.scores.dtype))

        newextra_catalog = None
        if self.extra_catalog is not None:
            newextra_catalog = encode_table(np.zeros( (nold + nnew, ),
                                   dtype=self.extra_catalog.dtype))

        newwave = {}
        newflux = {}
        newivar = {}
        
        newmask = None
        if add_mask or self.mask is not None:
            newmask = {}
        
        newres = None
        newR = None
        if add_res or self.resolution_data is not None:
            newres = {}
            newR = {}

        newextra = None
        if add_extra or self.extra is not None:
            newextra = {}

        for b in bands:
            nwave = None
            if b in self.bands:
                nwave = self.wave[b].shape[0]
                newwave[b] = self.wave[b]
            else:
                nwave = other.wave[b].shape[0]
                newwave[b] = other.wave[b].astype(self._ftype)
            newflux[b] = np.zeros( (nold + nnew, nwave), dtype=self._ftype)
            newivar[b] = np.zeros( (nold + nnew, nwave), dtype=self._ftype)
            if newmask is not None:
                newmask[b] = np.zeros( (nold + nnew, nwave), dtype=np.uint32)
                #newmask[b][:,:] = specmask["NODATA"]
            if newres is not None:
                newres[b] = np.zeros( (nold + nnew, ndiag[b], nwave), dtype=self._ftype)
            if newextra is not None:
                newextra[b] = {}

        # Copy the old data

        if nold > 0:
            # We have some data (i.e. we are not starting with an empty Spectra)
            for newtable, original_table in zip([newfmap, newscores, newextra_catalog],
                                           [self.fibermap, self.scores, self.extra_catalog]):
                if original_table is not None:
                    newtable[:nold] = original_table

            for b in self.bands:
                newflux[b][:nold,:] = self.flux[b]
                newivar[b][:nold,:] = self.ivar[b]
                if self.mask is not None:
                    newmask[b][:nold,:] = self.mask[b]
                elif add_mask:
                    newmask[b][:nold,:] = 0
                if self.resolution_data is not None:
                    newres[b][:nold,:,:] = self.resolution_data[b]
                if self.extra is not None:
                    for ex in self.extra[b].items():
                        newextra[b][ex[0]] = np.zeros( newflux[b].shape,
                            dtype=self._ftype)
                        newextra[b][ex[0]][:nold,:] = ex[1]

        # Update existing spectra

        for i, s in enumerate(indx_exists):
            row = indx_original[i]
            for b in other.bands:
                newflux[b][row,:] = other.flux[b][s,:].astype(self._ftype)
                newivar[b][row,:] = other.ivar[b][s,:].astype(self._ftype)
                if other.mask is not None:
                    newmask[b][row,:] = other.mask[b][s,:]
                else:
                    newmask[b][row,:] = 0
                if other.resolution_data is not None:
                    newres[b][row,:,:] = other.resolution_data[b][s,:,:].astype(self._ftype)
                if other.extra is not None:
                    for ex in other.extra[b].items():
                        if ex[0] not in newextra[b]:
                            newextra[b][ex[0]] = np.zeros(newflux[b].shape,
                                dtype=self._ftype)
                        newextra[b][ex[0]][row,:] = ex[1][s,:].astype(self._ftype)

        # Append new spectra

        if nnew > 0:
            for newtable, othertable in zip([newfmap, newscores, newextra_catalog],
                                           [other.fibermap, other.scores, other.extra_catalog]):
                if othertable is not None:
                    if newtable.dtype == othertable.dtype:
                        newtable[nold:] = othertable[indx_new]
                    else:
                    #- if table contents do not match, still merge what we can, based on key names
                    # (possibly with numpy automatic casting)
                        for k in set(newtable.keys()).intersection(set(othertable.keys())):
                            newtable[k][nold:] = othertable[k][indx_new]

            for b in other.bands:
                newflux[b][nold:,:] = other.flux[b][indx_new].astype(self._ftype)
                newivar[b][nold:,:] = other.ivar[b][indx_new].astype(self._ftype)
                if other.mask is not None:
                    newmask[b][nold:,:] = other.mask[b][indx_new]
                else:
                    newmask[b][nold:,:] = 0
                if other.resolution_data is not None:
                    newres[b][nold:,:,:] = other.resolution_data[b][indx_new].astype(self._ftype)
                if other.extra is not None:
                    for ex in other.extra[b].items():
                        if ex[0] not in newextra[b]:
                            newextra[b][ex[0]] = np.zeros(newflux[b].shape,
                                dtype=self._ftype)
                        newextra[b][ex[0]][nold:,:] = ex[1][indx_new].astype(self._ftype)

        # Update all sparse resolution matrices

        for b in bands:
            if newres is not None:
                newR[b] = np.array( [ Resolution(r) for r in newres[b] ] )

        # Swap data into place

        self._bands = bands
        self.wave = newwave
        self.fibermap = newfmap
        self.flux = newflux
        self.ivar = newivar
        self.mask = newmask
        self.resolution_data = newres
        self.R = newR
        self.extra = newextra
        self.scores = newscores
        self.extra_catalog = newextra_catalog

        return