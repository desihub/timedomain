#!/usr/bin/env python
"""Apply the DESITrIP CNN classifier to observed spectra,
chosen by tile ID and date.
"""

from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra
from desitarget.cmx.cmx_targetmask import cmx_mask

from desitrip.preproc import rebin_flux, rescale_flux

from astropy.io import fits
from astropy.table import Table, vstack, hstack

from glob import glob
from datetime import date
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from tensorflow import keras

p = ArgumentParser(description='DESITrIP data processing',
                   formatter_class=ArgumentDefaultsHelpFormatter)
p.add_argument('--tile', type=int, default=0,
               help='Tile ID for processing.')
p.add_argument('--date', default=date.today().strftime('%Y%m%d'),
               help='Date of observation [YYYYMMDD]')
p.add_argument('--tfmodel', default=None,
               help='TensorFlow model HDF5 definition')
args = p.parse_args()

# Access redux folder.
redux='/global/project/projectdirs/desi/spectro/redux/daily/tiles'
prefix_in='/'.join([redux, '{:05d}'.format(args.tile), args.date]) 
if not os.path.isdir(prefix_in):
    raise SystemExit('{} does not exist.'.format(prefix_in))

# Set up BGS target bit selection.
cmx_bgs_bits = '|'.join([_ for _ in cmx_mask.names() if 'BGS' in _])

# List zbest and coadd files.
zbfiles = sorted(glob('{}/zbest*.fits'.format(prefix_in)))
cafiles = sorted(glob('{}/coadd*.fits'.format(prefix_in)))

if args.tfmodel is not None:
    classifier = keras.models.load_model(args.tfmodel)

# Loop through zbest and coadd files for each petal.
# Extract the fibermaps, ZBEST tables, and spectra.
# Keep only BGS targets passing basic event selection.
allzbest = None
allfmap = None
allwave = None
allflux = None
allivar = None
allmask = None
allres  = None

for cafile, zbfile in zip(cafiles, zbfiles):
    # Access data per petal.
    zbest = Table.read(zbfile, 'ZBEST')
    fibermap = Table.read(zbfile, 'FIBERMAP')
    pspectra = read_spectra(cafile)

    # Apply standard event selection.
    isTGT = fibermap['OBJTYPE'] == 'TGT'
    isGAL = zbest['SPECTYPE'] == 'GALAXY'
    isBGS = fibermap['CMX_TARGET'] & cmx_mask.mask(cmx_bgs_bits) != 0
    select = isTGT & isGAL & isBGS

    # Accumulate spectrum data.
    if allzbest is None:
        allzbest = zbest[select]
        allfmap = fibermap[select]
        allwave = pspectra.wave['brz']
        allflux = pspectra.flux['brz'][select]
        allivar = pspectra.ivar['brz'][select]
        allmask = pspectra.mask['brz'][select]
        allres  = pspectra.resolution_data['brz'][select]
    else:
        allzbest = vstack([allzbest, zbest[select]])
        allfmap = vstack([allfmap, fibermap[select]])
        allflux = np.vstack([allflux, pspectra.flux['brz'][select]])
        allivar = np.vstack([allivar, pspectra.ivar['brz'][select]])
        allmask = np.vstack([allmask, pspectra.mask['brz'][select]])
        allres  = np.vstack([allres, pspectra.resolution_data['brz'][select]])

# Apply the DESITrIP preprocessing to selected spectra.
rewave, reflux, reivar = rebin_flux(allwave, allflux, allivar, allzbest['Z'],
                                    minwave=2500., maxwave=9500., nbins=150,
                                    log=True, clip=True)
rsflux = rescale_flux(reflux)

# Run the classification.
if args.tfmodel is not None:
    pred = classifier.predict(rsflux)

# Create output: selected target spectra.
selected_spectra = Spectra(bands=['brz'],
                           wave={'brz' : allwave},
                           flux={'brz' : allflux},
                           ivar={'brz' : allivar},
                           mask={'brz' : allmask},
                           resolution_data={'brz' : allres},
                           fibermap=allfmap)

write_spectra('selected-{}-{}.fits'.format(args.tile, args.date), selected_spectra)

# Append preprocess spectra to output.
hx = fits.HDUList()

hdu_rewave = fits.PrimaryHDU(rewave)
hdu_rewave.header['EXTNAME'] = 'REWAVE'
hdu_rewave.header['BUNIT'] = 'Angstrom'
hdu_rewave.header['AIRORVAC'] = ('vac', 'Vacuum wavelengths')
hx.append(hdu_rewave)

hdu_reflux = fits.ImageHDU(reflux)
hdu_reflux.header['EXTNAME'] = 'REFLUX'
hx.append(hdu_reflux)

hdu_rsflux = fits.ImageHDU(rsflux)
hdu_rsflux.header['EXTNAME'] = 'RSFLUX'
hx.append(hdu_rsflux)

hdu_classify = fits.ImageHDU(pred)
hdu_classify.header['EXTNAME'] = 'OBJCLASS'
hx.append(hdu_classify)

hx.append(fits.BinTableHDU(allzbest))
hx.writeto('reduced-{}-{}.fits'.format(args.tile, args.date), overwrite=True)
