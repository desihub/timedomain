import os
from glob import glob
import numpy as np
from . import config

redux='/global/project/projectdirs/desi/spectro/redux/'
filts = ['b','r','z']
panels = np.arange(10).astype('str')

def recentRelease(redux,tile, date, panel, trunk):
    ans = None
    for testdir in config.releases:
        name  = os.path.join(redux,testdir,'tiles',tile, date, '{}-{}-{}-{}.fits'.format(trunk,panel,tile,date))
        if os.path.exists(name):
            ans = testdir
    if ans is None:
#         name  = os.path.join(redux,'daily','tiles',tile, date, '{}-{}-{}-{}.fits'.format(trunk,panel,tile,date))
#         if os.path.exists(name):
        ans = 'daily'
    return ans
        

def fitsfile(tile, date, panel, subdir='recent',trunk='coadd'):
    
    if subdir == 'recent':
        subdir = recentRelease(redux,tile,date, panel, trunk)
    name  = os.path.join(redux,subdir,'tiles',tile, date, '{}-{}-{}-{}.fits'.format(trunk,panel,tile,date))
    exists = os.path.exists(name)
    if exists:
        return os.path.join(redux,subdir,'tiles',tile, date, '{}-{}-{}-{}.fits'.format(trunk,panel,tile,date))
    else:
        return None

def tiledateToExposures(tile, date, subdir='andes'):

    dirname = os.path.join(redux,subdir,'tiles',tile,date)
    if not os.path.isdir(dirname):
        print('{} does not exist.'.format(dirname))    
    cafiles =  glob(os.path.join(dirname,'cframe-??-*.fits'))
    exposures=[]

    for cafile in cafiles:
        index0 = cafile.find('cframe-')
        index1 = cafile.find('.fits',index0+10)
        exposures.append(cafile[index0+10:index1])

    exposures = np.unique(exposures)
    return exposures

# tile = "70006"
# date = "20200305"    
# tn=tiledateToExposures(tile,date)

def dateToTiles(date, subdir='andes'):

    dirname = os.path.join(redux,subdir,'tiles')
    if not os.path.isdir(dirname):
        print('{} does not exist.'.format(dirname))
    cafiles =  glob(dirname+'/*/'+date)
    tiles=[]

    for cafile in cafiles:
        index0 = cafile.find('/'+date)
        index1 = cafile.rfind('/',0,index0-1)
        tiles.append(cafile[index1+1:index0])

    return tiles

# date = "20200305"
# print(dateToTiles(date))

def tileToDates(tile, subdir='andes'):

    dirname = os.path.join(redux,subdir,'tiles',tile)
    if not os.path.isdir(dirname):
        print('{} does not exist.'.format(dirname))

    dates = [f.path for f in os.scandir(dirname) if f.is_dir()]

    for i in range(len(dates)):
        dates[i] = dates[i][dates[i].rfind('/')+1:]

    return dates

# tile = "70006"
# print(tileToDates(tile))
