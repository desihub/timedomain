import os
from glob import glob
import numpy as np
from . import config

redux='/global/project/projectdirs/desi/spectro/redux/'
filts = ['b','r','z']
panels = np.arange(10).astype('str')

# this set of images officially dead to DESI
bad = [["20210228","80726"],["20210228","80726"],["20210228","80740"],["20210228","80741"]]

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

def useSubdir(subdir,date="99999999"):
    if subdir == 'recent':
        if date <= config.last_release_date:
            return  config.last_release
        else:
            return "daily"   
    return subdir
        

def fitsfile(tile, date, panel, subdir='recent',trunk='coadd'):
    
    if [date,tile] in bad:
        return None
    
    usesubdir=useSubdir(subdir,date)

    name  = os.path.join(redux,usesubdir,'tiles',tile, date, '{}-{}-{}-{}.fits'.format(trunk,panel,tile,date))
    exists = os.path.exists(name)
    if exists:
        return os.path.join(redux,usesubdir,'tiles',tile, date, '{}-{}-{}-{}.fits'.format(trunk,panel,tile,date))
    else:
        return None

def tiledateToExposures(tile, date, subdir='andes'):
    
    usesubdir=useSubdir(subdir,date)


    dirname = os.path.join(redux,usesubdir,'tiles',tile,date)
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
    usesubdir=useSubdir(subdir,date)

    dirname = os.path.join(redux,usesubdir,'tiles')
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

    usesubdir=useSubdir(subdir)

    dirname = os.path.join(redux,usesubdir,'tiles',tile)
    if not os.path.isdir(dirname):
        print('{} does not exist.'.format(dirname))

    dates = [f.path for f in os.scandir(dirname) if f.is_dir()]

    for i in range(len(dates)):
        dates[i] = dates[i][dates[i].rfind('/')+1:]

    return dates

# tile = "70006"
# print(tileToDates(tile))
