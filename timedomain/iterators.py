from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

import numpy as np
import copy

from . import fs_utils

"""
Class that returns Tiles that were observed on a date
"""
class Date_Tile_Iterator:
    
    def __init__(self, date, subdir='andes'):
        self.date = date
        self.tiles = fs_utils.dateToTiles(date, subdir=subdir)
        self.index = 0
        
    def __iter__(self):
        
        if (len(self.tiles) == 0):
            self.index=1
        else:
            self.index = 0
        return self
    
    def __next__(self):
        
        if self.index < len(self.tiles):
            ans = self.tiles[self.index]
            print("Tile {}".format(ans))
            self.index = self.index+1
            return ans
        else:
            raise StopIteration
            
    def hasNext(self):
        if self.index < len(self.tiles):
            return True
        else:
            return False
            
"""
Class that returns all spectra frome a date
"""
class Date_Spectra_Iterator:
       
    def __init__(self, date, subdir='andes', trunk='coadd'):
        self.date = date
        self.subdir = subdir
        self.trunk = trunk
        self.date_tile_iterator = Date_Tile_Iterator(self.date, subdir=subdir)
        self.panel_iterator = np.nditer(fs_utils.panels) # Panel_Iterator()
        
        self.tile = None
        self.panel = None

    def __iter__(self):
        
        self.date_tile_iterator = Date_Tile_Iterator(self.date, subdir=self.subdir)
        self.panel_iterator = np.nditer(fs_utils.panels) # Panel_Iterator()
        self.tile = None
        self.panel = None
        return self
    
    def __next__(self):
        
        while True:
            if self.tile is None:
                self.tile = self.date_tile_iterator.__next__()
                
            try:
                self.panel = self.panel_iterator.__next__()
                fname = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                if fname is not None:
                    break
                    
            except StopIteration:
                if self.date_tile_iterator.hasNext():
                    self.tile = self.date_tile_iterator.__next__()
                    self.panel_iterator = np.nditer(fs_utils.panels) #Panel_Iterator()
                    self.panel = self.panel_iterator.__next__()
                    fname = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                    if fname is not None:
                        break
                else:
                    raise StopIteration

        print("Filename {}".format(fname))
        return read_spectra(fname), fname
    
    @staticmethod
    def test():
        date = "20201223"
        
        it = Date_Spectra_Iterator(date,subdir='daily', trunk='coadd')
        for a in it:
            print(a)

"""
Class that returns subspectra by target id from spectra
"""
class Spectra_Subspectra_Iterator:
    def __init__(self, spectrum):
        self.spectrum = spectrum
        self.titerator = np.nditer(np.array(self.spectrum.target_ids()))

    def __iter__(self):
        self.titerator = np.nditer(np.array(self.spectrum.target_ids()))
        return self
    
    def __next__(self):
        
        try:
            ans = self.spectrum.select(targets=np.array(self.titerator.__next__()))
            return ans
        except StopIteration:
            raise StopIteration
    

## Difference Iterators

"""
Class that returns all Spectra pairs from Spectra
"""
class Spectra_Pairs_Iterator:
    
    def __init__(self, spectra):
        self.spectra = spectra
        self.i=0
        self.j=1
        self.spectrai=None


    def __iter__(self):
        self.i=0
        self.j=1
        self.spectrai=None
        return self
    
    def subspectra(self, index):
        ans  = copy.deepcopy(self.spectra)
        for band in ans.bands:
            ans.flux[band] = self.spectra.flux[band][index,:]
            ans.ivar[band] = self.spectra.ivar[band][index,:]
            ans.mask[band] =   self.spectra.mask[band][index,:]
        ans.fibermap = ans.fibermap[[index]]

        return ans
            
    def __next__(self):

        while self.i< self.spectra.num_spectra():
        
            if self.spectrai is None:
                self.spectrai = self.subspectra(self.i)    
                
            if self.j < self.spectra.num_spectra():
                ans = (self.spectrai,self.subspectra(self.j))
                print("Targets {} {} of {}".format(self.i, self.j,self.spectra.num_spectra()))
                self.j = self.j+1
                break
            else:
                self.i = self.i+1
                self.j=self.i+1
                self.spectrai = None
        else:
            raise StopIteration
            
        return ans

class PairCoadds:

    def __init__(self, tile, subdir='andes',dates=None):
        self.tile = tile
        self.subdir=subdir
        if dates is None:
            self.dates = fs_utils.tileToDates(tile,subdir=subdir)
        else:
            self.dates = dates
        if (len(self.dates) <=1):
            self.a=1
        else:
            self.a = 0
        self.b = 1
        self.panel = 0
        
    def __iter__(self):
        if (len(self.dates) <=1):
            self.a=1
        else:
            self.a = 0
        self.b = 1
        self.panel = 0
        return self
    
    def __next_index__(self):
        if self.panel+1 < len(fs_utils.panels):
            self.panel +=1
        else:
            self.panel=0
            if self.b +1 < len(self.dates):
                self.b+=1
            else:
                self.a+=1
                self.b = self.a+1

    def __next__(self):
           
        while self.a+1 < len(self.dates):
            df0 = fs_utils.fitsfile(self.tile,self.dates[self.a],fs_utils.panels[self.panel])
            df1 = fs_utils.fitsfile(self.tile,self.dates[self.b],fs_utils.panels[self.panel])

            if (df0 is not None and df1 is not None):
                ans =  (read_spectra(df0),read_spectra(df1))
                self.__next_index__()
                break
        
            self.__next_index__()

        else:
            raise StopIteration
        
        return ans
    
    #staticmethod
    def test():
        tile = "70006"
        for c in PairCoadds(tile):
            print (c)
            
class DailyPairs:

    def __init__(self, date, subdir='daily'):
        self.date = date
        self.tiles = fs_utils.dateToTiles(date, subdir='daily')
        self.subdir=subdir
        
        self.spec = None
        self.fibermap = None
        self.tile_index = 0
        self.panel = 0
        
    def __iter__(self):
        self.panel = 0
        self.tile_index = 0
        self.df = None
        return self

    def __next__(self):
        
        # if need a new file
        if self.panel == 0:
            df = fs_utils.fitsfile(self.tiles[self.tile_index],self.date,fs_utils.panels[self.panel],trunk='spectra')
            print(df)
            self.spec = read_spectra(df)
            self.fibermap = self.spec.fibermap
            
        print(self.fibermap)
        wef
            
        while self.a+1 < len(self.dates):
            df0 = fs_utils.fitsfile(tile,self.dates[self.a],fs_utils.panels[self.panel])
            df1 = fs_utils.fitsfile(tile,self.dates[self.b],fs_utils.panels[self.panel])

            if (df0 is not None and df1 is not None):
                ans = (read_spectra(df0), read_spectra(df1))
#                 ans =  (df0,df1)
                self.__next_index__()
                break
        
            self.__next_index__()

        else:
            raise StopIteration
            
#         while self.a+1 < len(self.dates):
#             df0 = fs_utils.fitsfile(tile,self.dates[self.a],fs_utils.panels[self.panel])
#             df1 = fs_utils.fitsfile(tile,self.dates[self.b],fs_utils.panels[self.panel])

#             if (df0 is not None and df1 is not None):
#                 ans = (read_spectra(df0), read_spectra(df1))
# #                 ans =  (df0,df1)
#                 self.__next_index__()
#                 break
        
#             self.__next_index__()

#         else:
#             raise StopIteration
        
        print(self.tile_index,self.panel)
        return ans

    #staticmethod
    def test():
        for c in PairSingle("70006","20200305"):
            print (c)

# PairSingle.test()

class TileDate:
    
    def __init__(self, tile, date, subdir='andes'):
        dirname = os.path.join(redux,subdir,'tiles',tile,date)
        if not os.path.isdir(dirname):
            print('{} does not exist.'.format(dirname))    
        self.tile = tile
        self.date = date
        self.subdir = subdir
        cafiles =  glob(os.path.join(dirname,'cframe-??-*.fits'))
        exposures=[]

        for cafile in cafiles:
            index0 = cafile.find('cframe-')
            index1 = cafile.find('.fits',index0+10)
            exposures.append(cafile[index0+10:index1])

        exposures = np.unique(exposures)
        print(exposures)
        
    def __iter__(self):
        self.a = 1
        return self

    def __next__(self):
        x = self.a
        self.a += 1
        return x
    
# tile = "70006"
# date = "20200305"

