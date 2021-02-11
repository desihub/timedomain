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
Class that returns Tiles that were observed on a date
"""
class TileDate_PreDate_Iterator:
    
    def __init__(self, tile, date, subdir='andes'):
        self.tile = tile
        self.date = date
        self.subdir= subdir
        
        self.it0 = None
        
    def __iter__(self):
        
        self.it0 = None
        return self
    
    def __next__(self):
        
        if self.it0 is None:
            predates = fs_utils.tileToDates(self.tile, subdir=self.subdir)
            predates = np.array(predates,dtype='str')
            w = predates < self.date
            predates = predates[w]
            if len(predates)==0:
                raise StopIteration
            self.it0 = predates.flat
            return self.it0.__next__()
        
        try:
            return self.it0.__next__()
        except StopIteration:
            raise StopIteration

    @staticmethod
    def test():
        tile="80623"
        date = "20201223"
        it = TileDate_PreDate_Iterator(tile,date,subdir='daily')
        for a in it:
            print(type(a),a.shape,a.astype('str'))
            
"""
Class that returns all spectra frome a date
"""
class Date_Spectra_Iterator:
       
    def __init__(self, date, subdir='andes', trunk='coadd', verbose=False):
        self.date = date
        self.subdir = subdir
        self.trunk = trunk
        self.verbose = verbose
        
        self.date_tile_iterator = Date_Tile_Iterator(self.date, subdir=subdir)
        self.panel_iterator = fs_utils.panels.flat #panel_Iterator()
        
        self.tile = None
        self.panel = None

    def __iter__(self):
        self.date_tile_iterator = Date_Tile_Iterator(self.date, subdir=self.subdir)
        self.panel_iterator = fs_utils.panels.flat # Panel_Iterator()
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
                    self.panel_iterator = fs_utils.panels.flat #Panel_Iterator()
                    self.panel = self.panel_iterator.__next__()
                    fname = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                    if fname is not None:
                        break
                else:
                    raise StopIteration

        if self.verbose:
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
    def __init__(self, spectrum, verbose=False):
        self.spectrum = spectrum
        self.titerator = np.array(self.spectrum.target_ids()).flat
        self.verbose = verbose
        
    def __iter__(self):
        self.titerator = np.array(self.spectrum.target_ids()).flat
        return self
    
    def __next__(self):
        
        try:
            ans = self.spectrum.select(targets=np.array(self.titerator.__next__()))
            if self.verbose:
                print("TARGETID {}".format(ans.fibermap["TARGETID"][0]))
            return ans
        except StopIteration:
            raise StopIteration
    

## Difference Iterators
"""

Given a Spectra return all pairs

"""
class Spectra_Pairs_Iterator:
    
    def __init__(self, spectra, verbose=False):
        self.spectra = spectra
        self.i=0
        self.j=1
        self.spectrai=None
        self.verbose = verbose


    def __iter__(self):
        self.i=0
        self.j=1
        self.spectrai=None
        return self
    
    def subspectra(self, index):
        ans  = copy.deepcopy(self.spectra)
        for band in ans.bands:
            ans.flux[band] = self.spectra.flux[band][index,:][np.newaxis,:]
#             print(self.spectra.flux[band].shape,self.spectra.flux[band][index,:].shape,ans.flux[band].shape)
            ans.ivar[band] = self.spectra.ivar[band][index,:][np.newaxis,:]
            ans.mask[band] =   self.spectra.mask[band][index,:][np.newaxis,:]
        ans.fibermap = ans.fibermap[[index]]

        return ans
            
    def __next__(self):

        while self.i< self.spectra.num_spectra():
        
            if self.spectrai is None:
                self.spectrai = self.subspectra(self.i)    
                
            if self.j < self.spectra.num_spectra():
                ans = (self.spectrai,self.subspectra(self.j))
#                 print("Targets {} {} of {}".format(self.i, self.j,self.spectra.num_spectra()))
                self.j = self.j+1
                break
            else:
                self.i = self.i+1
                self.j=self.i+1
                self.spectrai = None
        else:
            raise StopIteration
        
        if self.verbose:
            print("EXPID 1: {}  EXPID 2: {}".format(ans[0].fibermap["EXPID"].data, ans[1].fibermap["EXPID"].data))
        
        return ans

    
"""
All pairs of the same target taken on a night
"""
class Date_TargetPairs_Iterator:
       
    def __init__(self, date, subdir='daily', trunk='spectra', verbose=False):
        self.date = date
        self.subdir = subdir
        self.trunk = trunk
        self.verbose = verbose
        
        self.it0 = None    # Date_Spectra_Iterator
        self.it1 = None    # Spectra_Subspectra_Iterator
        self.it2 = None    # Spectra_Pairs_Iterator
        
#         self.spectra = None
#         self.fname = None

#         self.subspectra = None
        
    def __iter__(self):
        
        self.it0 = None    # Date_Spectra_Iterator
        self.it1 = None    # Spectra_Subspectra_Iterator
        self.it2 = None    # Spectra_Pairs_Iterator
        
#         self.spectra = None
#         self.fname = None
        return self
    
    def __next__(self):
        
        # the first
        
        if self.it0 is None:
            try:
                self.it0 = Date_Spectra_Iterator(self.date,subdir=self.subdir, trunk=self.trunk, verbose=self.verbose)
                self.spectra, self.fname  = self.it0.__next__()
                self.it1 = Spectra_Subspectra_Iterator(self.spectra, verbose=False)
#                 self.subspectra = self.it1.__next__()
                self.it2 = Spectra_Pairs_Iterator(self.it1.__next__(), verbose=False)
                ans = self.it2.__next__()
            except:
                raise StopIteration
        
        while True:

            try:
                ans = self.it2.__next__()
                break
            except StopIteration:
                try:
                    self.it2 = Spectra_Pairs_Iterator(self.it1.__next__(), verbose=False)
                except StopIteration:
                    try:
                        self.spectra, self.fname  = self.it0.__next__()
                        self.it1 = Spectra_Subspectra_Iterator(self.spectra, verbose=False)
                        self.it2 = Spectra_Pairs_Iterator(self.it1.__next__(), verbose = False)
                    except StopIteration:
                        raise StopIteration

        return ans


    @staticmethod
    def test():
        it = Date_TargetPairs_Iterator("20201223",verbose=True) 
        for a in it:
            pass

"""

Given a date, return all pairs of spectra from that date and preceeding dates

"""
class Date_SpectraPairs_Iterator:
       
    def __init__(self, date, subdir='daily', trunk='coadd', verbose =False):
        self.date = date
        self.subdir = subdir
        self.trunk = trunk
        self.verbose = verbose
        
        # less I/O if the latter two were swapped
        self.it0 = None     # Date_Tile_Iterator
        self.it1 = None     # TileDate_PreDate_Iterator 
        self.it2 = None     #  np.nditer(fs_utils.panels)
            
        #current store
        
        self.tile = None
        self.pdate = None
        self.panel = None
        

    def __iter__(self):
        self.it0 = None     # Date_Tile_Iterator
        self.it1 = None     # TileDate_PreDate_Iterator 
        self.it2 = None     #  np.nditer(fs_utils.panels)        
        return self
    
    def __next__(self):
        
        filename = None
        filename2 = None
        
        #handle the first case
        if self.it0 is None:
            try:
                self.it0 =  Date_Tile_Iterator(self.date, subdir=self.subdir)
                self.tile = self.it0.__next__()
                self.it1 =  TileDate_PreDate_Iterator(self.tile, self.date, subdir=self.subdir)
                self.pdate = self.it1.__next__()
                self.it2 = fs_utils.panels.flat
                self.panel=self.it2.__next__()
                filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                filename2 = fs_utils.fitsfile(self.tile, self.pdate, self.panel, subdir=self.subdir,trunk=self.trunk)
            except StopIteration:
                raise StopIteration
                
        while filename is None or filename2 is None:
            try:
                self.panel = self.it2.__next__()                
            except StopIteration:
                try:
                    self.pdate=self.it1.__next__()
                    self.it2 = fs_utils.panels.flat
                    self.panel=self.it2.__next__()
                except StopIteration:
                    try:
                        self.tile = self.it0.__next__()
                        self.it1 =  TileDate_PreDate_Iterator(self.tile, self.date, subdir=self.subdir)
                        self.pdate=self.it1.__next__()
                        self.it2 = fs_utils.panels.flat
                        self.panel=self.it2.__next__()
                    except StopIteration:
                        raise StopIteration

            filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
            filename2 = fs_utils.fitsfile(self.tile, self.pdate, self.panel, subdir=self.subdir,trunk=self.trunk)
            
        if self.verbose:
            print("Iterator: Tile {}, Panel {}, Date {}, Date 2 {}".format(self.tile, self.panel,self.date,self.pdate))
            
        return (read_spectra(filename) , read_spectra(filename2))

    # staticmethod
    def test():
        it = Date_SpectraPairs_Iterator("20201223",verbose=True)
        for a in it:
            pass
            
"""

Given a date, return all pairs of spectra from that date and preceeding dates

"""
class TileDate_SpectraPairs_Iterator:
       
    def __init__(self, tile, date, subdir='daily', trunk='coadd', verbose =False):
        
#         self.list = tdlist
        self.list = [[x,y] for x,y in zip(tile,date)]
        self.subdir = subdir
        self.trunk = trunk
        self.verbose = verbose
        
        # less I/O if the latter two were swapped
        self.it0 = None
        self.it1 = None     # TileDate_PreDate_Iterator 
        self.it2 = None     #  np.nditer(fs_utils.panels)
            
        #current store
        self.tile = None
        self.date = None        
        self.pdate = None
        self.panel = None
        

    def __iter__(self):
        self.it1 = None     # TileDate_PreDate_Iterator 
        self.it2 = None     #  np.nditer(fs_utils.panels)        
        return self
    
    def __next__(self):
        
        filename = None
        filename2 = None
        
        #handle the first case
        if self.it0 is None:
            try:
                self.it0 = iter(self.list)
                (self.tile,self.date)= self.it0.__next__()
                self.it1 =  TileDate_PreDate_Iterator(self.tile, self.date, subdir=self.subdir)
                self.pdate = self.it1.__next__()
                self.it2 = fs_utils.panels.flat
                self.panel=self.it2.__next__()
                filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                filename2 = fs_utils.fitsfile(self.tile, self.pdate, self.panel, subdir=self.subdir,trunk=self.trunk)
            except StopIteration:
                raise StopIteration
                
        while filename is None or filename2 is None:
            try:
                self.panel = self.it2.__next__()                
            except StopIteration:
                try:
                    self.pdate=self.it1.__next__()
                    self.it2 = fs_utils.panels.flat
                    self.panel=self.it2.__next__()
                except StopIteration:
                    try:
                        (self.tile,self.date)= self.it0.__next__()
                        self.it1 =  TileDate_PreDate_Iterator(self.tile, self.date, subdir=self.subdir)
                        self.pdate = self.it1.__next__()
                        self.it2 = fs_utils.panels.flat
                        self.panel=self.it2.__next__()
                    except:
                        raise StopIteration

            filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
            filename2 = fs_utils.fitsfile(self.tile, self.pdate, self.panel, subdir=self.subdir,trunk=self.trunk)
            
        if self.verbose:
            print("Iterator: Tile {}, Panel {}, Date {}, Date 2 {}".format(self.tile, self.panel,self.date,self.pdate))
            
        return (read_spectra(filename) , read_spectra(filename2))

    # staticmethod
    def test():
        it = TileDate_SpectraPairs_Iterator(["80613","80622"],["20201223","20201223"],verbose=True)
        for a in it:
            pass

"""

Given a date, return all pairs of spectra from that date and preceeding dates

"""
class TileDate_TargetPairs_Iterator:
       
    def __init__(self, tile, date, subdir='daily', trunk='spectra', verbose =False):
        
#         self.list = tdlist
        self.list = [[x,y] for x,y in zip(tile,date)]
        self.subdir = subdir
        self.trunk = trunk
        self.verbose = verbose

        
        self.it0 = None    # 
        self.it1 = None    # Spectra_Subspectra_Iterator
        self.it2 = None    # Spectra_Pairs_Iterator
        self.it3 = None
        
    def __iter__(self):
        
        self.it0 = None    # Date_Spectra_Iterator
        self.it1 = None    # Spectra_Subspectra_Iterator
        self.it2 = None    # Spectra_Pairs_Iterator
        self.it3 = None
        return self
    
    def __next__(self):
        
        # the first
        
        if self.it0 is None:
            try:
                self.it0 = iter(self.list)
                (self.tile,self.date)= self.it0.__next__()
                self.it1 = fs_utils.panels.flat
                self.panel = self.it1.__next__()
                filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                self.spectra = read_spectra(filename)
                self.it2 = Spectra_Subspectra_Iterator(self.spectra, verbose=True)
                self.it3 = Spectra_Pairs_Iterator(self.it2.__next__(), verbose=True)
                ans = self.it3.__next__()
            except:
                raise StopIteration
        
        while True:

            try:
                ans=self.it3.__next__()
                break
            except StopIteration:
                try:
                    self.it3 = Spectra_Pairs_Iterator(self.it2.__next__(), verbose=True)
                except StopIteration:
                    try:
                        self.panel = self.it1.__next__()
                        filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                        self.spectra = read_spectra(filename)
                        self.it2 = Spectra_Subspectra_Iterator(self.spectra, verbose=True)
                        self.it3 = Spectra_Pairs_Iterator(self.it2.__next__(), verbose=True)
                    except StopIteration:
                        try:
                            (self.tile,self.date)= self.it0.__next__()
                            self.it1 = fs_utils.panels.flat
                            self.panel = self.it1.__next__()
                            filename = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                            self.spectra = read_spectra(filename)
                            self.it2 = Spectra_Subspectra_Iterator(self.spectra, verbose=True)
                            self.it3 = Spectra_Pairs_Iterator(self.it2.__next__(), verbose=True)
                        except StopIteration:
                            raise StopIteration
        if self.verbose:
            print(self.tile, self.date)
        return ans

    # staticmethod
    def test():
        it = TileDate_TargetPairs_Iterator(["80613","80622"],["20201223","20201223"],verbose=True)
        for a in it:
            pass        
# class PairCoadds:

#     def __init__(self, tile, subdir='andes',dates=None):
#         self.tile = tile
#         self.subdir=subdir
#         if dates is None:
#             self.dates = fs_utils.tileToDates(tile,subdir=subdir)
#         else:
#             self.dates = dates
#         if (len(self.dates) <=1):
#             self.a=1
#         else:
#             self.a = 0
#         self.b = 1
#         self.panel = 0
        
#     def __iter__(self):
#         if (len(self.dates) <=1):
#             self.a=1
#         else:
#             self.a = 0
#         self.b = 1
#         self.panel = 0
#         return self
    
#     def __next_index__(self):
#         if self.panel+1 < len(fs_utils.panels):
#             self.panel +=1
#         else:
#             self.panel=0
#             if self.b +1 < len(self.dates):
#                 self.b+=1
#             else:
#                 self.a+=1
#                 self.b = self.a+1

#     def __next__(self):
           
#         while self.a+1 < len(self.dates):
#             df0 = fs_utils.fitsfile(self.tile,self.dates[self.a],fs_utils.panels[self.panel])
#             df1 = fs_utils.fitsfile(self.tile,self.dates[self.b],fs_utils.panels[self.panel])

#             if (df0 is not None and df1 is not None):
#                 ans =  (read_spectra(df0),read_spectra(df1))
#                 self.__next_index__()
#                 break
        
#             self.__next_index__()

#         else:
#             raise StopIteration
        
#         return ans
    
#     #staticmethod
#     def test():
#         tile = "70006"
#         for c in PairCoadds(tile):
#             print (c)
            
# class DailyPairs:

#     def __init__(self, date, subdir='daily'):
#         self.date = date
#         self.tiles = fs_utils.dateToTiles(date, subdir='daily')
#         self.subdir=subdir
        
#         self.spec = None
#         self.fibermap = None
#         self.tile_index = 0
#         self.panel = 0
        
#     def __iter__(self):
#         self.panel = 0
#         self.tile_index = 0
#         self.df = None
#         return self

#     def __next__(self):
        
#         # if need a new file
#         if self.panel == 0:
#             df = fs_utils.fitsfile(self.tiles[self.tile_index],self.date,fs_utils.panels[self.panel],trunk='spectra')
#             print(df)
#             self.spec = read_spectra(df)
#             self.fibermap = self.spec.fibermap
            
#         print(self.fibermap)
#         wef
            
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
            
# #         while self.a+1 < len(self.dates):
# #             df0 = fs_utils.fitsfile(tile,self.dates[self.a],fs_utils.panels[self.panel])
# #             df1 = fs_utils.fitsfile(tile,self.dates[self.b],fs_utils.panels[self.panel])

# #             if (df0 is not None and df1 is not None):
# #                 ans = (read_spectra(df0), read_spectra(df1))
# # #                 ans =  (df0,df1)
# #                 self.__next_index__()
# #                 break
        
# #             self.__next_index__()

# #         else:
# #             raise StopIteration
        
#         print(self.tile_index,self.panel)
#         return ans

#     #staticmethod
#     def test():
#         for c in PairSingle("70006","20200305"):
#             print (c)

# # PairSingle.test()

# class TileDate:
    
#     def __init__(self, tile, date, subdir='andes'):
#         dirname = os.path.join(redux,subdir,'tiles',tile,date)
#         if not os.path.isdir(dirname):
#             print('{} does not exist.'.format(dirname))    
#         self.tile = tile
#         self.date = date
#         self.subdir = subdir
#         cafiles =  glob(os.path.join(dirname,'cframe-??-*.fits'))
#         exposures=[]

#         for cafile in cafiles:
#             index0 = cafile.find('cframe-')
#             index1 = cafile.find('.fits',index0+10)
#             exposures.append(cafile[index0+10:index1])

#         exposures = np.unique(exposures)
#         print(exposures)
        
#     def __iter__(self):
#         self.a = 1
#         return self

#     def __next__(self):
#         x = self.a
#         self.a += 1
#         return x
    
# # tile = "70006"
# # date = "20200305"

