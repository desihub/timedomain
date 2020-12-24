from desispec.io import read_spectra, write_spectra
from desispec.spectra import Spectra

from . import fs_utils

##

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
            self.index = self.index+1
            return ans
        else:
            raise StopIteration
            
    def hasNext(self):
        if self.index < len(self.tiles):
            return True
        else:
            return False
            
    
class Panel_Iterator:
    def __init__(self):
        self.index = 0
        
    def __iter__(self):
        self.index = 0
        return self
    
    def __next__(self):
        
        if self.index < len(fs_utils.panels):
            ans = fs_utils.panels[self.index]
            self.index = self.index+1
            return ans
        else:
            raise StopIteration
            
    def hasNext(self):
        if self.index < len(fs_utils.panels):
            return True
        else:
            return False
            
            
class Date_Spectra_Iterator:
       
    def __init__(self, date, subdir='andes', trunk='coadd'):
        self.date = date
        self.subdir = subdir
        self.trunk = trunk
        self.date_tile_iterator = Date_Tile_Iterator(self.date, subdir=subdir)
        self.panel_iterator = Panel_Iterator()
        
        self.tile = None
        self.panel = None

    def __iter__(self):
        
        self.date_tile_iterator = Date_Tile_Iterator(self.date, subdir=self.subdir)
        self.panel_iterator = Panel_Iterator()
        self.tile = None
        self.panel = None
        return self
    
    def __next__(self):
        
        while True:
            if self.tile is None:
                self.tile = self.date_tile_iterator.__next__()

            if self.panel_iterator.hasNext():
                self.panel = self.panel_iterator.__next__()
                fname = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                if fname is not None:
                    break

            elif self.date_tile_iterator.hasNext():
                self.tile = self.date_tile_iterator.__next__()
                self.panel_iterator = Panel_Iterator()
                self.panel = self.panel_iterator.__next__()
                fname = fs_utils.fitsfile(self.tile, self.date, self.panel, subdir=self.subdir,trunk=self.trunk)
                if fname is not None:
                    break
            else:
                raise StopIteration
        
        return read_spectra(fname), fname
                    

## Difference Iterators


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

