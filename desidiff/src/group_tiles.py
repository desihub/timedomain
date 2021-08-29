import glob
from astropy.table import Table, vstack, unique, SortedArray
import numpy

# arguments:
#  yyyymmdd of the night of interest
# outputs:
#  tile_petal: subsets of tile/petals whose RA/DEC targets are not contained in other tile/petals from that night
#  group_tid: the set of targetid's associated with each RA/DEC
#  group_tp:  the set of tile/petals associated with each RA/DEC

def getMatchedTileid(yyyymmdd):
    
    # read in and store in one place all the fibermap information in the spectra files
    dats=[]
    for filename in glob.glob(f"/global/project/projectdirs/desi/spectro/redux/daily/tiles/cumulative/*/{yyyymmdd}/spectra-*.fits"):
        t = Table.read(filename, format='fits',hdu=1, memmap=True) \
            [['TARGETID','TARGET_RA','TARGET_DEC','TILEID','OBJTYPE','PETAL_LOC','FIBERSTATUS','NIGHT']]
        
#         ##### DAVE'S ADDITION ##############
#         targetcols = [i for i in t.colnames if i[-7:] =='_TARGET']
#         nonzerocheck = [True in k for k in [[j != 0 for j in t[targetcols][i]] for i in range(len(t))]]
#         #a really ugly line, basically generates a list of bools, 
#         #True if there is at least one nonzero element in all columns ending in _TARGET
#         t.remove_rows([i for i, val in enumerate(nonzerocheck) if not val])
#         #This gets the index of all False values from the previous list and removes those rows
#         t = t['TARGETID','TARGET_RA','TARGET_DEC','TILEID','OBJTYPE','PETAL_LOC','FIBERSTATUS']
#         ######## END DAVE'S ADDITION ############
        
        t=t[numpy.logical_and(t['OBJTYPE']=='TGT', t['FIBERSTATUS']==0)]
        dats.append(t)
    dats=vstack(dats, join_type='inner',metadata_conflicts='silent')

    # group all the observations by TARGET_RA and TARGET_DEC
    # note that this is more reliable than grouping by TARGETID as TARGETID is NOT a unique identifier of RA and DEC
    dats_group = dats.group_by(['TARGET_RA','TARGET_DEC'])
    
    # for each group make a list of all targetid's with that ra/dec
    group_tid=[]
    group_night=[]
    # for each group make a tuple containing all tileid/petal combinations that have that ra/dec
    group_tp=[]    
    for g in dats_group.groups:
        group_tid.append(numpy.unique(g['TARGETID'].data))
        group_night.append(numpy.unique(g['NIGHT'].data))
#     for g in dats_group.groups:
        gu = unique(g,keys=['TILEID'])
        dum=[]
        for a,b in zip(gu['TILEID'],gu['PETAL_LOC']):
            dum.append((a,b))
        group_tp.append(tuple(dum))

    # compress things down to the unique tile/petal combinations
    ans = list(set(group_tp))
    
    # union sets that have intersecting tile/petal combinations
    for i in range(len(ans)-1,0,-1):
        for j in range(i-1,-1,-1):
            if len(set(ans[i]) & set(ans[j])) !=0:
                ans[j]=ans[i]+ans[j]
                del ans[i]
                break

    # unique of each set
    for i in range(len(ans)):
        ans[i]=list(set(ans[i]))
        
    tile_petal=ans
        
    # for each set the TARGET_RA TARGET_DEC associated
        
    return(tile_petal, group_tid, group_tp, group_night)