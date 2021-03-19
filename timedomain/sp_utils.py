import requests
import time
from astropy.time import Time
import json
from desispec.coaddition import coadd_cameras
import speclite
import numpy

class SkyPortal:

    #skyportal location and 
    url = "http://desi2.lbl.gov:5000"
    
    #skyportal token im a secret file
#     secret_file = "/global/cscratch1/sd/akim/secrets/desi_sp.txt"
    secret_file = "/global/cfs/cdirs/desi/science/td/secrets/desi_sp.txt"
    with open(secret_file, 'r') as file:
        token = file.read().replace('\n', '')
#     token='5e126952-11b8-4f37-af2c-2e0056635c09'
    headers = {'Authorization': f'token {token}'}
    
    instrument_name = 'DESI'
    telescope_name ='Kitt Peak Mayall 4-m Telescope'
    
    inst_id = None
    tel_id = None
    filt_id = {}
    grp_id={}
    strm_id = {}

    @staticmethod
    def instrument_id():
        
        if SkyPortal.inst_id is None:
            response = SkyPortal.api('GET', '{}/api/instrument'.format(SkyPortal.url))
            data = response.json()['data']
            theOne = list(filter(lambda datum: datum['name']==SkyPortal.instrument_name,data))
            if len(theOne) !=0:
                SkyPortal.inst_id = theOne[0]['id']
            else:
                raise NameError('Instrument Not Defined')
        
        return SkyPortal.inst_id

    @staticmethod
    def telescope_id():
        
        if SkyPortal.tel_id is None:       
            response = SkyPortal.api('GET', '{}/api/telescope'.format(SkyPortal.url))
            data = response.json()['data']
            theOne = list(filter(lambda datum: datum['name']==SkyPortal.telescope_name,data))
            if len(theOne) !=0:
                SkyPortal.tel_id = theOne[0]['id']
            else:
                raise NameError('Telescope Not Defined')
        
        return SkyPortal.tel_id 
    
    @staticmethod
    def filter_id(filter_name):
        if filter_name in SkyPortal.filt_id:
            return SkyPortal.filt_id[filter_name]
        else:
            response = SkyPortal.api('GET', '{}/api/filters'.format(SkyPortal.url))
            data = response.json()['data']
            theOne = list(filter(lambda datum: datum['name']==filter_name,data))

            if len(theOne) !=0:
                SkyPortal.filt_id[filter_name] = theOne[0]['id']
            else:
                raise NameError('Filter Not Defined')
            return SkyPortal.filt_id[filter_name]
            
    @staticmethod
    def group_id(group_name):
        if group_name in SkyPortal.grp_id:
            return SkyPortal.grp_id[group_name]
        else:
            response = SkyPortal.api('GET', '{}/api/groups'.format(SkyPortal.url))
            data = response.json()['data']['all_groups']
            theOne = list(filter(lambda datum: datum['name']==group_name,data))
            if len(theOne) !=0:
                SkyPortal.grp_id[group_name] = theOne[0]['id']
            else:
                raise NameError('Group Not Defined')
            return SkyPortal.grp_id[group_name]           
        
    @staticmethod
    def stream_id(name):
        if name in SkyPortal.strm_id:
            return SkyPortal.strm_id[name]
        else:
            response = SkyPortal.api('GET', '{}/api/streams'.format(SkyPortal.url))
            data = response.json()['data']
            theOne = list(filter(lambda datum: datum['name']==name,data))
            if len(theOne) !=0:
                SkyPortal.strm_id[name] = theOne[0]['id']
            else:
                raise NameError('Stream Not Defined')
            return SkyPortal.strm_id[name]   
        
    @staticmethod
    def api(method, endpoint, data=None):
        response = requests.request(method, endpoint, json=data, headers=SkyPortal.headers)
        return response

    @staticmethod
    def nukeCandidate(targetid, filt):
        
        print(targetid)
        # Remove annotations related to this object
        response = SkyPortal.api('GET', '{}/api/sources/DESI{}/annotations'.format(SkyPortal.url,targetid))
        for an in response.json()['data']:
            r2 = SkyPortal.api('DELETE', '{}/api/annotation/{}'.format(SkyPortal.url,an['id']))
            print(f'JSON response: {r2.json()}')

        # Remove photometry related to this object
        response = SkyPortal.api('GET', '{}/api/sources/DESI{}/photometry'.format(SkyPortal.url,targetid))
        for an in response.json()['data']:
            r2 = SkyPortal.api('DELETE', '{}/api/photometry/{}'.format(SkyPortal.url,an['id']))
            print(f'JSON response: {r2.json()}')

        # Remove photometry related to this object
        response = SkyPortal.api('GET', '{}/api/sources/DESI{}/spectroscopy'.format(SkyPortal.url,targetid))
        for an in response.json()['data']:
            r2 = SkyPortal.api('DELETE', '{}/api/spectroscopy/{}'.format(SkyPortal.url,an['id']))
            print(f'JSON response: {r2.json()}')
                        
            
        filterid = SkyPortal.filter_id(filt)        
        response = SkyPortal.api('DELETE', '{}/api/candidates/DESI{}/{}'.format(SkyPortal.url,targetid,filterid))
        print(f'JSON response: {response.json()}')
        
    
    @staticmethod
    def postCandidate(index,fibermap,filt, data_override=None):
                
        filterid = SkyPortal.filter_id(filt)

        # information we want to save from fibermap
        fibermap_keys=['FIBER','TILEID','EXPID','PETAL_LOC','MJD']
        fiber_dict = dict()
        for key in fibermap_keys:
            if key in fibermap.keys():
                fiber_dict[key] = fibermap[key].data[index].astype('str')
            else:
                print('missing ',key)

        altdata = {'fibermap': fiber_dict}

        data = {
            "ra": fibermap['TARGET_RA'].data[index],
            "dec": fibermap['TARGET_DEC'].data[index],
            "id": "DESI{}".format(fibermap['TARGETID'].data[index].astype('str')),
            "filter_ids": [
            filterid
            ],
            "passing_alert_id": filterid,
            "passed_at": time.strftime("20%y-%m-%d",time.gmtime()),
            "altdata": altdata
            }

        if data_override is not None:
            for k,v in data_override.items():
                if isinstance(v,dict) and k in data:
                    data[k].update(v)
                    data_override[k]=data[k]

            data.update(data_override)

        response = SkyPortal.api('POST', '{}/api/candidates'.format(SkyPortal.url),data=data)

        print(f'HTTP code: {response.status_code}, {response.reason}')
        if response.status_code == 400:
            print(f'JSON response: {response.json()}')

        data = {
            "obj_id": "DESI{}".format(fibermap['TARGETID'].data[index].astype('str')),
            "origin": 'postCandidate',
            "data": altdata['fibermap']
            }
        response = SkyPortal.api('POST', '{}/api/annotation'.format(SkyPortal.url),data=data)
        print(f'HTTP code: {response.status_code}, {response.reason}')
        if response.status_code == 400:
            print(f'JSON response: {response.json()}')

    

    @staticmethod            
    def postSpectra(target_id, spectra_in, data_override=None,coadd_camera=False):
        spectra = spectra_in.select(targets=[int(target_id)])

        # coadd_cameras does not have an MJD.  Rely on getting this from the original.
        fibermap = spectra.fibermap
        
        # At this point there should only be one Spectrum
        # This should be the case of everything derived from "coadd" or from
        # SpectraPairsIterator

        if spectra.num_spectra() >1:
            raise IndexError
                
        # combines the arms into one spectrum
        if coadd_camera:
            spectra = coadd_cameras(spectra)

        objID = 'DESI{}'.format(target_id)
        
        # this for statement should be superfluous as there is only one index
        for index, mjd in enumerate(fibermap.iterrows('MJD')):
            t = Time(mjd[0], format='mjd')
            for band in spectra.bands:
                data = {
                  "obj_id": objID,
                  "group_ids": [SkyPortal.group_id('DESI')],
                  "fluxes": spectra.flux[band][index,:].tolist(),
    #                   "errors": (1./np.sqrt(spectra.ivar[band][index,:])).tolist(),  
                  "wavelengths": spectra.wave[band].tolist(),
                  "observed_at": t.isot,
                  "instrument_id": SkyPortal.instrument_id()
                }                    

                if data_override is not None:
                    for k,v in data_override.items():
                        if isinstance(v,dict) and k in data:
                            data[k].update(v)
                            data_override[k]=data[k]

                    data.update(data_override)

                response = SkyPortal.api('POST', '{}/api/spectrum'.format(SkyPortal.url),data=data)
                print(f'HTTP code: {response.status_code}, {response.reason}')
                if response.status_code == 400:
                    print(f'JSON response: {response.json()}')
 
    @staticmethod            
    def postPhotometry(target_id, spectra_in, data_override=None,coadd_camera=False):
        spectra = spectra_in.select(targets=[int(target_id)])

        # At this point there should only be one Spectrum
        # This should be the case of everything derived from "coadd" or from
        # SpectraPairsIterator
        if spectra.num_spectra() >1:
            raise IndexError

        # The bands of interest
        bands = speclite.filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z')
    
            
        # coadd_cameras does not have an MJD.  Rely on getting this from the original.
        fibermap = spectra.fibermap
                
        # combines the arms into one spectrum
        if coadd_camera:
            spectra = coadd_cameras(spectra)
            if 'brz' not in spectra.bands:
                print('postPhotometry failed because brz not available but carrying on')
                return
            
        objID = 'DESI{}'.format(target_id)
        
        # this for statement should be superfluous as there is only one index
        for index, mjd in enumerate(fibermap.iterrows('MJD')):
            fluxes = []
            for b in bands:
                padspec = b.pad_spectrum(spectra.flux[spectra.bands[0]][index,:],spectra.wave[spectra.bands[0]])
                fluxes.append(b.get_ab_maggies(padspec[0],padspec[1])*1e-17)
                    
            data = {
                "mjd": mjd[0],
                "obj_id": objID,
                "instrument_id": SkyPortal.instrument_id(),
                "origin": "DESI",
                "magsys": 'ab',
                "filter": ['sdssg','sdssr','sdssz'],
                "group_ids": [SkyPortal.group_id('DESI') ],
                "flux":     fluxes,
                "fluxerr": (numpy.array(fluxes)/100).tolist(),
                "zp": 0
                }
            
            if data_override is not None:
                for k,v in data_override.items():
                    if isinstance(v,dict) and k in data:
                        data[k].update(v)
                        data_override[k]=data[k]

                data.update(data_override)

            print(data)
            response = SkyPortal.api('POST', '{}/api/photometry'.format(SkyPortal.url),data=data)
            print(f'HTTP code: {response.status_code}, {response.reason}')
            if response.status_code == 400:
                print(f'JSON response: {response.json()}')
            
    @staticmethod
    def test():

        
        # response = SkyPortal.api('DELETE', '{}/api/candidates/{}'.format(SkyPortal.url,"35191288252861933"))

        # print(f'HTTP code: {response.status_code}, {response.reason}')
        # if response.status_code in (200, 400):
        #     print(f'JSON response: {response.json()}')
            
        filename = fitsfile("70006","20200305", "6", subdir='andes',trunk='coadd') 
        # fibermap = Table.read(filename, 'FIBERMAP')
        # index = np.where(fibermap['TARGETID'].data== 35191288252861933)[0]
        # SkyPortal.postCandidate(index[0],fibermap)
        
        spectra = read_spectra(filename)
#         SkyPortal.postSpectra("35191288252861933",spectra)

        
        
# SkyPortal.instrument_id()

# SkyPortal.test()

#     @staticmethod
#     def test2():
#         print(SkyPortal.filter_id("DESI Difference CV"))
#         print(SkyPortal.filter_id("AMPEL.HU_RANDOM"))
#         print(SkyPortal.filter_id("DESIDIFF_CV_daily"))
#         print(SkyPortal.filter_id("xyz"))
        
# SkyPortal.test2()