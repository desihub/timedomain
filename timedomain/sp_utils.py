import requests
import time
from astropy.time import Time
import json
from desispec.coaddition import coadd_cameras

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
    def postCandidate(index,fibermap,filt, data_override=None):
        
        # Does the candidate exist?
        response = SkyPortal.api('GET', '{}/api/candidates/DESI{}'.format(SkyPortal.url,fibermap['TARGETID'].data[index]))
        
        # if it doesn't exist save
        if response.status_code != 200:
        
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
                "origin": "DESI",
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
                "data": altdata
                }            
            response = SkyPortal.api('POST', '{}/api/annotation'.format(SkyPortal.url),data=data)

            print(f'HTTP code: {response.status_code}, {response.reason}')
            if response.status_code == 400:
                print(f'JSON response: {response.json()}')        

        else:
            # clear out the spectra for now
            response = SkyPortal.api('GET', '{}/api/sources/DESI{}/spectra'.format(SkyPortal.url,fibermap['TARGETID'].data[index]))
            for s in response.json()['data']['spectra']:
                response = SkyPortal.api('DELETE','{}/api/spectrum/{}'.format(SkyPortal.url,s['id']))
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
                  "origin": "DESI",
                  "group_ids": [23, 24],
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