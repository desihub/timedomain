import requests
import time
from astropy.time import Time


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
    def api(method, endpoint, data=None):
        response = requests.request(method, endpoint, json=data, headers=SkyPortal.headers)
        return response
    
    
    @staticmethod
    def postCandidate(index,fibermap,filt):
    	filterid = SkyPortal.filter_id(filt)

        data = {
          "ra": fibermap['TARGET_RA'].data[index],
          "dec": fibermap['TARGET_DEC'].data[index],
          "id": "DESI{}".format(fibermap['TARGETID'].data[index].astype('str')),
          "origin": "DESI",
          "filter_ids": [
            filterid
          ],
          "passing_alert_id": 0,
          "passed_at": time.strftime("20%y-%m-%d",time.gmtime())
        }

        response = SkyPortal.api('POST', '{}/api/candidates'.format(SkyPortal.url),data=data)

        print(f'HTTP code: {response.status_code}, {response.reason}')
        if response.status_code in (200, 400):
            print(f'JSON response: {response.json()}')
    
    

    @staticmethod            
    def postSpectra(target_id, spectra_in):
        spectra = spectra_in.select(targets=[int(target_id)])
        fibermap = spectra.fibermap

        objID = 'DESI{}'.format(target_id)
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
                response = SkyPortal.api('POST', '{}/api/spectrum'.format(SkyPortal.url),data=data)
                print(f'HTTP code: {response.status_code}, {response.reason}')
                if response.status_code in (200, 400):
                    print(f'JSON response: {response.json()}')
        

#     @staticmethod
#     def postMayallTelescope():
#         data = {
#           "name": SkyPortal.telescope_name,
#           "nickname": "KPNO 4-m",
#           "lat": 31.9634,
#           "lon": -111.6,
#           "elevation": 2120,
#           "diameter": 4,
#           "robotic": False,
#           "fixed_location": True
#         }
        
#         print(data)
#         response = SkyPortal.api('POST', '{}/api/telescope'.format(SkyPortal.url),data=data)
#         print(f'HTTP code: {response.status_code}, {response.reason}')
#         if response.status_code in (200, 400):
#             print(f'JSON response: {response.json()}')
            
#     @staticmethod
#     def postDESIInstrument():
#         data = {
#           "name": SkyPortal.instrument_name,
#           "type": "spectrograph",
#           "band": "optical",
#           "telescope_id": SkyPortal.telescope_id()
#         }
        
#         print(data)
#         response = SkyPortal.api('POST', '{}/api/instrument'.format(SkyPortal.url),data=data)
#         print(f'HTTP code: {response.status_code}, {response.reason}')
#         if response.status_code in (200, 400):
#             print(f'JSON response: {response.json()}')
            
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