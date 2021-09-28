import requests
import time
from astropy.time import Time
import json
from desispec.coaddition import coadd_cameras
import speclite
import numpy
import json
from desiutil.log import get_logger
log = get_logger()

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class ResponseError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, response):
        self.response = response
        
class UniqueViolationError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, response):
        self.response = response
        
        
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
#             if response.status_code in (200, 400):
#                 print(f'JSON response: {response.json()}')
            if response.status_code ==400:
                raise Exception('Something wrong')
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
        try:
            response = requests.request(method, endpoint, json=data, headers=SkyPortal.headers)
        except:
            raise
        if response.status_code == 400:
            if "psycopg2.errors.UniqueViolation" in response.json()['message']:
                raise UniqueViolationError(response)
            else:
                raise ResponseError(response)
        return response

    @staticmethod
    def nukeCandidate(targetid, filt):
        # Remove annotations related to this object
        response = SkyPortal.api('GET', '{}/api/sources/DESI{}/annotations'.format(SkyPortal.url,targetid))
        for an in response.json()['data']:
            r2 = SkyPortal.api('DELETE', '{}/api/annotation/{}'.format(SkyPortal.url,an['id']))
            log.info(f'Deleting Annotation JSON response: {r2.json()}')

        # Remove photometry related to this object
        response = SkyPortal.api('GET', '{}/api/sources/DESI{}/photometry'.format(SkyPortal.url,targetid))
        for an in response.json()['data']:
            r2 = SkyPortal.api('DELETE', '{}/api/photometry/{}'.format(SkyPortal.url,an['id']))
            log.info(f'Deleting Photometry JSON response: {r2.json()}')

        # Remove spectra related to this object
        response = SkyPortal.api('GET', '{}/api/sources/DESI{}/spectra'.format(SkyPortal.url,targetid))
        if response.status_code == 200:
            for an in response.json()['data']['spectra']:
                r2 = SkyPortal.api('DELETE', '{}/api/spectrum/{}'.format(SkyPortal.url,an['id']))
                log.info(f'Deleting Spectra JSON response: {r2.json()}')
                             
        filterid = SkyPortal.filter_id(filt)        
        response = SkyPortal.api('DELETE', '{}/api/candidates/DESI{}/{}'.format(SkyPortal.url,targetid,filterid))
        log.info(f'Deleting Candidat JSON response: {response.json()}')
        
    
    @staticmethod
    def postCandidate(index,fibermap,filt, data_override=None):

        filterid = SkyPortal.filter_id(filt)

        # information we want to save from fibermap
        fibermap_keys=['FIBER','TILEID','EXPID','PETAL_LOC','MJD','LAST_MJD']
        fiber_dict = dict()
        for key in fibermap_keys:
            if key in fibermap.keys():
                fiber_dict[key] = fibermap[key].data[index].astype('str')

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

        try:
            response = SkyPortal.api('POST', '{}/api/candidates'.format(SkyPortal.url),data=data)
        except:
            raise
            
#         log.info(f'HTTP code: {response.status_code}, {response.reason}')
#         if response.status_code == 400:
#             log.warning(print(f'JSON response: {response.json()}'))

    def postAnnotation(index,fibermap, data_override=None):

        # information we want to save from fibermap
        fibermap_keys=['FIBER','TILEID','EXPID','PETAL_LOC','MJD','LAST_MJD']
        fiber_dict = dict()
        for key in fibermap_keys:
            if key in fibermap.keys():
                fiber_dict[key] = fibermap[key].data[index].astype('str')

        altdata = {'fibermap': fiber_dict}
        
        data = {
            "obj_id": "DESI{}".format(fibermap['TARGETID'].data[index].astype('str')),
            "origin": 'postCandidate',
            "data": altdata['fibermap']
            }
        
        if data_override is not None:
            for k,v in data_override.items():
                if isinstance(v,dict) and k in data:
                    data[k].update(v)
                    data_override[k]=data[k]

            data.update(data_override)

        try:
            response = SkyPortal.api('POST', '{}/api/annotation'.format(SkyPortal.url),data=data)
        except:
            raise
            
#         log.info(f'HTTP code: {response.status_code}, {response.reason}')
#         if response.status_code == 400:
#             log.warning(print(f'JSON response: {response.json()}'))

    

    @staticmethod            
    def postSpectra(target_id, spectra_in, data_override=None,coadd_camera=False):
            
        # fix to an upstream bug in spectra
        keep = numpy.equal(spectra_in.fibermap['TARGETID'].data, int(target_id))
        spectra = spectra_in[keep]

        # At this point there should only be one Spectrum
        # This should be the case of everything derived from "coadd" or from
        # SpectraPairsIterator

        if spectra.num_spectra() !=1:
            log.error(target_id, int(target_id) in spectra_in.fibermap['TARGETID'].data)
            log.error(spectra.num_spectra())
            raise IndexError
        
#         spectra = spectra_in.select(targets=[int(target_id)])

        # coadd_cameras does not have an MJD.  Rely on getting this from the original.
        fibermap = spectra.fibermap
        
        if 'MJD' in fibermap.keys():
            mjd = fibermap['MJD'][0]
        elif 'LAST_MJD' in fibermap.keys():
            mjd = fibermap['LAST_MJD'][0]
        else:
            mjd=-999
            log.error('we are screwed')
        t = Time(mjd, format='mjd')
            
        # combines the arms into one spectrum
        if coadd_camera:
            spectra = coadd_cameras(spectra)

        objID = 'DESI{}'.format(target_id)
        
        # this for statement should be superfluous as there is only one index
#         for index, mjd in enumerate(fibermap.iterrows('MJD')):
        index=0

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

            try:
                response = SkyPortal.api('POST', '{}/api/spectrum'.format(SkyPortal.url),data=data)
            except:
                raise
                
#             log.info(f'HTTP code: {response.status_code}, {response.reason}')
#             if response.status_code == 400:
#                 log.error(f'JSON response: {response.json()}')
 
    @staticmethod            
    def putPhotometry(target_id, spectra_in, data_override=None,coadd_camera=False):
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
        
        if 'MJD' in fibermap.keys():
            mjd = fibermap['MJD'][0]
        elif 'LAST_MJD' in fibermap.keys():
            mjd = fibermap['LAST_MJD'][0]
        else:
            log.error('we are screwed')
            
                
        # combines the arms into one spectrum
        if coadd_camera:
            spectra = coadd_cameras(spectra)
            if 'brz' not in spectra.bands:
                log.warning('postPhotometry failed because brz not available but carrying on')
                return
            
        objID = 'DESI{}'.format(target_id)
        
        # this for statement should be superfluous as there is only one index
#         for index, mjd in enumerate(fibermap.iterrows('MJD')):
        index=0
        fluxes = []
        fluxerrs = []
        for b in bands:
            ok = numpy.logical_and(spectra.mask['brz'][index,:] ==0,spectra.ivar[spectra.bands[0]][index,:])
            padspec = b.pad_spectrum(spectra.flux[spectra.bands[0]][index,ok],spectra.wave[spectra.bands[0]][ok])
            fluxes.append(b.get_ab_maggies(padspec[0],padspec[1])*1e-17)
            padspec = b.pad_spectrum(spectra.ivar[spectra.bands[0]][index,ok],spectra.wave[spectra.bands[0]][ok])
            fluxerrs.append(numpy.sqrt(b.get_ab_maggies(1/padspec[0],padspec[1]))*1e-17)
            
        filters = ['sdssg','sdssr','sdssz']

        w=numpy.where(numpy.isfinite(fluxes))[0]
        
        if len(w)!=0:

            data = {
                "mjd": str(mjd),
                "obj_id": objID,
                "instrument_id": SkyPortal.instrument_id(),
                "origin": "DESI",
                "magsys": 'ab',
                "filter": filters[w[0]],
                "group_ids": [SkyPortal.group_id('DESI') ],
                "flux":     fluxes[w[0]],
                "fluxerr": fluxerrs[w[0]],
                "zp": 0
                }

            if data_override is not None:
                for k,v in data_override.items():
                    if isinstance(v,dict) and k in data:
                        data[k].update(v)
                        data_override[k]=data[k]
                data.update(data_override)

            try:
                response = SkyPortal.api('PUT', '{}/api/photometry'.format(SkyPortal.url),data=data)
            except ResponseError as err:
                raise
#             log.info(f'HTTP code: {response.status_code}, {response.reason}')                                                                         
#             if response.status_code == 400:
#                 log.warning(f'JSON response: {response.json()}')

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

    @staticmethod            
    def deleteSpectra(target_id, altdata):
        try:
            response=SkyPortal.api('GET',  '{}/api/sources/DESI{}/spectra'.format(SkyPortal.url,target_id))
        except:
            raise

        spec = response.json()['data']['spectra']
        for s in spec:
            if (altdata == s['altdata']):
                try:
                    response=SkyPortal.api('DELETE',  '{}/api/spectrum/{}'.format(SkyPortal.url,s['id']))
                except:
                    raise
        
# SkyPortal.instrument_id()

# SkyPortal.test()

#     @staticmethod
#     def test2():
#         print(SkyPortal.filter_id("DESI Difference CV"))
#         print(SkyPortal.filter_id("AMPEL.HU_RANDOM"))
#         print(SkyPortal.filter_id("DESIDIFF_CV_daily"))
#         print(SkyPortal.filter_id("xyz"))
        
# SkyPortal.test2()