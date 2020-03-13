import argparse
from math import *

parser = argparse.ArgumentParser(description='Client API to query ANTARES alert database')
parser.add_argument("--ra", default=None, type=float, help="RA (J2000)")
parser.add_argument("--dec", default=None, type=float, help="Dec (J2000)")
parser.add_argument("--tobs", default=None, type=float, help="Time of observation (MJD)")

args = parser.parse_args()
ra = args.ra
dec = args.dec
tobs = args.tobs
if ra==None or dec==None or tobs==None:
    print ('Usage: query_alert.py --ra --dec --tobs')
    exit()
    

dra = 1.5*cos(dec/180.*pi)
ddec = 1.5
#dt=3
#print (tobs-dt)

query = {
    "query": {
        "bool": {
            "must": [
                {
                    "match":{
                        "properties.snfilter_last_proc_status": "*DESI*"
                    }
                },
                {
                    "range": {
                        "ra": {
                            "gte": ra-dra,
                            "lte": ra+dra,
                        }
                    }
                },
                {
                    "range": {
                        "dec": {
                            "gte": dec-ddec,
                            "lte": dec+ddec,
                        }
                    }
                },
                {
                    "range": {
                        "mjd": {
                            "gte": tobs,
                        }
                    }
                }
             ]
        }
    }
}

from antares_client.search import search
result_set = search(query)

from antares_client.search import download
result_set = download(query, "results.csv", output_format="csv", decompress=True)
