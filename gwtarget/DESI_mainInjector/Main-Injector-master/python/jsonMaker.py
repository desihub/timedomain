import numpy as np

# Support:
#
# How to use JSONs
#   import json
#   from pprint import pprint
#   json_data=open('json_data')
#   data = json.load(json_data)
#   pprint(data)
#   json_data.close()
# 
# What is format of JSON for CTIO?
#   http://dns.ctio.noao.edu/noao/content/Preparing-Observing-Scripts
#   {"count": 3,
#       "expType": "object",
#       "object": "Sgr",
#       "filter": "g",
#       "RA": 281.0,
#       "dec": -29.866,
#       "expTime": 60
#       },

#
# This code wants a list ra,decs and will write a JSON file
# to cause the Blanco to observe them.
#
def writeJson(ra,dec,id, seqid="none", seqnum=0, seqtot=0,
        exposureList = [90,90,90], 
        filterList = ["i","z","z"],
        tilingList = [1, 5],
        trigger_id = "LIGO/Virgo", trigger_type = "bright", propid='propid',
        skymap="bayestar.fits", jsonFilename="des-gw.json") :
    offsets = tileOffsets()
    fd = open(jsonFilename,"w")
    fd.write("[\n")

# we're going to implement tiling by running Francisco's code
# git clone git@github.com:paztronomer/gw_helper.git
# which duplicates json files with changes-
# here we underfill the slot with just 1 tiling,
# then use gw_helper to duplicate the json with a different
# dither/tile
    tilingList = [tilingList[0],]

    size = ra.size
    nexp = np.size(exposureList)
    ntiles = np.size(tilingList) 
    seqtot= seqtot*nexp
    for i in range(0,size) :
        for k in range(0,ntiles) :
            for j in range(0,nexp) :
                seqnum +=1
                tiling = tilingList[k]
                filter = filterList[j]
                exp = exposureList[j]
                offsets[tiling]
                delRa = offsets[tiling][0]
                delDec = offsets[tiling][1]
                tra = ra[i]
                tdec = dec[i]
                tdec = tdec+delDec
                tra = tra + delRa/np.cos(tdec*2*np.pi/360.)
                if tra < 0 : tra = tra+360.
                if tra > 360. : tra = tra-360.
                #comment = "DESGW: LIGO {} event {}: {} of {}, hex {} tiling {}".format(
                #    trigger_type, seqid, seqnum, seqtot, id[i], 9)
                comment = "{} strategy {} on {}: image {} of {}, filter {}, tiling {} in {}".format(
                    trigger_id, trigger_type, skymap, j+1, nexp, filter, tiling, tilingList)
                object = trigger_id
    
                fd.write("{")
                fd.write(" \"expType\" : \"object\",\n")
                fd.write("  \"object\" : \"{}\",\n".format(object))
                #fd.write("  \"seqid\" : \"{}\",\n".format(seqid))
                #fd.write("  \"seqnum\" : \"{:d}\",\n".format(int(seqnum)))
                #fd.write("  \"seqtot\" : \"{:d}\",\n".format(int(seqtot)))
                fd.write("  \"expTime\" : {:d},\n".format(int(exp)))
                fd.write("  \"wait\" : \"False\",\n")
                #fd.write("  \"count\" : \"1\",\n")
                #fd.write("  \"note\" : \"Added to queue from desgw json file, not obstac\",\n")
                fd.write("  \"filter\" : \"{}\",\n".format(filter))
                fd.write("  \"program\" : \"des gw\",\n")
                fd.write("  \"RA\" : {:.6f},\n".format(tra))
                fd.write("  \"dec\" : {:.5f},\n".format(tdec))
                fd.write("  \"propid\" : \"{}\",\n".format(propid))
                fd.write("  \"comment\" : \"{}\"\n".format(comment)) 
                # note lack of comma for end
                fd.write("}")
                if (i == size-1) and ( j == nexp-1) and ( k == ntiles-1) :
                    pass
                else :
                    fd.write(",")
                fd.write("\n")
    fd.write("]\n")
    fd.close()

# production offsets are in DES docdb 7269
# production offsets are one based, not zero based:
def tileOffsets() :
    offsets = dict()
# blessed offsets 0
    offsets[0] = 0.000, 0.000
# we'll use productions 9,10 (as if we were doing DES year 5)
# production offsets: docdb 7269, offsets_cross.txt
    offsets[1] = 0.000, 0.000

# real
#    offsets[2] = -0.76668, 0.473424 
#    offsets[3] = -0.543065, -0.828492
#    offsets[4] = 0.0479175, 0.777768 
#    offsets[5] = 0.06389, 0.287436 
#    offsets[6] = -0.4632025, 0.490332 
#    offsets[7] = 0.9423775, 0.405792 
#    offsets[8] =-0.2395875, -0.135264 
#    offsets[9] = 0.76668, 0.4227
#    offsets[10] = -0.0479175, 0.388884
# for GW work
    offsets[2] = 0.06389, 0.287436 
    offsets[3] =-0.2395875, -0.135264 
    offsets[4] = -0.0479175, 0.388884
    offsets[5] = -0.76668, 0.473424 
    offsets[6] = -0.4632025, 0.490332 
    offsets[7] = 0.9423775, 0.405792 
    offsets[8] = -0.543065, -0.828492
    offsets[9] = 0.76668, 0.4227
    offsets[10] = 0.0479175, 0.777768 
# production offsets 9
#    offsets[9] = 0.76668, 0.4227
# production offsets 10
#    offsets[10] = -0.0479175, 0.388884
# production offsets 10
    offsets[11] = -0.5257, 0.7222
# production offsets 17,18
    offsets[17] = -1.1388, 0.0166
    offsets[18] =  0.0484, -0.6725
# fake offset 19
    offsets[19] =  0.9423775, -0.405792 
#
    return offsets



def test3(ra,dec, tile1, tile2, seqid="G184098-lmc-2", seqnum=0, seqtot=14) :
    exposureList = [5.,5.] 
    filterList = ["i","i"]
    jsonFilename="des-gw.json"
    offsets = tileOffsets()
    fd = open(jsonFilename,"w")
    fd.write("[\n")

    size = ra.size
    nexp = np.size(exposureList)
    for i in range(0,size) :
        for j in range(0,nexp) :
            seqnum +=1
#not clobbered
#count, seqid, seqnum, seqtot,note,comment
            if j == 0 : tiling = tile1[i]
            if j == 1 : tiling = tile2[i]
            filter = filterList[j]
            exp = exposureList[j]
            offsets[tiling]
            #print tiling, offsets[tiling]
            delRa = offsets[tiling][0]
            delDec = offsets[tiling][1]
            tra = ra[i]
            tdec = dec[i]
            tdec = tdec+delDec
            tra = tra + delRa/np.cos(tdec*2*np.pi/360.)
            comment = "DESGW: LIGO event {}: {} of {}".format(seqid, seqnum, seqtot)
            intra = np.int(np.round(np.int(ra[i]*10.)/10.)*10.)
            intdec = np.int(np.round(np.int(dec[i]*10.)/10.)*10.)
            signDec = "+"; 
            if tdec < 0: signDec = "-"
            object = "DESGW: LIGO event {}: {} of {}".format(seqid, seqnum, seqtot)

            fd.write("{")
            fd.write(" \"expType\" : \"object\",\n")
            fd.write("  \"object\" : \"{}\",\n".format(object))
            fd.write("  \"seqid\" : \"{}\",\n".format(seqid))
            fd.write("  \"seqnum\" : \"{:d}\",\n".format(int(seqnum)))
            fd.write("  \"seqtot\" : \"{:d}\",\n".format(int(seqtot)))
            fd.write("  \"expTime\" : {:d},\n".format(int(exp)))
            fd.write("  \"wait\" : \"False\",\n")
            fd.write("  \"count\" : \"1\",\n")
            fd.write("  \"note\" : \"Added to queue from desgw json file, not obstac\",\n")
            fd.write("  \"filter\" : \"{}\",\n".format(filter))
            fd.write("  \"program\" : \"des gw\",\n")
            fd.write("  \"RA\" : {:.6f},\n".format(tra))
            fd.write("  \"dec\" : {:.5f},\n".format(tdec))
            fd.write("  \"comment\" : \"{}\"\n".format(comment)) 
            # note lack of comma for end
            fd.write("}")
            if (i == size-1) and ( j == nexp-1) :
                pass
            else :
                fd.write(",")
            fd.write("\n")
    fd.write("]\n")
    fd.close()

def test2(ra,dec,tile1,tile2) :
    off = tileOffsets()
    new_ra, new_dec = np.array([]),np.array([])
    for i in range(0,ra.size) :
        t1 = np.int(tile1[i])
        t2 = np.int(tile2[i])
        o1 = off[t1]
        o2 = off[t2]
        cosd = np.cos(o1[1]*np.pi/180.)
        new_ra = np.append( new_ra, ra[i]+o1[0]/cosd)
        new_ra = np.append( new_ra, ra[i]+o2[0]/cosd)
        new_dec = np.append( new_dec, dec[i]+o1[1])
        new_dec = np.append( new_dec, dec[i]+o2[1])

    return new_ra, new_dec

# testing the jsons
#
#
def test() :
    import matplotlib.pyplot as plt
    import json
    mjdList = ("UTC-2015-9-20-6:29:00","UTC-2015-9-20-6:59:00","UTC-2015-9-20-7:29:00","UTC-2015-9-20-7:59:00","UTC-2015-9-20-8:29:00","UTC-2015-9-20-8:59:00")

    dir = "./night4/"
    f1="G184098-6-UTC-2015-9-20-6:29:00.json"
    f2="G184098-7-UTC-2015-9-20-6:59:00.json"
    f3="G184098-8-UTC-2015-9-20-7:29:00.json"
    f4="G184098-9-UTC-2015-9-20-7:59:00.json"
    f5="G184098-10-UTC-2015-9-20-8:29:00.json"
    f6="G184098-11-UTC-2015-9-20-8:59:00.json"

    j=open(dir+f1); 
    data6 = json.load(j);x=[(e["RA"],e["dec"]) for e in data6]; ra6,dec6 = zip(*x)
    j.close()
    j=open(dir+f2); 
    data7 = json.load(j);x=[(e["RA"],e["dec"]) for e in data7]; ra7,dec7 = zip(*x)
    j.close()
    j=open(dir+f3); 
    data8 = json.load(j);x=[(e["RA"],e["dec"]) for e in data8]; ra8,dec8 = zip(*x)
    j.close()
    j=open(dir+f4); 
    data9 = json.load(j);x=[(e["RA"],e["dec"]) for e in data9]; ra9,dec9 = zip(*x)
    j.close()
    j=open(dir+f5); 
    data10 = json.load(j);x=[(e["RA"],e["dec"]) for e in data10]; ra10,dec10 = zip(*x)
    j.close()
    j=open(dir+f6); 
    data11 = json.load(j);x=[(e["RA"],e["dec"]) for e in data11]; ra11,dec11 = zip(*x)
    j.close()


    plt.clf();plt.scatter(ra6,dec6,c="cyan");
    plt.scatter(ra7,dec7,c="b");
    plt.scatter(ra8,dec8,c="g");
    plt.scatter(ra9,dec9,c="orange");
    plt.scatter(ra10,dec10,c="r");
    plt.scatter(ra11,dec11,c="pink")

