import numpy as np

#===================================================================
#
# neutron star tidal remnant models
#
# Example:
#    models = getModels()
#    time  = 1.0 ;# days past Time Zero
#    models_at_t = modelsAtTimeT (models, time)
#    abs_mag = models_at_t[0]
#
#
# tell me what the max absMag is for the disk B+K model
def maxAbsMag(models) :
    time1,mags1=models["model1"]
    a=np.argmin(mags1)
    return mags1[a], time1[a]

def getModels() :
    models = readModels()
    models = interpolateModels(models)
    return models

def modelsAtTimeT (models, time_in_days) :
    mags1_interp = models["mags1_interp"] 
    mags2_interp = models["mags2_interp"] 
    mags3_interp = models["mags3_interp"] 
    mags4_interp = models["mags4_interp"] 

    mag1 = mags1_interp(time_in_days)
    mag2 = mags2_interp(time_in_days)
    mag3 = mags3_interp(time_in_days)
    mag4 = mags4_interp(time_in_days)
    return mag1, mag2, mag3, mag4

def interpolateModels (models) :
    from scipy.interpolate import interp1d
    time1, mags1 = models["model1"]
    time2, mags2 = models["model2"]
    time3, mags3 = models["model3"]
    time4, mags4 = models["model4"]

    mags1_interp = interp1d(time1, mags1, kind="cubic")
    mags2_interp = interp1d(time2, mags2, kind="cubic")
    mags3_interp = interp1d(time3, mags3, kind="cubic")
    mags4_interp = interp1d(time4, mags4, kind="cubic")

    models["mags1_interp"] = mags1_interp
    models["mags2_interp"] = mags2_interp
    models["mags3_interp"] = mags3_interp
    models["mags4_interp"] = mags4_interp

    maxAbs, time_maxAbs = maxAbsMag(models)
    models["scale"] = maxAbs, time_maxAbs
    return models

def readModels (dir = "/data/des30.a/data/annis/des-gw/ligo/fourModels/") :
    import cPickle as pickle

    file=open(dir+"model-m1.pickle","rb"); model1 = pickle.load(file); file.close()
    file=open(dir+"model-m2.pickle","rb"); model2 = pickle.load(file); file.close()
    file=open(dir+"model-m3.pickle","rb"); model3 = pickle.load(file); file.close()
    file=open(dir+"model-m4.pickle","rb"); model4 = pickle.load(file); file.close()

    modelDistance = 100.
    distanceModulus = -5*np.log10(modelDistance * 1e6/10.)

    mags_d1=model1["disk-barnes"]
    mags_d2=model2["disk-barnes"]
    mags_d3=model3["disk-barnes"]
    mags_d4=model4["disk-barnes"]
    time_d1=model1["disk-time"];
    time_d2=model2["disk-time"];
    time_d3=model3["disk-time"];
    time_d4=model4["disk-time"];
    disk1=model1["disk"];
    disk2=model2["disk"];
    disk3=model3["disk"];
    disk4=model4["disk"];

    mags_w3=model3["wind-barnes"]
    time_w3=model3["wind-time"];
    wind3=model3["wind"]
    mags_3 = combine_mags(mags_d3, mags_w3, time_d3, time_w3)
    time3 = time_d3

    mags = mags_d1, mags_d2, mags_3, mags_d4
    times = time_d1, time_d2, time3, time_d4
    disk = disk1, disk2, disk3, disk4
    wind = wind3

    filter = "i"
    mags1, mags2, mags3, mags4 = mags[0], mags[1], mags[2], mags[3]
    time1, time2, time3, time4 = times[0], times[1], times[2], times[3]
    model1 = time1, mags1["{}".format(filter)]+distanceModulus
    model2 = time2, mags2["{}".format(filter)]+distanceModulus
    model3 = time3, mags3["{}".format(filter)]+distanceModulus
    model4 = time4, mags4["{}".format(filter)]+distanceModulus

    # time is in days
    ans = dict()
    ans["model1"] = model1
    ans["model2"] = model2
    ans["model3"] = model3
    ans["model4"] = model4
    return ans


def combine_mags(disk_mags, wind_mags, time_d, time_w) :
    from scipy.interpolate import interp1d
    sum_mags = dict()
    for filter in ["g", "r", "i", "z", "y"] :
        mag_d = disk_mags["{}".format(filter)]
        mag_w = wind_mags["{}".format(filter)]
        mag_w_interp = interp1d(time_w, mag_w, kind="cubic")
        sum = -2.5*np.log10( 10**(-0.4*mag_d) + 10**(-0.4*mag_w_interp(time_d)))
        sum_mags[filter] = sum
    return sum_mags

