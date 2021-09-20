def fixer(y,m, threshold): # at the moment, takes in a flux dict, int number of surrounding data points, int binarization threshold
    y_out = y.copy() # So we don’t overwrite y
    #threshold = 3 # 7 binarization threshold.
    #m=2 # used to select 2 m + 1 data points around our spike
    spikes = dict()
    for k in y.keys():
        spikes[k] = abs(numpy.array(modified_z_score(numpy.diff(y[k]))) > threshold)
        for i in numpy.arange(len(spikes[k])):
            if spikes[k][i] != 0: # if a spike is in position i
                w = numpy.arange(i-m,i+1+m) # select 2 m + 1 points around our spike
                w2 = w[spikes[k][w] == 0] # From such interval, choose the ones which are not spikes
                y_out[k][i] = numpy.mean(y[k][w2]) # average their values
    return y_out
                
def modified_z_score(b):
    median_int = numpy.median(b)
    mad_int = numpy.median([numpy.abs(b - median_int)])
    modified_z_scores = 0.6745 * (b - median_int) / mad_int #The multiplier 0.6745 is the 0.75th quartile of the standard normal distribution, to which the MAD converges
    return modified_z_scores

