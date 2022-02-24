#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 20:00:24 2021

@author: cl
"""


import numpy as np
import pandas as pd

#First I read in the data
filepath = '/Users/clepart/Downloads/desidiff/2021/'

fName_lux = filepath + 'difflux.csv'
fName_WL = filepath + 'difwave.csv'
data_lux = pd.read_csv(fName_lux, delimiter=',', header=None)
data_WL = pd.read_csv(fName_WL, delimiter=',', header=None)

#transpose pandas dataframe from rows to columns
dataY = data_lux.T
dataX = data_WL.T

#pull out sub-bands after converting from a dataframe to a numpy array
dataX = dataX.to_numpy()
dataY = dataY.to_numpy()

BdataX = dataX[:,0]
BdataY = dataY[:,0]

RdataX = dataX[:,1]
RdataY = dataY[:,1]

ZdataX = dataX[:,2]
ZdataY = dataY[:,2]

#eliminate nan values, note this changes the sizes of the arrays (except B)
BdataX = BdataX[~np.isnan(BdataX)]
RdataX = RdataX[~np.isnan(RdataX)]
ZdataX = ZdataX[~np.isnan(ZdataX)]

BdataY = BdataY[~np.isnan(BdataY)]
RdataY = RdataY[~np.isnan(RdataY)]
ZdataY = ZdataY[~np.isnan(ZdataY)]

#Now ready to try the numpy.diff on each band, starting with B (old code commented out)
 #yc = y.copy()
    #newmask = dict(mask)
    #for b in yc.keys():
    #    if b == 'z':
            #yb = BdataY
#Here's the main part that isn't just reformatting the data:
Bmask = np.ones(len(BdataY)) # setup the mask array with all ones and same size as BdataY
ybs = BdataY  # get a copy of the data
dely = abs(np.diff(ybs)) # find the diff
# go through dely values looking for those >9 and set mask to zero for those values only
for i in range(len(dely)):
    if dely[i] > 10:
        Bmask[i+1] = 0

newBdataY = np.multiply(BdataY, Bmask) #multiply each element of array by mask

Rmask = np.ones(len(RdataY)) # setup the mask array with all ones and same size as BdataY
ybs = RdataY  # get a copy of the data
dely = abs(np.diff(ybs)) # find the diff
# go through dely values looking for those >9 and set mask to zero for those values only
for i in range(len(dely)):
    if dely[i] > 8:
        Rmask[i+1] = 0

newRdataY = np.multiply(RdataY, Rmask) #multiply each element of array by mask

Zmask = np.ones(len(ZdataY)) # setup the mask array with all ones and same size as BdataY
ybs = ZdataY  # get a copy of the data
dely = abs(np.diff(ybs)) # find the diff
# go through dely values looking for those >9 and set mask to zero for those values only
for i in range(len(dely)):
    if dely[i] > 3:
        Zmask[i+1] = 0

newZdataY = np.multiply(ZdataY, Zmask) #multiply each element of array by mask

import matplotlib.pyplot as plt

plt.plot(BdataX, newBdataY)
plt.plot(RdataX, newRdataY)
plt.plot(ZdataX, newZdataY)