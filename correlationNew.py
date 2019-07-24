# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:24:37 2019

@author: UNR Math Stat
"""

import pandas   # data analysis
import numpy    # n-dim arrays
import math
import matplotlib
from matplotlib import pyplot
from statsmodels.graphics.gofplots import qqplot
import scipy
from scipy import stats


path = '/Users/Olga Rumyantseva/Desktop/Python biomass/'
df = pandas.read_csv('data.csv')
# df = pandas.read_csv('biodata.csv')
data = df.values  # work with values only
NOBSERV = numpy.size(data, 0) # number of rows in our dataframe (number of observations)


patch = data[:,0]   # plot ID (1st column)
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years converted in integer from float
biomass = data[:,2]
shadow = data[:,3]


value = biomass # this is what we consider now: biomass or shadow
logvalue = numpy.array([]) # here will be logarithms of value for steps
lognext = numpy.array([])  # increments of logarithms 
interval = numpy.array([]) # increments of years

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        Y0 = year[observ]
        Y1 = year[observ+1]
        if not (Y0 == 1982 or Y0 == 2004 or Y1 == 1982 or Y1 == 2004): 
            logvalue = numpy.append(logvalue, math.log(value[observ]))
            lognext = numpy.append(lognext, math.log(value[observ+1]))
            interval = numpy.append(interval, Y1 - Y0)

Years = numpy.unique(year)
NYears = numpy.size(Years)

YearAvBiomass = numpy.array([])
for j in range(NYears):
    if not (Years[j] == 1982 or Years[j] == 2004):
        logbio = math.log(numpy.mean(value[year == Years[j]]))
        YearAvBiomass = numpy.append(YearAvBiomass, logbio)
   
    

value = shadow # this is what we consider now: biomass or shadow
logvalue1 = numpy.array([]) # here will be logarithms of value for steps
lognext1 = numpy.array([])  # increments of logarithms 
interval1 = numpy.array([]) # increments of years

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        Y0 = year[observ]
        Y1 = year[observ+1]
        if not (Y0 == 1982 or Y0 == 2004 or Y1 == 1982 or Y1 == 2004): 
            logvalue1 = numpy.append(logvalue1, math.log(value[observ]))
            lognext1 = numpy.append(lognext1, math.log(value[observ+1]))
            interval1 = numpy.append(interval1, Y1 - Y0)

YearAvBiomass1 = numpy.array([])
for j in range(NYears):
    if not (Years[j] == 1982 or Years[j] == 2004):
        logbio = math.log(numpy.mean(value[year == Years[j]]))
        YearAvBiomass1 = numpy.append(YearAvBiomass1, logbio)
          


print(numpy.corrcoef(logvalue, logvalue1)[0][1])
print(numpy.corrcoef(YearAvBiomass, YearAvBiomass1)[0][1])
