# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:48:35 2019

@author: Olga Rumyantseva
"""



import pandas   # data analysis
import numpy    # n-dim arrays
import math
import matplotlib
from matplotlib import pyplot
from statsmodels.graphics.gofplots import qqplot
import scipy    # lin regress
import matplotlib.pyplot as plt
from pandas.tools.plotting import table

import random

##################################################################
#################### Preparation Block: ##########################
##################################################################


path = '/Users/Olga Rumyantseva/Desktop/Python biomass/'
# df = pandas.read_csv(path + 'biodata.csv')
df = pandas.read_csv(path + 'biodata1.csv') # header added as the 

data = df.values  # work with values only
NOBSERV = numpy.size(data, 0) # number of rows in our dataframe (number of observations)
print(NOBSERV) 

patch = data[:,0]   # plot ID (1st column)
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years converted in integer from float
biomass = data[:,2]
shadow = data[:,3]
value = biomass # this is what we consider now: biomass or shadow

mu = 0.021   # mean   numpy.mean(DeltasYearAvBiomass)
sigma = 0.512 # standard deviation   numpy.std(DeltasYearAvBiomass)

from math import e

# Years = numpy.unique(year)
# m0 = numpy.mean(value[year == Years[0]]  # Years[0] == 1970

#  means and variances based on existing data:
ExpectationsB = numpy.array([]) 
VariancesB = numpy.array([])
ExpectationsS = numpy.array([]) 
VariancesS = numpy.array([])



#  means and variances based on existing data:
for t in range(38): 
       meanB = numpy.mean(biomass[year == 1970 + t])
       meanB = round(meanB, 2)
       ExpectationsB = numpy.append(ExpectationsB, meanB)
       
       vB = numpy.var(biomass[year == 1970 + t])
       vB = round(math.sqrt(vB), 2)
       VariancesB = numpy.append(VariancesB, vB)
       
       meanS = numpy.mean(shadow[year == 1970 + t])
       meanS = round(meanS, 2)
       ExpectationsS = numpy.append(ExpectationsS, meanS)
       
       vS = numpy.var(shadow[year == 1970 + t])
       vS = round(math.sqrt(vS), 2)
       VariancesS = numpy.append(VariancesS, vS)
       
       nobserv = numpy.size(biomass[year == 1970 + t])
       
       print(1970 + t,'&', nobserv, '&', meanB, '&', vB,'&',  meanS, '&', vS, '\\')



#  simulated means and variances:
ExpB = numpy.array([]) 
VarB = numpy.array([])
ExpS = numpy.array([]) 
VarS = numpy.array([])

for t in range(13): 
       m0 = numpy.mean(biomass[year == 2007])
       meanB = round(m0*(e**(t*mu + t*(sigma**2)*(1/2))), 1)
       meanB = round(meanB, 2)
       ExpB = numpy.append(ExpB, meanB)
       
       vB = (m0**2)*(e**(2*t*mu + t*(sigma**2)*(1/2)))*(e**(t*(sigma**2)) - 1)
       vB = round(math.sqrt(vB), 2)
       VarB = numpy.append(VarB, vB)
       
       m0 = numpy.mean(shadow[year == 1970])
       meanS = round(m0*(e**(t*mu + t*(sigma**2)*(1/2))), 2)
       meanS = round(meanS, 2)
       
       ExpS = numpy.append(ExpS, meanS)
       vS = (m0**2)*(e**(2*t*mu + t*(sigma**2)*(1/2)))*(e**(t*(sigma**2)) - 1)
       vS = round(math.sqrt(vS), 1)
       VarS = numpy.append(VarS, vS)
       
       print(2007+t,'&',  meanB, '&', vB,'&',  meanS, '&', vS, '\\')
       















