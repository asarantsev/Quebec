# -*- coding: utf-8 -*-
"""
Created on Sun May  5 15:41:21 2019

@author: UNR Math Stat
"""
import numpy
import math
import scipy
import pandas   
import matplotlib.pyplot as plt


from numpy import random
from numpy import mean
from matplotlib import pyplot
from statsmodels.graphics.gofplots import qqplot
from scipy import stats


def inverseGamma(alpha, beta):
    return (1/numpy.random.gamma(alpha, 1/beta))



def Histograms(meansB, varsB, t):
    msB = numpy.array([])
    vsB = numpy.array([])

    for j in range(NSIMS):
        msB = numpy.append(msB, meansB[t][j])      
        vsB = numpy.append(vsB, varsB[t][j])
    time = Years[t]
    
    plt.hist(msB, bins = 50) 
    plt.title("histogram for " + str(NSIMS) + " simulated means \n of ln(biomass) in " + str(time))
    plt.show()
    
    plt.hist(vsB, bins = 50) 
    plt.title("histogram for " + str(NSIMS) + " simulated variances \n of ln(biomass) in " + str(time))
    plt.show()


path = '/Users/Olga Rumyantseva/Desktop/Python biomass/'
df = pandas.read_csv(path + 'biodata.csv')

# df = pandas.read_csv('data.csv')


data = df.values  
NOBSERV = numpy.size(data, 0) # number of rows in our dataframe (number of observations)
print(NOBSERV) 

patch = data[:,0]   # plot ID 
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years 
biomass = data[:,2]
shadow = data[:,3]
value = biomass # this is what we consider now: biomass or shadow


Patches  = numpy.unique(patch)
NPatches = numpy.size(Patches)

Years = numpy.unique(year)
# 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980,
#       1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991,
#       1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
#       2003, 2004, 2005, 2006, 2007
NYears = numpy.size(Years)


LogValuePatchYear = [[0] * NPatches for i in range(NYears)] 
# [[p1,..., pm]  [p1,..., pm] .... [p1,..., pm]]
#     year_1       year_2             year_n


#   LogValuePatchYear[i][p] = log(biomass at plot p in year i): 
for k in range(NOBSERV):
    indY = numpy.where(Years == year[k])
    i = indY[0][0]     # return i such that year[k] == Years[i]
    indP = numpy.where(Patches == patch[k])
    j = indP[0][0]  # return j such that patch[k] == Patches[j]
    LogValuePatchYear[i][j] = math.log(value[k])
#   i - year, j - patch

#  LogValuePatchYear[year][patch]
 
############################################################################
############################################################################
############################################################################  

NSIMS = 1000 # Number of simulations of Bayesian estimates

# means simulated by Bayes:
meansB =  [[0] * NSIMS for i in range(NYears)] 
# [[m_1,..., m_1000]  .... [m_1,..., m_1000]]
#     year_1                  year_38
#  meansB[year][simulation]


# vars simulated by Bayes:
varsB = [[0] * NSIMS for i in range(NYears)]
#  varsB[year][simulation]


for t in range(NYears):    
    LogValuesYear = numpy.array([])  # biomasses vector x_{1}(t), ... , x_{#p}(t)
    for p in range(NPatches):    
        LogValuesYear = numpy.append(LogValuesYear, LogValuePatchYear[t][p])
    LogValuesYear = LogValuesYear[numpy.nonzero(LogValuesYear)] #  biomass data doesn't  contain 0 records, so we can remove zeroes:
    print('year = ', t)
    empMean = numpy.mean(LogValuesYear) # mean of biomasses in year t
    empVar = numpy.var(LogValuesYear) # variance of biomasses in year t
    print('mean = ', empMean)
    print('var = ', empVar)
    n = len(LogValuesYear) # number of patches observed in year t
    print('number = ', n)
    for j in range(NSIMS):
        varsB[t][j] = inverseGamma((n-1)/2, (n*empVar)/2) 
        meansB[t][j] = random.normal(empMean, varsB[t][j]/n)


Histograms(meansB, varsB, 5)
        
         
# g^hat - estimate:
summ = 0
for i in range(NSIMS):
    summ = summ + meansB[NYears-1][i]-meansB[0][i]  

empG = summ/(NSIMS*NYears)
print('global mean = ', empG)

# sigma^hat - estimate:
summ = 0
for i in range(NSIMS):
    for t in range(NYears):
        summ = summ + (meansB[t][i]-meansB[t-1][i] - empG)**2

empGlobalV = summ/(NSIMS*NYears)    
print('global stdev = ', empGlobalV) 


# G, Sigma - estimates, updated by Bayes:
Sigma = inverseGamma((NSIMS*NYears-1)/2, (NSIMS*NYears*empGlobalV)/2)
G = random.normal(empG, Sigma/(NSIMS*NYears))

