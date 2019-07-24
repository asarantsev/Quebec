# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 11:37:59 2019

@author: UNR Math Stat
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:08:50 2019

@author: UNR Math Stat
"""

#We do both analysis of yearly average 
#and regression for individual plots, logarithms of biomass or shadow,
#NOT CENTERED by subtracting log of yearly avg, for all years

import pandas   # data analysis
import numpy    # n-dim arrays
import math
import matplotlib
from matplotlib import pyplot
from statsmodels.graphics.gofplots import qqplot
import scipy
from scipy import stats

def Ev(a, n):
    if (a == -1) or (a == 1):
        return math.sqrt(n)
    else:
        b = ( 1-a**(2*n) ) / ( 1-a**2 )  
        return math.sqrt(b)

def Od(a, n): 
    if (a == 1):
        return n
    else:
        return (1-a**(n))/(1-a)   
     
# This is our regression with given value of a 
def regression(a, logvalue, lognext, interval, NCHANGES):
    V = numpy.array([])
    A = numpy.array([])
    for j in range(NCHANGES): 
         d = interval[j]
         diff = lognext[j]-(a**d)*logvalue[j]
         V = numpy.append(V, diff/Ev(a, d))  
         A = numpy.append(A, Od(a, d)/Ev(a, d))
    r = numpy.dot(V, A)/numpy.dot(A, A)
    sigma = math.sqrt((1/(NCHANGES-1))*numpy.dot(V-r*A, V-r*A))
    return V, A, r, sigma

# QQ plot of residuals for this regression
# compare to standard normal distribution, 'r' - regression line    
def qq(a, logvalue, lognext, interval, NCHANGES):
    V, A, r, sigma = regression(a, logvalue, lognext, interval, NCHANGES)
    centralized = numpy.array([V[k] - r*A[k] for k in range(NCHANGES)])
    qqplot(centralized, line = 'r') 
    pyplot.show()
    


path = '/Users/Olga Rumyantseva/Desktop/Python biomass/'
df = pandas.read_csv(path + 'biodata.csv')
# df = pandas.read_csv('biodata.csv')
data = df.values  # work with values only
NOBSERV = numpy.size(data, 0) # number of rows in our dataframe (number of observations)
print('Number of all data points = ', NOBSERV) 



patch = data[:,0]   # plot ID (1st column)
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years converted in integer from float
biomass = data[:,2]
shadow = data[:,3]
value = biomass # this is what we consider now: biomass or shadow
#Change to shadow if you need


logvalue = numpy.array([]) # here will be logarithms of value for steps
lognext = numpy.array([])  # increments of logarithms 
interval = numpy.array([]) # increments of years

#next loop collects records from the same patch: 
#initial , change in logs, time interval 
#For example, if the same patch is observed in 1980 and 1987, 
#collect logvalue in 1980, logchange and time interval

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        logvalue = numpy.append(logvalue, math.log(value[observ]))
        lognext = numpy.append(lognext, math.log(value[observ+1]))
        interval = numpy.append(interval, year[observ+1] - year[observ])

       
#Number of pair-observations
NCHANGES = numpy.size(logvalue)  
print('Number of pair-observations = ', NCHANGES)
Years = numpy.unique(year)
NYears = numpy.size(Years)

#Average biomass or shadow in a given year
YearAvBiomass = numpy.array([])
for j in range(NYears):
    logbio = math.log(numpy.mean(value[year == Years[j]]))
    YearAvBiomass = numpy.append(YearAvBiomass, logbio)
    
#Increments of average biomass (or shadow)
DeltasYearAvBiomass = numpy.array([])  
for i in range(numpy.size(YearAvBiomass)-1):
    delt = YearAvBiomass[i+1]-YearAvBiomass[i]
    DeltasYearAvBiomass = numpy.append(DeltasYearAvBiomass, delt)  
    

print('mean change of avg year biomass = ', numpy.mean(DeltasYearAvBiomass))
#Mean of increments of average yearly biomass
print('stdev of this change = ', numpy.std(DeltasYearAvBiomass))
#Standard deviation of these increments




Sigmas = numpy.array([])    
Rs = numpy.array([])

x = numpy.arange(-2, 2, 0.01)
for a in x:
    V, A, r, sigma = regression(a, logvalue, lognext, interval, NCHANGES)
    Sigmas = numpy.append(Sigmas, sigma)
    Rs = numpy.append(Rs, r)
#Plot of parameter a vs standard error
pyplot.plot(x, Sigmas)
pyplot.show()
r = Rs[numpy.argmin(Sigmas)]#Corresponding r

print('Non-centered data for all years')
print('min std = ', numpy.min(Sigmas))
print('value of a = ', x[numpy.argmin(Sigmas)])
print('value of m = ', r)
print('qq plot for residuals for random walk case')
qq(1, logvalue, lognext, interval, NCHANGES)
#residuals = [(lognext[k] - logvalue[k])/math.sqrt(interval[k]) - m * math.sqrt(interval[k]) for k in range(NCHANGES)]

    