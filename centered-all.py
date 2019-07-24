# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:08:50 2019

@author: UNR Math Stat
"""


#We do both analysis of yearly average 
#and regression for individual plots, logarithm of biomass or shadow
#centered by subtracting logarithm of yearly average, for all years

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
     
#This is our regression with given value of parameter a
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

#QQ plot for residuals for the regression above
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
print(NOBSERV) 

patch = data[:,0]   # plot ID (1st column)
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years converted in integer from float
biomass = data[:,2]
shadow = data[:,3]
value = biomass # this is what we consider now: biomass or shadow
#change to shadow if you wish



logvalue = numpy.array([]) # here will be logarithms of value for steps
lognext = numpy.array([])  # increments of logarithms 
interval = numpy.array([]) # increments of years

#next loop collects records from the same patch: 
#initial, change in logs, time interval 
#For example, if observed in 1980 and 1987, 
#then collect logvalue in 1980, logchange and time interval

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        logvalue = numpy.append(logvalue, math.log(value[observ]))
        lognext = numpy.append(lognext, math.log(value[observ+1]))
        interval = numpy.append(interval, year[observ+1] - year[observ])

       
NCHANGES = numpy.size(logvalue)  
print('Number of observation pairs for same patch = ', NCHANGES)
Years = numpy.unique(year)
NYears = numpy.size(Years)

# Average biomass or shadow in a given year 
YearAvBiomass = numpy.array([])
for j in range(NYears):
    logbio = math.log(numpy.mean(value[year == Years[j]]))
    YearAvBiomass = numpy.append(YearAvBiomass, logbio)    


#Increments of average biomass (or shadow)
DeltasYearAvBiomass = numpy.array([])  
for i in range(numpy.size(YearAvBiomass)-1):
    delt = YearAvBiomass[i+1]-YearAvBiomass[i]
    DeltasYearAvBiomass = numpy.append(DeltasYearAvBiomass, delt)  

      
MeanDeltasAvg = numpy.mean(DeltasYearAvBiomass)
StdDeltasAvg = numpy.std(DeltasYearAvBiomass)
print('Mean of yearly average change = ', MeanDeltasAvg)
print('Stdev of yearly average change = ', StdDeltasAvg)
print('QQ plot for changes in yearly average')
qqplot(DeltasYearAvBiomass, line = 'r')
pyplot.show()
print('shapiro-wilk test = ', stats.shapiro(DeltasYearAvBiomass))
# This shows that these yearly changes obey normal distributions
# The test statistic 0.9816, 
# p-value 0.7867 - big, so we can't reject the
# null-hypothesis => normally distributed

####################################################################
####################################################################
Patches  = numpy.unique(patch)
NPatches = numpy.size(Patches)

# Years = numpy.unique(year)
# NYears = numpy.size(Years)

n = NYears  
m = NPatches

LogValuePatchYear = [[0] * m for i in range(n)] 
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
      

# Z[i][p] = log(biomass at plot p in year i) - Log(YearAvBiomass(i)): 
Z = [[0] * m for i in range(n)] 
for j in range(m):
    for i in range(n):
        Z[i][j] = LogValuePatchYear[i][j] - YearAvBiomass[i]

#####################################################################
############## Step 3  Analysis for Z instead of y_p(t):  ################
#####################################################################        

logvalueCentered = numpy.array([]) # here will be logarithms of value for steps
lognextCentered = numpy.array([])  # increments of logarithms 
intervalCentered = numpy.array([]) # increments of years

#next loop collects records from the same patch: 
# initial , 
# change in logs, 
# time interval :


for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        j = numpy.where(Patches == patch[observ])[0][0] 
        i = numpy.where(Years == year[observ])[0][0] 
        i1 = numpy.where(Years == year[observ+1])[0][0]
        logvalueCentered = numpy.append(logvalueCentered, Z[i][j])
        lognextCentered = numpy.append(lognextCentered, Z[i1][j])
        intervalCentered = numpy.append(intervalCentered, year[observ+1] - year[observ])

       
#  NCHANGES = numpy.size(logvalue) 


SigmasCentered = numpy.array([])    
RsCentered = numpy.array([])

x = numpy.arange(-2, 2, 0.01)
for a in x:
    V, A, r, sigma = regression(a, logvalueCentered, lognextCentered, intervalCentered, NCHANGES)
    SigmasCentered = numpy.append(SigmasCentered, sigma)
    RsCentered = numpy.append(RsCentered, r)

pyplot.plot(x, SigmasCentered)
pyplot.show()

print('Centered forest patch data for all years')
print('min std = ', numpy.min(SigmasCentered))
print('value of a = ', x[numpy.argmin(SigmasCentered)])
print('value of m = ', RsCentered[numpy.argmin(SigmasCentered)])
print('qq plot for residuals of this regression')
qq(1, logvalueCentered, lognextCentered, intervalCentered, NCHANGES)

    