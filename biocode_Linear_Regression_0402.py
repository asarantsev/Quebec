# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import pandas   # data analysis
import numpy    # n-dim arrays
import math
import matplotlib
from matplotlib import pyplot
from statsmodels.graphics.gofplots import qqplot
import scipy    # lin regress

path = '/Users/Olga Rumyantseva/Desktop/Python biomass/'
# df = pandas.read_csv(path + 'biodata.csv')


df = pandas.read_csv('biodata.csv') # header added as the 
# very first row doesn't read anyway and gets missed reading with biodata.csv

# checking that there is no duplicates in data frame:
#numpy.size(df, 0) - numpy.size(df.drop_duplicates(keep=False), 0) 



data = df.values  # work with values only
print(data[0][0]) # 1st row 1st column
NOBSERV = numpy.size(data, 0) # number of rows in our dataframe (number of observations)
print(NOBSERV) 



patch = data[:,0]   # plot ID (1st column)
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years converted in integer from float
biomass = data[:,2]
shadow = data[:,3]
value = shadow # this is what we consider now: biomass or shadow



logvalue = numpy.array([]) # here will be logarithms of value for steps
lognext = numpy.array([])  # increments of logarithms 
interval = numpy.array([]) # increments of years

#next loop collects records from the same patch: initial , change in logs, time interval (if observed
# in 1980 and 1987, collect logvalue in 1980, logchange and time interval)

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        logvalue = numpy.append(logvalue, math.log(value[observ]))
        lognext = numpy.append(lognext, math.log(value[observ+1]))
        interval = numpy.append(interval, year[observ+1] - year[observ])

       

NCHANGES = numpy.size(logvalue)  
#quantity of all pairs of years for the same plot

counted = [numpy.count_nonzero([gap == k for gap in interval]) for k in range(max(year)-min(year))]
print(counted)
#by a certain number of years
#quantity of patches with pairs of observations with interval k years, for k = 0, ..., 39

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
     


Sigmas = numpy.array([])    
Ms = numpy.array([])

def regression(a):
    V = numpy.array([])
    for j in range(NCHANGES): 
         d = interval[j]
         diff = lognext[j]-(a**d)*logvalue[j]
         V = numpy.append(V, diff/Ev(a, d))  

    A = numpy.array([])
    for j in range(NCHANGES):
         d = interval[j]
         A = numpy.append(A, Od(a, d)/Ev(a, d))

    m = numpy.dot(V, A)/numpy.dot(A, A)
    sigma = (1/(NCHANGES-1))*numpy.dot(V-m*A, V-m*A)
    return V, A, m, sigma

def qq(a):
    V, A, m, sigma = regression(a)
    centralized = numpy.array([V[k] - m*A[k] for k in range(NCHANGES)])
    qqplot(centralized, line = 'r')
    pyplot.show()

x = numpy.arange(-2, 2, 0.01)    
for a in x:
    V, A, m, sigma = regression(a)
    Ms = numpy.append(Ms, m)
    Sigmas = numpy.append(Sigmas, sigma)


# print(Sigmas)
print('all data, no centering')
print('min std = ', numpy.min(Sigmas))
n = numpy.argmin(Sigmas)
print('value of a = ', x[n])
print('value of m = ', Ms[n])

pyplot.plot(x, Sigmas)
pyplot.show()

print('a = ', x[n])
qq(x[n])

#Sigmas = numpy.array([])    
#Ms = numpy.array([])
#
#y = [0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00]
#for a in y:
#    V, A, m, sigma = regression(a)
#    Ms = numpy.append(Ms, m)
#    Sigmas = numpy.append(Sigmas, sigma)
#    
## print(Sigmas)
#print('all data, no centering')
#print('min std = ', numpy.min(Sigmas))
#m = numpy.argmin(Sigmas)
#print('value of a = ', y[m])
#print('value of m = ', Ms[m])
#
#pyplot.plot(y, Sigmas)
#pyplot.show()
#
#print('a = ', y[m])
#qq(y[m])
#
#print('a = ', 1)
#qq(1)

Years = numpy.unique(year)
NYears = numpy.size(Years)


YearAvBiomass = numpy.array([])

for j in range(NYears):
    logbio = math.log(numpy.mean(value[year == Years[j]]))
    YearAvBiomass = numpy.append(YearAvBiomass, logbio)
    

DeltasYearAvBiomass = numpy.array([])  

for i in range(numpy.size(YearAvBiomass)-1):
    delt = YearAvBiomass[i+1]-YearAvBiomass[i]
    DeltasYearAvBiomass = numpy.append(DeltasYearAvBiomass, delt)  
    

numpy.mean(DeltasYearAvBiomass)
numpy.std(DeltasYearAvBiomass)

# matplotlib.pyplot.plot(DeltasYearAvBiomass)

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

logvalCentered = numpy.array([]) # here will be logarithms of value for steps
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
        logvalCentered = numpy.append(logvalCentered, Z[i][j])
        lognextCentered = numpy.append(lognextCentered, Z[i1][j])
        intervalCentered = numpy.append(intervalCentered, year[observ+1] - year[observ])

       
#  NCHANGES = numpy.size(logvalue) 


SigmasCentered = numpy.array([])    
MsCentered = numpy.array([])

x = numpy.arange(-2, 2, 0.01)
for a in x:
    VCentered = numpy.array([])
    for j in range(NCHANGES): 
         d = interval[j]
         diff = lognextCentered[j]-(a**d)*logvalCentered[j]
         VCentered = numpy.append(VCentered, diff/Ev(a, d))  

    A = numpy.array([])
    for j in range(NCHANGES):
         d = interval[j]
         A = numpy.append(A, Od(a, d)/Ev(a, d))

    m = numpy.dot(VCentered, A)/numpy.dot(A, A)
    sigma = (1/(NCHANGES-1))*numpy.dot(VCentered-m*A, VCentered-m*A)
         
    MsCentered = numpy.append(MsCentered, m)
    SigmasCentered = numpy.append(SigmasCentered, sigma)


# print(SigmasCentered)

#print(numpy.min(SigmasCentered))

#print(x[numpy.argmin(SigmasCentered)]) # a = 1 precisely!! so
# this is random walk
print('All data, centered')
print('min std = ', numpy.min(SigmasCentered))
print('value of a = ', x[numpy.argmin(SigmasCentered)])
print('value of m = ', Ms[numpy.argmin(SigmasCentered)])

pyplot.plot(x, SigmasCentered)
pyplot.show()

def qqCentered(a):
    centralized = numpy.array([VCentered[k] - m*A[k] for k in range(NCHANGES)])
    qqplot(centralized, line = 'r')

##############################################################
##############  Step 4  Linear Regression  ###################
##############################################################        

# there is no records for a certain patch with one year difference d:
# there is 1 record for d = 2,
# there are 4 records for d = 3:

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        d = year[observ+1] - year[observ]
        if d == 1:
           print(patch[observ], observ)



#for observ in range(NOBSERV-1):
#    if patch[observ] == patch[observ+1]:
#        d = year[observ+1] - year[observ]
#        if d == 1:
#           j = numpy.where(Patches == patch[observ])[0][0] 
#           i = numpy.where(Years == year[observ])[0][0] 
#           deltaZ = numpy.append(deltaZ, Z[i+1][j] - Z[i][j])
        


##############################################################
##############  Step 5  ######################################
##############################################################        

# Throw away all data for all years for which we have less than 100
# observations: 

# count[i] = number of observations for the year with id i
count = [0] * NYears 
   
for observ in range(NOBSERV-1):
    j = numpy.where(Years == year[observ])[0][0]
    count[j] = count[j] + 1
    if year[observ] == year[observ+1]:
        count[j] = count[j] + 1       



## drop data for years having less than 100 observations:
#df1 = df.values  # keep df back up      
#for i in range(NYears):
#    if count[i] < 100:
#       df1 = df1[df1[i, 0] != Years[i]]
#
#data1 = df1
#
## compare sizes of original data and reduced data:
#numpy.size(data, 0) - numpy.size(data1, 0)   


############################################################
############### Repeat steps 1,3,4: ########################
############################################################

data1 = data

NOBSERV = numpy.size(data, 0) # number of rows in our dataframe (number of observations)
print(NOBSERV) 

patch = data[:,0]   # plot ID (1st column)
year = numpy.array([int(data[k,1]) for k in range(NOBSERV)]) #years converted in integer from float
biomass = data[:,2]
shadow = data[:,3]

value = biomass # this is what we consider now: biomass or shadow
threshold = 100
censoredpatch = []
censoredvalue = []
censoredyear = []
for k in range(NOBSERV):
    if not (year[k] == 1982 or year[k] == 2004):
        censoredpatch.append(patch[k])
        censoredvalue.append(value[k])
        censoredyear.append(year[k])
        
patch = censoredpatch
value = censoredvalue
year = censoredyear

NOBSERV = numpy.size(patch)
print('new NOBSERV = ', NOBSERV)

logvalue = numpy.array([]) # here will be logarithms of value for steps
lognext = numpy.array([])  # increments of logarithms 
interval = numpy.array([]) # increments of years

#next loop collects records from the same patch: initial , change in logs, time interval (if observed
# in 1980 and 1987, collect logvalue in 1980, logchange and time interval)

for observ in range(NOBSERV-1):
    if patch[observ] == patch[observ+1]:
        logvalue = numpy.append(logvalue, math.log(value[observ]))
        lognext = numpy.append(lognext, math.log(value[observ+1]))
        interval = numpy.append(interval, year[observ+1] - year[observ])

       
NCHANGES = numpy.size(logvalue)  




Sigmas = numpy.array([])    
Ms = numpy.array([])

x = numpy.arange(-2, 2, 0.01)
for a in x:
    V = numpy.array([])
    for j in range(NCHANGES): 
         d = interval[j]
         diff = lognext[j]-(a**d)*logvalue[j]
         V = numpy.append(V, diff/Ev(a, d))  

    A = numpy.array([])
    for j in range(NCHANGES):
         d = interval[j]
         A = numpy.append(A, Od(a, d)/Ev(a, d))

    m = numpy.dot(V, A)/numpy.dot(A, A)
    sigma = (1/(NCHANGES-1))*numpy.dot(V-m*A, V-m*A)
         
    Ms = numpy.append(Ms, m)
    
    Sigmas = numpy.append(Sigmas, sigma)


# print(Sigmas)
print('Censored data, not centered')
print('sigma = ', numpy.min(Sigmas))
print('a = ', x[numpy.argmin(Sigmas)])
print('m = ', Ms[numpy.argmin(Sigmas)])

pyplot.plot(x, Sigmas)
pyplot.show()


Years = numpy.unique(year)
NYears = numpy.size(Years)


YearAvBiomass = numpy.array([])

for currentyear in Years:
    temp = 0
    num = 0
    for k in range(NOBSERV):
        if (year[k] == currentyear):
            temp = temp + value[k]
            num = num + 1
    YearAvBiomass = numpy.append(YearAvBiomass, temp/num)
    

DeltasYearAvBiomass = numpy.array([])  

for i in range(numpy.size(YearAvBiomass)-1):
    delt = YearAvBiomass[i+1]-YearAvBiomass[i]
    DeltasYearAvBiomass = numpy.append(DeltasYearAvBiomass, delt)  
    

numpy.mean(DeltasYearAvBiomass)
numpy.std(DeltasYearAvBiomass)

# matplotlib.pyplot.plot(DeltasYearAvBiomass)


##############################################################
############# Step 3 repetition: #############################
##############################################################

logvalCentered = numpy.array([]) # here will be logarithms of value for steps
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
        logvalCentered = numpy.append(logvalCentered, Z[i][j])
        lognextCentered = numpy.append(lognextCentered, Z[i1][j])
        intervalCentered = numpy.append(intervalCentered, year[observ+1] - year[observ])

       
#  NCHANGES = numpy.size(logvalue) 


SigmasCentered = numpy.array([])    
MsCentered = numpy.array([])

x = numpy.arange(-2, 2, 0.01)
for a in x:
    VCentered = numpy.array([])
    for j in range(NCHANGES): 
         d = interval[j]
         diff = lognextCentered[j]-(a**d)*logvalCentered[j]
         VCentered = numpy.append(VCentered, diff/Ev(a, d))  

    A = numpy.array([])
    for j in range(NCHANGES):
         d = interval[j]
         A = numpy.append(A, Od(a, d)/Ev(a, d))

    m = numpy.dot(VCentered, A)/numpy.dot(A, A)
    sigma = (1/(NCHANGES-1))*numpy.dot(VCentered-m*A, VCentered-m*A)
         
    MsCentered = numpy.append(MsCentered, m)
    SigmasCentered = numpy.append(SigmasCentered, sigma)


# print(SigmasCentered)
print('Censored and centered data')
print('stderr = ', numpy.min(SigmasCentered))
print('a = ', x[numpy.argmin(SigmasCentered)]) 
print('m = ', MsCentered[numpy.argmin(SigmasCentered)])
# a = 1 precisely!! so
# this is random walk
pyplot.plot(x, SigmasCentered)
pyplot.show()
