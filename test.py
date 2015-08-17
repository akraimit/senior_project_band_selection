import glob
import math
import os
import numpy
import time 
import itertools
import numpy.linalg
import scipy.linalg

os.chdir('tss_spectra/')
tssAll = numpy.zeros((60,len(glob.glob('*.txt'))))
tssFilenames = glob.glob('*.txt')
for i in range(0,len(tssFilenames)):
   tssAll[:,i] = numpy.loadtxt(tssFilenames[i])[::2]
tssMean = numpy.mean(tssAll, 1)
os.chdir('../')

print 'tssAll: '
print tssAll.shape
print tssAll
print (tssAll.T).shape

cov = numpy.cov(tssAll)
print cov
print cov.shape


print 'jm distance: \n\n\n'
orange = numpy.array([37, 27, 45, 30, 57, 48, 34, 50, 20, 53, 33, 25, 51])
lemon = numpy.array([12, 17, 20, 32, 16, 30, 30, 37, 25, 42, 13, 56, 13])
pear = numpy.array([41, 19, 15, 12, 15, 55, 33, 37, 40, 40, 43, 46, 54])
apple = numpy.array([38, 39, 12, 60, 34, 47, 13, 24, 30, 19, 57, 54, 55])

sample1 = numpy.array([1362, 1411, 1457, 1735, 1621, 1621, 1791, 1863, 1863, 1838])
sample2 = numpy.array([1354, 1458, 1458, 1458, 1550, 1145, 1428, 1573, 1573, 1657])

v1 = sample1
v2 = sample2

m1 = numpy.mean(v1)
c1 = numpy.asmatrix(numpy.cov(v1.T))
d1 = v1 - m1

m2 = numpy.mean(v2)
c2 = numpy.asmatrix(numpy.cov(v2.T))
d2 = v2 - m2

meanDif = numpy.asmatrix(m1 - m2)
cAv = (c1 + c2) / 2.0
# # test way in R
# a =  numpy.asarray(0.125 * (numpy.asmatrix(oDif)).T * (1./ numpy.asmatrix(meanDif)))

# b =  + 0.5 * numpy.log (numpy.linalg.det ( covAv ) / numpy.sqrt(numpy.linalg.det( oCov ) * numpy.linalg.det ( lCov )))
# bh = a + b 
# print 2 * ( 1 - numpy.exp ( -bh ) )

# part1 = 0.125 * (meanDif) * 1./cAv * meanDif.T 
# part2 = (0.5) * numpy.log10( numpy.abs(cAv) / (numpy.dot(numpy.sqrt(numpy.abs(c1)),numpy.sqrt(numpy.abs(c2)))))

# a =  part1 + part2 


print ' \n'

# # MATLAB way
# d = (m1-m2)/(numpy.linalg.cholesky(cAv).T)
# a= 0.125*d*d.T + 0.5*numpy.log(numpy.abs(numpy.linalg.det(cAv/ scipy.linalg.sqrtm(c1*c2))));
# #a = 0.125 * d *d.T + 0.5 * numpy.log10(numpy.linalg.det(cAv/numpy.linalg.cholesky(c1*c2).T))
# #a = 0.125 * d * d.T + 0.25 *  numpy.log(numpy.linalg.det((c1+c2)/2)**2 / (numpy.linalg.det(c1) * numpy.linalg.det(c2)))

# jm = 2.0 * (1-numpy.exp(-a))
# print 'jm: '
# print jm






bh = 0.125 * (m1-m2) * numpy.linalg.inv(cAv) * (m1-m2).T + 0.5*numpy.log(numpy.abs(numpy.linalg.det(cAv/ scipy.linalg.sqrtm(c1*c2))));

jm = 2.0 * (1-numpy.exp(-bh))

print 'bh: '
print bh
print 'jm: '
print jm


