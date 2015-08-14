import glob
import math
import os
import numpy
import time 
import itertools
import numpy.linalg

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

v1 = orange
v2 = lemon

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

print ' \n\n\n'

part1 = 0.125 * (meanDif) * 1./cAv * meanDif.T 
part2 = (0.5) * numpy.log10( numpy.abs(cAv) / (numpy.dot(numpy.sqrt(numpy.abs(c1)),numpy.sqrt(numpy.abs(c2)))))
   
a =  part1 + part2 

jm = 2.0 * (1-numpy.exp(-a))
print jm
print ' '
print 'c1'
print c1
   


