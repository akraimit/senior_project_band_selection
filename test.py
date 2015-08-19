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


###############  TESTING BY HAND 
x = numpy.array([[2,2,3,4,5,5],[2,3,4,3,5,4]])
m1 = numpy.mean(x,1)
c1 = numpy.asmatrix(numpy.cov(x))
print 'mean 1: '
print m1
print 'cov 1: '
print c1
print 'cov 1 T: '
print numpy.asmatrix(c1).T

y = numpy.array([[1,0,5,4,5,6],[3,3,5,2,1,4]])
m2 = numpy.mean(y,1)
c2 = numpy.asmatrix(numpy.cov(y))
print 'mean 2: '
print m2
print 'cov 2: '
print c2
print 'cov 2 T: '
print numpy.asmatrix(c2).T


# Divergence 
print 'Trace of cov 1: '
print numpy.trace(c1)
print 'Trace of cov 2: '
print numpy.trace(c2)

print 'Inv of cov 1: '
c1i = numpy.linalg.inv(c1)
print numpy.linalg.inv(c1)
print 'Inv of cov 2: '
c2i = numpy.linalg.inv(c2)
print numpy.linalg.inv(c2)

print 'part 1 of divergence calculation: '
cDif = c1 - c2
cInvDif = c1i - c2i
div1 = 0.5 * numpy.trace(numpy.dot(cDif,cInvDif))

cInvSum = c1i + c2i
print cInvSum
mDif = numpy.matrix(m1 - m2).T
print 'dot product of cInvSum and mDif: '
dp1 = numpy.dot(cInvSum,mDif)
dp2 = numpy.dot(dp1,mDif.T)
print 'divergence part 2 calculation: '
div2 = 0.5 * numpy.trace(dp2)
print div2

print 'total divergence value: '
div = div1 + div2
print div

# Bhattacharyya
cAv = (c1 + c2)/2.0
print 'average cov matrix: '
print cAv
cAvInv = numpy.linalg.inv(cAv)
dp1 = numpy.dot(mDif.T,cAvInv) 
bh1 = 0.125 * numpy.dot(dp1,mDif)
print 'part 1 of BH dist: '
print bh1

cAvDet = numpy.linalg.det(cAv)
c1Det = numpy.linalg.det(c1)
c2Det = numpy.linalg.det(c2)

bh2 = 0.5 * numpy.log(cAvDet/((c1Det**0.5)*(c2Det**0.5)))
print 'part 2 of BH dist: '
print bh2

bh = bh1 + bh2
print 'Bhattacharyya distance: ', str(bh)

# JM 
jm = 2.0 * (1-numpy.exp(-bh))
print 'jm distance: ', str(jm)


# bh = 0.125 * (m1-m2) * numpy.linalg.inv(cAv) * (m1-m2).T + 0.5*numpy.log(numpy.abs(numpy.linalg.det(cAv/ scipy.linalg.sqrtm(c1*c2))));

# jm = 2.0 * (1-numpy.exp(-bh))


