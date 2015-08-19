import glob
import math
import os
import numpy
import time 
import itertools
import numpy.linalg
import code


def bh_jm_div(m1, m2, c1, c2): 
 
   m1 = numpy.matrix(m1)
   m2 = numpy.matrix(m2) 
   c1 = numpy.asmatrix(c1)
   c2 = numpy.asmatrix(c2)

   # meanDif = m1 - m2
   # cAv = (c1 + c2) / 2.0

   # # Bhattacharyya computation
   # bh1 = 0.125 * (m1-m2) * numpy.linalg.inv(cAv) * (m1-m2).T 
   # bh2 = 0.5*numpy.log(numpy.abs(numpy.linalg.det(cAv/ scipy.linalg.sqrtm(c1*c2))))
   # bh = bh1 + bh2

   # # Jeffries-Matusita (JM) computation 
   # jm = 2.0 * (1-numpy.exp(-bh)) 

   # # divergence computation
   # c1i = numpy.linalg.inv(c1)
   # c2i = numpy.linalg.inv(c2)
   # div1 = 0.5 * numpy.trace((c1-c2)*(c1i-c2i)) 
   # div2 = 0.5 * numpy.trace((c1i+c2i)*(m1-m2).T*(m1-m2))
   # div = div1 + div2


   # Divergence 
   c1i = numpy.linalg.inv(c1)
   c2i = numpy.linalg.inv(c2)

   cDif = c1 - c2
   cInvDif = c1i - c2i
   div1 = 0.5 * numpy.trace(numpy.dot(cDif,cInvDif))

   cInvSum = c1i + c2i
   mDif = numpy.matrix(m1 - m2).T
   dp1 = numpy.dot(cInvSum,mDif)
   dp2 = numpy.dot(dp1,mDif.T)
   div2 = 0.5 * numpy.trace(dp2)
   div = div1 + div2
   #print 'divergence: ', str(div)

   # Bhattacharyya
   cAv = (c1 + c2)/2.0
   cAvInv = numpy.linalg.inv(cAv)
   dp1 = numpy.dot(mDif.T,cAvInv) 
   bh1 = 0.125 * numpy.dot(dp1,mDif)

   cAvDet = numpy.linalg.det(cAv)
   c1Det = numpy.linalg.det(c1)
   c2Det = numpy.linalg.det(c2)
   bh2 = 0.5 * numpy.log(cAvDet/((c1Det**0.5)*(c2Det**0.5)))
   bh = bh1 + bh2
   #print 'Bhattacharyya distance: ', str(bh)

   # JM 
   jm = 2.0 * (1-numpy.exp(-bh))
   #print 'jm distance: ', str(jm)

   return bh, jm, div



def nCr(n,r):
   f = math.factorial
   return f(n) / f(r) / f(n-r)


if __name__ == '__main__':

   startTime = time.clock()

   # for each constituent directory, create array and read in spectra
   # save ever-other data point yielding 60-sample spectra

   # Class 1: total suspended sediments/solids/materials
   os.chdir('tss_spectra/')
   tssAll = numpy.zeros((60,len(glob.glob('*.txt'))))
   tssFilenames = glob.glob('*.txt')
   for i in range(0,len(tssFilenames)):
      tssAll[:,i] = numpy.loadtxt(tssFilenames[i])[::2]
   m1 = numpy.mean(tssAll, 1)
   os.chdir('../')

   # Class 2: chlorophyll A
   os.chdir('chlA_spectra/')
   chlAll = numpy.zeros((60,len(glob.glob('*.txt'))))
   chlFilenames = glob.glob('*.txt')
   for i in range(0,len(chlFilenames)):
      chlAll[:,i] = numpy.loadtxt(chlFilenames[i])[::2]
   os.chdir('../')

   # Class 3: color-dissolved organic material 
   os.chdir('cdom_spectra/')
   cdomAll = numpy.zeros((60,len(glob.glob('*.txt'))))
   cdomFilenames = glob.glob('*.txt')
   for i in range(0,len(cdomFilenames)):
      cdomAll[:,i] = numpy.loadtxt(cdomFilenames[i])[::2]
   os.chdir('../')

   # define number of spectral samples
   wavelengths = numpy.loadtxt('wavelength.txt')[::2]
   N = len(wavelengths) 
   print 'Number of samples/wavelengths:', str(N)

   # compute mean vector per class
   m1 = numpy.mean(tssAll, 1)
   m2 = numpy.mean(chlAll, 1)
   m3 = numpy.mean(cdomAll, 1)

   # compute covariance matrix, C, for each class
   c1 = numpy.cov(tssAll)
   c2 = numpy.cov(chlAll)
   c3 = numpy.cov(cdomAll)

   # determine total number of band combinations 
   numberOfComputations =  nCr(60,6) * nCr(3,2)   
   print 'Number of Computations to perform: ', str(numberOfComputations)

   # create array to store band combinations (columns 0-5) and distances 
   # (columns 6-8) for every possible band combination (row dimension) and 
   # combination of 2 classes (band dimension. 3 choose 2 = 3)
   data = numpy.zeros([numberOfComputations,9,3])

   # loop through every 6-band combination
   bands = numpy.arange(0,60)

   r = 0 # current band combination counter 
   print 'Band combination ', str(r), ' of ', str(numberOfComputations)
   for b1 in bands:
      for b2 in range(b1+1,60):
         for b3 in range(b2+1,60):
            for b4 in range(b3+1,60):
               for b5 in range(b4+1,60):
                  for b6 in range(b5+1,60): 

                     #print 'Band combination ', str(r), ' of ', str(numberOfComputations)
                     # save the current 6-band combo in an the array
                     sixBands = [b1, b2, b3, b4, b5, b6]
                     data[r,0:6,0] = sixBands

   	               # subset the mean vectors based on current 6-band combination
                     m1Subset = m1[sixBands]
                     m2Subset = m2[sixBands]
                     m3Subset = m3[sixBands]

                     # subset covariance matrices
                     c1Subset = numpy.zeros([6,6])
                     c2Subset = numpy.zeros([6,6])
                     c3Subset = numpy.zeros([6,6])
                     for i in range(0,6):
                        for j in range(0,6):
                           c1Subset[i,j] = c1[sixBands[i],sixBands[j]]
                           c2Subset[i,j] = c2[sixBands[i],sixBands[j]]
                           c3Subset[i,j] = c3[sixBands[i],sixBands[j]]

                     # compute measures of seperability distances and store in array 
                     
                     try:
                        # Classes 1 & 2
                        data[r,6:9,0] = bh_jm_div(m1Subset,m2Subset,c1Subset,c2Subset)

                        # Classes 2 & 3
                        data[r,6:9,1] = bh_jm_div(m2Subset,m3Subset,c2Subset,c3Subset)

                        # Classes 1 & 3
                        data[r,6:9,2] = bh_jm_div(m1Subset,m3Subset,c1Subset,c3Subset)

                     except numpy.linalg.linalg.LinAlgError as err:
                        if 'Singular matrix' in err.message:
                           #print 'Singular matrix for band combination ', str(r),
                           #print ' No distances were calculated for it.'
                           pass
                        else:
                           pass
                        


                     #time.sleep(1)
                     r += 1
   elapsedTime = time.clock() - startTime
   print 'Elapsed time = %s [s]' % elapsedTime
   code.interact(local=locals())



                     



