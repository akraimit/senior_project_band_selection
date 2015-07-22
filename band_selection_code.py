import glob
import math
import os
import numpy
import time 
import itertools

def cov(vector1,vector2): 
   m1 = numpy.mean(vector1)
   m2 = numpy.mean(vector2)
   dev1 = ((numpy.matrix(vector1 - m1)))
   dev2 = ((numpy.matrix(vector2 - m2)))
   N = len(vector1)
   print 'cov method 1: '
   print dev1.T * dev2 / (N - 1)
   # return 


   X = numpy.column_stack([vector1, vector2])
   X -= X.mean(axis=0) 
   fact = N - 1 
   by_hand = numpy.dot(X.conj(), X.T) / fact
   print 'cov method 2'
   print by_hand

   print 'cov method 3, numpy'
   print numpy.cov(vector1,vector2)

   print 'cov method 4, stacking and transposing with numpy'
   z = numpy.vstack((vector1,vector2))
   print numpy.cov(z.T)


   time.sleep(5)

   return covariance


def jm_dist(m1, m2, c1, c2): 
   """ 
   title::
      jm_dist
   description::
      Computes the Jeffries-Matusita distance between a pair of classes. 
   attributes::
      m1, m2
         Mean vectors of the two classes, of type numpy array. 
      c1, c2
         Covariance matrices of the two classes, of type numpy matrix. 
   returns::
      dist - JM distance between two classes
   """
   m1 = numpy.matrix(m1)
   m2 = numpy.matrix(m2) 

   part1 = (m1 - m2).T * ((( c1 + c2 ) / 2) ** (-1)) * (m1 - m2)
   #part1 = (m1 - m2).T * numpy.linalg.inv((( c1 + c2 ) / 2)) * (m1 - m2)
   
   part2 = (1/2) * numpy.log( (1/2) * numpy.abs(c1+c2) / numpy.sqrt (numpy.cross(numpy.abs(c1),numpy.abs(c2))))
   
   a = (1/8) * part1 + part2 

   dist = numpy.sqrt(2 * (1-numpy.exp(-a)))

   return dist


def bhattacharyya_dist():

   return dist


def nCr(n,r):
   f = math.factorial
   return f(n) / f(r) / f(n-r)


if __name__ == '__main__':
   

   print 'TESTING: '
   m1 = numpy.array([3,6,4])
   m2 = numpy.array([7,12,9])
   print 'cov 1,2: '
   print cov(m1,m2)
   print 'cov 1,1: '
   c1 = cov(m1,m1)
   print c1

   print 'cov 2,2: '
   c2 = cov(m2,m2)
   print c2



   print nCr(60,6) * nCr(3,2)

   # read in spectral text files for 3 classes, and wavelengths
   chlMean = numpy.loadtxt('chlA_spectra_avg.txt')[::2]
   cdomMean = numpy.loadtxt('cdom_spectra_avg.txt')[::2]
   tssMean = numpy.loadtxt('tss_spectra_avg.txt')[::2]
   wavelengths = numpy.loadtxt('wavelength.txt')[::2]

   # for each constituent directory, create array and read in spectra
   # save ever-other data point yielding 60-sample spectra
   # also compute mean vector per class

   # total suspended solids
   os.chdir('tss_spectra/')
   tssAll = numpy.zeros((60,len(glob.glob('*.txt'))))
   tssFilenames = glob.glob('*.txt')
   for i in range(0,len(tssFilenames)):
      tssAll[:,i] = numpy.loadtxt(tssFilenames[i])[::2]
   tssMean = numpy.mean(tssAll, 1)
   os.chdir('../')

   # chlorophyll A
   os.chdir('chlA_spectra/')
   chlAll = numpy.zeros((60,len(glob.glob('*.txt'))))
   chlFilenames = glob.glob('*.txt')
   for i in range(0,len(chlFilenames)):
      chlAll[:,i] = numpy.loadtxt(chlFilenames[i])[::2]
   chlMean = numpy.mean(chlAll, 1)
   os.chdir('../')

   # color-dissolved organic material 
   os.chdir('cdom_spectra/')
   cdomAll = numpy.zeros((60,len(glob.glob('*.txt'))))
   cdomFilenames = glob.glob('*.txt')
   for i in range(0,len(cdomFilenames)):
      cdomAll[:,i] = numpy.loadtxt(cdomFilenames[i])[::2]
   cdomMean = numpy.mean(cdomAll, 1)
   os.chdir('../')

   # define number of spectral samples
   N = len(wavelengths) 

   # compute covariance matrix, C, for each class
   tssC = cov(tssMean, tssMean)
   chlC = cov(chlMean, chlMean)
   cdomC = cov(cdomMean, cdomMean)

   # create array to store 6-band combinations and distances
   print 'Number of Computations to perform: ', str(nCr(60,6) * nCr(3,2))

   numberOfComputations =  nCr(60,6) * nCr(3,2)

   # create array to store band combinations (columns 0-6) and distances 
   # (columns 6-7) for every possible band combination (row dimension) and 
   # combination of 2 classes (n choose k band dimension)
   data = numpy.zeros([numberOfComputations,8,15])

   # loop through every 6-band combination
   bands = numpy.arange(0,60)

   m1 = numpy.array([[2, 0, -9], [3, 4, 1]])
   m2 = numpy.array([[5, 2, 6], [-4, 4, 9]])
   c1 = cov(m1,m1)
   c2 = cov(m2,m2)
   jm = jm_dist(m1,m2,c1,c2)
   print jm


   for b1 in bands:
      for b2 in range(b1+1,60):
         for b3 in range(b2+1,60):
            for b4 in range(b3+1,60):
               for b5 in range(b4+1,60):
                  for b6 in range(b5+1,60): 

                     # save the current 6-band combo in an the array
                     sixBands = [b1, b2, b3, b4, b5, b6]
                     data[0,0:6] = sixBands
                     print data[0,:]




                     # loop throught all possible 2-band combinations, compute
                     # measures of separability



                     
   	               # subset the mean and cov per class based on current 
                     # 6-band combination. compute seperability distance metrics.
                     tssMeanSubset = tssMean[sixBands]
                     chlMeanSubset = chlMean[sixBands]
                     cdomMeanSubset = cdomMean[sixBands]

                     # FIGURE OUT HOW TO PULL THE CORRECT ELEMENTS OF THE COV MATRIX 
                     tssCovSubset = numpy.zeros([6,6])
                     chlCovSubset = numpy.zeros([6,6])
                     cdomCovSubset = numpy.zeros([6,6])
                     for i in range(0,6):
                        for j in range(0,6):
                           tssCovSubset[i,j] = tssC[sixBands[i],sixBands[j]]
                           chlCovSubset[i,j] = tssC[sixBands[i],sixBands[j]]
                           cdomCovSubset[i,j] = tssC[sixBands[i],sixBands[j]]




                     # Compute JM distance
                     data[6] = jm_dist()
                     



