import glob, os, numpy, time


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
tssM = numpy.mean(tssMean)
tssDev = ((numpy.matrix(tssMean - tssM)))
tssC = tssDev.T * tssDev / (N - 1)

chlM = numpy.mean(chlMean)
chlDev = ((numpy.matrix(chlMean-chlM)))
chlC = chlDev.T * chlDev / (N-1)

cdomM = numpy.mean(cdomMean)
cdomDev = ((numpy.matrix(cdomMean-cdomM)))
cdomC = cdomDev.T * cdomDev / (N-1)

time.sleep(5)

# create array to store 6-band combinations and distances
data = numpy.zeros([1,6])

# loop through every 6-band combination
bands = numpy.arange(0,60)
for b1 in bands:
   for b2 in range(b1+1,60):
      for b3 in range(b2+1,60):
         for b4 in range(b3+1,60):
            for b5 in range(b4+1,60):
               for b6 in range(b5+1,60): 

	          # create a filter with 1 at each band location, 0 elsewhere
                  testBands = numpy.zeros(60)
                  sixBands = [b1, b2, b3, b4, b5, b6]
		  testBands[sixBands] = 1
                  
                  # sample each class with the current band combination
                  # (multiply each class spectra by the filter array) 
                  tssMeanSampled = tssMean * testBands
                  chlMeanSampled = chlMean * testBands 
                  cdomMeanSampled = cdomMean * testBands 


