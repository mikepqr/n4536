#!/usr/bin/env python
#
#    Python module to calculates mean sky value from a set of imexam-m
# measurements on an image, using the median values.
#

import sys, random, math
import numpy as N


def median( x, sorted=0 ):
	"""median( x ): computes median of values in list, tuple, or numpy array x.
	"""
	
	if len(x) < 1:
		print "median: length of x is < 1!"
		return None
		
	if not sorted:
		x_sorted = N.sort(x)
	else:
		x_sorted = x

	n = len(x_sorted)
	if (n % 2) == 0:   # Even number of elements:
		med = ( x_sorted[n/2 - 1] + x_sorted[n/2] )/2.0
	else:              # odd number of elements:
		med = x_sorted[n/2]
	return med


def stdev( x, correctBias=True ):
	"""stdev( list ): computes standard deviation of values in 
	list, tuple, or numpy array x.
	
	By default, we use Bessel's correction to compute the standard deviation
	as the square root of the unbiased estimator of the variance (dividing by
	N - 1, not by N), which corrects for using the estimated mean instead of
	the true (usually unknown) mean.  The user can specify computation of the
	"naive" standard deviation (division by N) by setting correctBias to False.
	"""
	
	n = len(x)
	if (n < 2):
		return 0.0
	x_bar = N.mean(x)
	sum = 0.0
	for i in range(n):
		diff = x[i] - x_bar
		sum = sum + diff*diff
	if correctBias is True:
		Ntot = n - 1
	else:
		Ntot = n
	standard_deviation = math.sqrt(sum / Ntot)
	return standard_deviation


def BootstrapError(origData, nRounds=200):
	"""Calculate error on sky background value using bootstrap technique."""

	random.seed()
	nPts = len(origData)
	orderedIndices = range(nPts)

	means = []
	for i in xrange(nRounds):
		# Generate new dataset:
		new_indices = [random.randint(0, nPts - 1) for i in orderedIndices]
		newData = [origData[i] for i in new_indices]
		means.append(N.mean(newData))

	return stdev(means)


	
def main(argv=None):

	doBootstrap = False
	nRounds = 1000
	
	if len(argv) < 2:
		print "You must supply an input filename.\n"
		sys.exit(1)
	
	fname = argv[1]
	if len(argv) > 2:
		if (argv[2] == "--bootstrap"):
			doBootstrap = True
	
	lines = open(fname, 'r').readlines()
	# ignore lines starting with "#"
	datalines = [ line for line in lines if line[0] != "#" and len(line.strip()) > 1 ]
	# extract median values (4th column) as floating-point values
	medianVals = [ float(dline.split()[3]) for dline in datalines ]
	
	print "%d values; mean = %f, median = %f" % (len(medianVals), N.mean(medianVals), median(medianVals))
	if doBootstrap:
		stdev_bootstrap = BootstrapError(medianVals, nRounds)
		print "standard deviation from %d bootstrap samples = %f\n" % (nRounds, stdev_bootstrap)


if __name__ == '__main__':
	
	main(sys.argv)
