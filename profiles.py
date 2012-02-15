# Code for working with surface-brightness profiles -- or any other
# two-column data file, where the first column is some kind of position -- in various ways.
#
# Includes code for unfolding a profile about a specified index

import copy
import numpy as N
import spline
#import spline2 as spline
import ellipse



def ReadProfile( inputFile, pix=1.0, ZP=None, skip=None ):
	"""Read in a two-column profile, returning a tuple of (x, y),
	where x is first column and y is second column, with all values converted
	to floating-point numpy arrays.  Blank lines and lines beginning with '#' are skipped.
		Optionally, the first N lines will be skipped, if skip=N is set.
		The x-values are multiplied by pix, which can be used to define
	the arcsec/pixel scale. (The default value is 1.0, which leaves the x-values
	unchanged.)
		If ZP is specified, then it is treated as a magnitude zero point and the
	y-values become ZP - 2.5*N.log10(y).
	"""

	lines = open(inputFile).readlines()
	if skip is not None:
		lines = lines[skip:]
	datalines = [ line for line in lines if line[0] != "#" and len(line.strip()) > 0 ]
	x = [ line.split()[0] for line in datalines ]
	y = [ line.split()[1] for line in datalines ]
	x = pix * N.array([ float(i) for i in x ])
	y = N.array([ float(i) for i in y ])
	if ZP is not None:
		y = ZP - 2.5*N.log10(y)
	return (x, y)


def WriteProfile( x, y, outputFilename, header=None, errs=None, mask=False,
					maskIndices=None ):
	"""Take two input vectors x and y (integer or float) and write them
	to a text file:
		   x    y

	Optionally, a header line can be supplied; it will become the first line
	of the file.

	If errs = a vector, it will be written as the third column

	If mask == True, then a mask vector will be generated and stored
	(mask value = 1 for data with values of "nan" or "inf").
	In addition, a list of indices to be masked can be passed in via maskIndices.
	"""

	nX = len(x)
	nY = len(y)
	if (nX != nY):
		msg = "WARNING: number of elements in x (%d) not same as number" % nX
		msg += " of elements in y (%d)!\nNothing saved.\n" % nY
		print msg
		return
	if (errs is not None) and (len(errs) != nX):
		msg = "WARNING: number of elements in error vector (%d) not same as number" % len(errs)
		msg += " of elements in x (%d)!\nNothing saved.\n" % nX
		print msg
		return
	nPts = nX
	if (mask is True) and (errs is None):
		# make a fake error vector
		errs = N.zeros(nPts) + 1.0
		if maskIndices is None:
			maskIndices = []

	outf = open(outputFilename, 'w')
	if header is not None:
		outf.write(header + "\n")
	for i in range(nPts):
		line = "%g\t\t%g" % (x[i], y[i])
		if errs is not None:
			line += "\t\t%g" % errs[i]
			if mask is True:
				if (i in maskIndices) or (N.isnan(y[i])) or (N.isinf(y[i])):
					thisMask = 1
				else:
					# this value is OK
					thisMask = 0
				line += "\t\t%d" % thisMask
		outf.write(line + "\n")
	outf.close()
	return


def MakeAvg( vector, index ):
	if index == 0:
		return vector[1]
	elif index == len(vector) - 1:
		return vector[-2]
	else:
		return 0.5*(vector[index - 1] + vector[index + 1])

def AvgBadVals( vector ):
	"""Returns a copy of input vector with values = "nan" or "inf" replaced by
	average of neighboring points.
	"""
	nPts = len(vector)
	badPts = [i for i in range(nPts) if (N.isnan(vector[i]) or N.isinf(vector[i]))]
	newVect = copy.copy(vector)
	for i in badPts:
		newVect[i] = MakeAvg(newVect, i)
	return newVect



def SimpleFold( vector, exclude=None ):
	"""Very simple function to fold a vector about its midpoint.  If the
	vector has an even number of points, then the midpoint is assumed to
	lie halfway between the middle two points, otherwise the middle point
	of the vector is the pivot point (and becomes the first point of the
	folded vector).

	Optionally, a specific region can be excluded by supplying a two-element
	list with the start and stop indices of the excluded region;
	as long as these do not include the midpoint of the vector, then the
	specified region will be excluded (only the corresponding points on the
	other side of the center will be used).

	Returns: tuple of (x-array, folded_input), where x-array is a corresponding
	vector of indices (starting at 0 for odd-number input profiles, or at 0.5
	for even-numbered input profiles).  Both returned vectors are numpy arrays.
	"""

	Npts = len(vector)
	midPoint = Npts / 2

	if exclude is not None:
		if (exclude[0] < midPoint) and (exclude[1] > midPoint):
			msg = "ERROR: start and stop indices of the excluded region "
			msg += " (%d, %d) bracket the profile's " % (exclude[0], exclude[1])
			msg += "midpoint (%d)!" % midPoint
			print msg
			return (None, None)

	if (Npts % 2) == 0:
		# even number of points in input
		rhs_start = midPoint
		rhs = N.array(vector[rhs_start:])
		lhs = list(vector[0:rhs_start])
		lhs.reverse()
		lhs = N.array(lhs)
		if exclude is not None:
			x1 = abs(exclude[0] - midPoint)
			x2 = abs(exclude[1] - midPoint)
			print x1, x2
			if exclude[0] > midPoint:
				replacement = lhs[x1:x2]
				rhs[x1:x2] = replacement
			else:
				replacement = rhs[x2:x1]
				lhs[x2:x1] = replacement
		foldedVector = 0.5*(lhs + rhs)
		outputIndices = N.arange(0.0, len(foldedVector)) + 0.5
	else:
		# odd number of points in input
		rhs_start = midPoint + 1
		rhs = N.array(vector[rhs_start:])
		lhs = list(vector[0:rhs_start - 1])
		lhs.reverse()
		lhs = N.array(lhs)
		if exclude is not None:
			# NOTE: THIS BIT OF CODE HAS NOT BEEN PROPERLY TESTED!
			x1 = abs(exclude[0] - midPoint)
			x2 = abs(exclude[1] - midPoint)
			print x1, x2
			if exclude[0] > midPoint:
				replacement = lhs[x1:x2]
				rhs[x1:x2] = replacement
			else:
				replacement = rhs[x2:x1]
				lhs[x2:x1] = replacement
		mainFoldedVector = 0.5*(lhs + rhs)
		foldedVector = [vector[midPoint]]
		foldedVector.extend(mainFoldedVector)
		outputIndices = N.arange(0.0, len(foldedVector))

	return (outputIndices, N.array(foldedVector))



def FoldProfileIDL( x, y, midIndex, offset ):
	"""Fold a profile about a user-specified midpoint.  User must specify
	the midpoint as follows: index value of the closest *leftmost* point to
	the midpoint, plus offset between that point and the midpoint.

	Returns a tuple of (x_folded, y_folded)

	Spline interpolation is used to regrid the profile so the specified
	midpoint value becomes a pixel center.  Then SimpleFold() is called to
	generate the final folded profile.

	Example: a vector of 100 points, index values starting at 0.0 [normal case,
	where the desired midpoint as at 35.7:
		foldeddata = FoldProfileIDL(x, y, 36, 0.7)
	because [36] = 35.0 and thus x[36] + 0.7 = 35.7
	"""

	Npts = len(x)
	# define vector of indices with offset
	offset_x = N.array(x) + offset
	# generate new y-values via spline interpolation
	spline_func = spline.Spline(x, y)
	y_new = [ spline_func(xx) for xx in offset_x ]

	print midIndex, Npts, Npts/2
	# Now extract a symmetric subvector centered on the new midpoint:
	if (midIndex < (Npts/2)):
		startIndex = 0
		endIndex = 2*midIndex
	else:
		startIndex = Npts - 2*(Npts - 1 - midIndex)
		endIndex = Npts - 2
	print startIndex, endIndex
	# note: we add 1 to the endIndex value 'cause this is Python
	newx = offset_x[startIndex:endIndex + 1]
	newy = y_new[startIndex:endIndex + 1]

	# Now fold the subvector
	folded_x, folded_y = SimpleFold(newy)

	return SimpleFold(newy)



def FoldProfile( x, y, foldIndex, pivot=False ):
	"""Fold a profile about foldIndex, such that the left half of
	the profile is y[0:foldIndex] and the right half is y[foldIndex + 1:].

	If pivot = True, then the foldIndex point is taken to be the unique pivot
	(r=0 point), so that left half = y[0:foldIndex + 1] and right half
	is y[foldIndex:].
	"""

	if pivot is True:
		leftEnd = foldIndex
		rightStart = foldIndex + 1
	else:
		leftEnd = foldIndex
		rightStart = foldIndex
	print leftEnd, rightStart
	yy = copy.copy(y)
	yy = yy[0:leftEnd]
	yy.reverse()
	y_left = N.array(yy)
	y_right = N.array(y[rightStart:])
	print len(y_left), len(y_right)
	y_folded = 0.5*(y_left + y_right)
	if pivot is True:
		extended = []
		extended.append(y[foldIndex])
		for yval in y_folded:
			extended.append(yval)
		y_folded = N.array(extended)

	return y_folded



def ReadAndFoldProfile( fileName, midIndex=None, pix=1.0, ZP=None ):
	"""Read in a profile from specified file, folding it about the
	midpoint.
		If midIndex is None (the default), then we fold about the
	middle of the profile using SimpleFold().  If midIndex is specified,
	then it must be a 1-based location; this is then fed to FoldProfileIDL(),
	which uses spline interpolation.
		Returns tuple of (radii, intensities), both of which will be
	numpy arrays.  If ZP is specified, then intensities --> magnitudes.
	"""

	(rr, ii) = ReadProfile(fileName)

	if midIndex is None:
		(rfold, ifold) = SimpleFold(ii)
	else:
		(offset, midIndex_floor) = math.modf(midIndex)
		midIndex_int = int(midIndex_floor) - 1
		nPts_orig = len(ii)
		indices = range(nPts_orig)
		ifold = FoldProfileIDL(indices, ii, midIndex_int, offset)
		rfold = N.arange(len(ifold))
		print "Non-simple folding not yet implemented!"

	if ZP is not None:
		ifold = ZP - 2.5*N.log10(ifold)
	return (pix*rfold, ifold)




def MergeTwoProfiles( r1, y1, r2, y2, r ):
	"""Produce a combined profile consisting of a radius vector and a corresponding
	vector of some other quantity (e.g., intensity, ellipticity), given two input
	profiles (r1 and y1, r2 and y2) and a designated transition radius r.
	We assume the vectors are NumPy arrays!
	"""

	# check for bad inputs
	if (r < r1[0]) or (r > r1[-1]):
		print "Requested transition radius (%g) is outside boundaries of r1 (%g--%g)!" % (r,
					r1[0], r1[-1])
		return None
	if (r < r2[0]) or (r > r2[-1]):
		print "Requested transition radius (%g) is outside boundaries of r2 (%g--%g)!" % (r,
					r2[0], r2[-1])
		return None
	r1_border = ellipse.NearestIndex(r1, r, noprint=True)
	r2_border = ellipse.NearestIndex(r2, r, noprint=True)
	end1 = r1_border[1]
	start2 = r2_border[1]
	if (r2[start2] < r1[end1]):
		start2 += 1

	newR = N.concatenate((r1[0:end1], r2[start2:]))
	newY = N.concatenate((y1[0:end1], y2[start2:]))

	return (newR, newY)



def MirrorReplace( y, xcenter, xstart, xend ):
	"""Given an input vector y, replace values y[xstart:xend] with the
	corresponding values on the other side of the profile, assuming the
	profile is approximately symmetric about the center at xcenter.

	xend MUST be the final point to be replaced, NOT the Python-style
	n+1 slice limit (i.e., to replace points 45 through 52, specify xstart=45
	and xend=52, NOT xend=53).

	Returns a copy of the input y-values, with the bad region replaced
	by its mirror-image segment from the other side of the profile.
	"""

	# check for bad input
	npts = len(y)
	xmin, xmax = 0, npts
	if (xcenter <= xmin) or (xcenter >= xmax):
		print("WARNING: xcenter is outside range of [0, length of y]!")
		return None
	if (xstart <= xmin) or (xstart >= xmax):
		print("WARNING: xstart is outside range of [0, length of y]!")
		return None
	if (xend <= xmin) or (xend >= xmax):
		print("WARNING: xend (%d) is outside range of [0, length of y]!")
		return None
	if (xstart < xcenter) and (xend > xcenter):
		print("WARNING: requested replacement zone includes center of profile!")
		return None
	if (xstart > xend):
		print("WARNING: xstart > xend!")
		return None
	if (xstart > xend):
		print("WARNING: xstart > xend!")
		return None

	# which side of profile are we on
	if (xstart < xcenter):
		leftSide = True
	else:
		leftSide = False

	# assuming that xcenter = integer
	d1, d2 = abs(xcenter - xstart), abs(xcenter - xend)
	nearDist = min([d1, d2])
	farDist = max([d1, d2])
	if leftSide:
		xrep_start = xcenter + nearDist
		xrep_end = xcenter + farDist
	else:
		xrep_start = xcenter - farDist
		xrep_end = xcenter - nearDist

	print("Replacing y[%d--%d] with y[%d--%d]" % (xstart, xend, xrep_start, xrep_end))
	# final check for unusable values
	if (xrep_start < 0) or (xrep_end >= npts):
		print("WARNING: replacement vector would be at least partly outside input vector!")
		return None

	ynew = list(y)
	goodDataSegment = ynew[xrep_start:xrep_end + 1]
	goodDataSegment.reverse()
	ynew[xstart:xend + 1] = goodDataSegment

	if type(y) is not list:
		ynew = N.array(ynew)
	return ynew



