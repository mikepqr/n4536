# Some utility functions for plotting with matplotlab

import subprocess, inspect

import numpy as N
from matplotlib import rc
import matplotlib.pyplot as plt
import pyfits

#import astro_utils
import profiles


PATH_TO_EPSTOPDF = "/usr/texbin/epstopdf"

def addmag( mag1, mag2 ):
	"""Takes two magnitudes (or numpy vectors of magnitudes) and
	returns their "sum" -- that is, the magnitude corresponding to
	the flux of the first value + the flux of the second value.
	"""

	zp = 30.0   # nominal zero point for conversions
	f1 = 10**(0.4*(zp - N.array(mag1)))
	f2 = 10**(0.4*(zp - N.array(mag2)))
	return zp - 2.5*N.log10(f1 + f2)


def StartTex():
	rc("text", usetex=True)
	rc("ps", usedistiller="ghostscript")
	return


def SaveEPSPDF( filePath ):
	"""Utility function to save a figure as both EPS and PDF in the same directory.
	"""
	
	# remove ".eps" if use supplied it
	filePath_base = filePath.split(".eps")[0]
	eps_filename = filePath_base + ".eps"
	pdf_filename = filePath_base + ".pdf"
	plt.savefig(eps_filename)
	plt.savefig(pdf_filename)

# following is older version written when matplotlib PDF generation wasn't
# working reliably
#
# def SaveEPSPDF( filePath ):
# 	"""Utility function to save a figure as EPS and then generate a PDF
# 	version in the same directory.
# 	"""
# 	
# 	# remove ".eps" if use supplied it
# 	filePath_base = filePath.split(".eps")[0]
# 	# save EPS file
# 	eps_filename = filePath_base + ".eps"
# 	plt.savefig(eps_filename)
# 	# spawn subprocess to convert EPS file to PDF
# 	cmd = PATH_TO_EPSTOPDF + " " + eps_filename + "\n"
# 	spawnedProc = subprocess.Popen(cmd, shell=True)
# 	spawnedProc.wait()



def ReadAndPlotProfile( profile_filename, color=None, returnVals=False ):
	"""Quick-and-dirty function to read a profile from a file and plot
	it (plots second column [y] vs first column [x]).  If returnVals=True,
	then the two vectors are returned as a tuple of (x, y).
	"""
	(x, y) = profiles.ReadProfile(profile_filename)
	if color is None:
		plt.plot(x, y)
	else:
		plt.plot(x, y, color=color)
	if returnVals:
		return (x, y)
	else:
		return None


def PlotProfilesWithSum( xVals, yList, fmtList=["b--", "r--", "g--", "m--", "y--", "b:"], fmtSum="k", 
						mag=True, yrange=None ):
	"""Given a vector of x-values and a list containing numpy vectors of y-values,
	plot y vs x for each vector in yList, and in addition plot the *sum* of
	the y-values.
		fmtList = list of matplotlib format specifications for each profile
		fmtSum = matplotlib format specification for the summed profile
		mag=True --> profiles in yList are in magnitudes (default)
	"""

	nX = len(xVals)
	nProfs = len(yList)
	if (nProfs > len(fmtList)):
		print "Warning: not enough formats in fmtList to match profiles in yList!"
		return None
	
	for i in range(nProfs):
		plt.plot(xVals, yList[i], fmtList[i])
	
	# make combined profile
	combinedVals = N.zeros(nX)
	if mag is True:
		for i in range(1, nProfs):
			combinedVals += addmag(yList[i - 1], yList[i])
	else:
		for i in range(nProfs):
			combinedVals += yList[i]
	
	# plot combined profile
	plt.plot(xVals, combinedVals, fmtSum)
	
	if (yrange is not None):
		plt.ylim(yrange[0], yrange[1])


def PlotFuncsWithSum( xVals, funcList, paramList, fmtList=["b--", "r--", "g--", "m--", "y--", "b:"],
						fmtSum="k", mag=True, xlog=False, ylog=False, yrange=None,
						funcNameList=None, plotSum=True, xOffset=0.0 ):
	"""Given a vector of x-values (xVals) and two lists containing functions (funcList) 
	and corresponding parameter vectors (paramList), plot func(x) vs x for each of the
	functions, and in addition plot the *sum* of the y-values.
	Functions in funcList must have the signature func(xVals, params), where params is
	a 1-D list of parameter values.
		xVals = vector of x values
		funcList = list or tuple of functions of the form func(xVals, params)
		paramList = list parameter vectors (i.e., a list of lists, or tuple of tuples)
		fmtList = list of matplotlib format specifications for each profile
		fmtSum = matplotlib format specification for the summed profile
		mag=True --> parameters in paramList are magnitudes, and output should be in
			magnitudes as well
		xlog=True --> plot vs log(xVals)
		yrange = optional list or tuple of (lower, upper) values for ylim()
		funcNameList = optional list of strings describing the components (used to assign
			labels to individual line plots, accessible to e.g. legend() function)
		plotSum -- if False, the sum of the individual functions is *not* plotted
		xOffset -- optional x-value offset to be added to xVals *for plotting purpose only*
			(i.e., the values of func(x)
	"""

	nX = len(xVals)
	nFuncs = len(funcList)
	if funcNameList is None:
		funcNameList = [None for i in range(nFuncs)]
	if (nFuncs > len(fmtList)):
		print "Warning: not enough formats in fmtList to match functions in funcList!"
		return None
	
	yValList = []
	for i in range(nFuncs):
		thisFunc = funcList[i]
		# check to see if this function has a "mag" keyword
		if "mag" in inspect.getargspec(thisFunc)[0]:
			yVals = funcList[i](xVals, paramList[i], mag=mag, magOutput=mag)
		else:
			yVals = funcList[i](xVals, paramList[i])
		if (xlog is True):
			plt.semilogx(xVals, yVals, fmtList[i], label=funcNameList[i])
		else:
			plt.plot(xVals, yVals, fmtList[i], label=funcNameList[i])
		yValList.append(yVals)
	
	if plotSum is True:
		# make combined profile
		combinedVals = N.zeros(nX)
		if mag:
			combinedVals = yValList[0]
			for i in range(1, nFuncs):
				combinedVals = addmag(combinedVals, yValList[i])
		else:
			for i in range(nFuncs):
				combinedVals += yValList[i]
	
		# plot combined profile
		if xlog is True and ylog is False:
			plt.semilogx(xVals, combinedVals, fmtSum, label="Sum")
		elif xlog is False and ylog is True:
			plt.semilogy(xVals, combinedVals, fmtSum, label="Sum")
		elif xlog is True and ylog is True:
			plt.loglog(xVals, combinedVals, fmtSum, label="Sum")
		else:
			plt.plot(xVals, combinedVals, fmtSum, label="Sum")
	
	if (yrange is not None):
		plt.ylim(yrange[0], yrange[1])


def GetFuncSum( xVals, funcList, paramList, mag=True ):
	"""Given a vector of x-values and two lists containing functions and corresponding
	parameter vectors, return the sum of the functions.
	Functions in funcList must have the signature func(xVals, params), where params is
	a 1-D list of parameter values.
		xVals = vector of x values
		funcList = list or tuple of functions of the form func(xVals, params)
		paramList = list parameter vectors (i.e., a list of lists, or tuple of tuples)
		mag=True --> parameters in paramList are magnitudes, and output should be in
			magnitudes as well
	"""

	nX = len(xVals)
	nFuncs = len(funcList)
	
	yValList = []
	for i in range(nFuncs):
		thisFunc = funcList[i]
		# check to see if this function has a "mag" keyword
		if "mag" in inspect.getargspec(thisFunc)[0]:
			yVals = funcList[i](xVals, paramList[i], mag=mag)
		else:
			yVals = funcList[i](xVals, paramList[i])
		yValList.append(yVals)
	
	# make combined profile
	combinedVals = N.zeros(nX)
	if mag:
		combinedVals = yValList[0]
		for i in range(1, nFuncs):
			combinedVals = addmag(combinedVals, yValList[i])
	else:
		for i in range(nFuncs):
			combinedVals += yValList[i]
	
	return combinedVals



def PlotFit( xVals, yDataVals, fitFuncList, fitParamList, x_fitted=None, y_fitted=None,
				fmtList=["b--", "r--", "g--", "m--"], fmtSum="k", mag=True, 
				xlog=False, ylog=False, xrange=None, yrange=None, residrange=None, ms=5,
				ylabel=None, filterName=None, title=None ):
	"""Plots data w/ fit, including a separate panel for residuals.
		x_fitted = range of x-values to plot as filled circles (actual data fitted)
		y_fitted = range of y-values to plot as filled circles (actual data fitted)
	"""
	
	nX = len(xVals)
	nFuncs = len(fitFuncList)
	if (nFuncs > len(fitParamList)):
		print "Warning: not enough parameter lists in fitParamList to match functions in fitFuncList!"
		return None
	if (nFuncs > len(fmtList)):
		print "Warning: not enough formats in fmtList to match functions in fitFuncList!"
		return None
	
	if mag is True and ylabel is None:
		if filterName is not None:
			ylabel1 = "$\\mu_{%s}$" % filterName
		else:
			ylabel1 = "$\\mu$"
	else:
		if filterName is not None:
			ylabel1 = "$I_{%s}$" % filterName
		else:
			ylabel1 = "$I$"
	ylabel2 = "$\\Delta$" + ylabel1
	
	
	# I. Plot data with fit
	plt.clf()
	ax1 = plt.axes([0.1, 0.35, 0.8, 0.55])
	
	# 1.a Plot data
	if (xlog is True):
		plt.semilogx(xVals, yDataVals, 'o', mfc='None', mec='k', ms=ms)
	else:
		plt.plot(xVals, yDataVals, 'o', mfc='None', mec='k', ms=ms)
	# 1.b. Overplot data used in fit as solid circles, if supplied by user
	if (x_fitted is not None) and (y_fitted is not None):
		plt.plot(x_fitted, y_fitted, 'ko', ms=ms)
	
	# 2.a Plot individual-function sub-profiles
	yValList = []
	for i in range(nFuncs):
		thisFunc = fitFuncList[i]
		# check to see if this function accepts a "mag" keyword
		if "mag" in inspect.getargspec(thisFunc)[0]:
			yVals = fitFuncList[i](xVals, fitParamList[i], mag=mag)
		else:
			yVals = fitFuncList[i](xVals, fitParamList[i])
		if (xlog is True):
			plt.semilogx(xVals, yVals, fmtList[i])
		else:
			plt.plot(xVals, yVals, fmtList[i])
		yValList.append(yVals)
	
	# 2.b Make & plot combined-functions profile
	combinedVals = GetFuncSum(xVals, fitFuncList, fitParamList, mag=mag)
#	combinedVals = N.zeros(nX)
#	if mag:
#		combinedVals = yValList[0]
#		for i in range(1, nFuncs):
#			combinedVals = addmag(combinedVals, yValList[i])
#	else:
#		for i in range(nFuncs):
#			combinedVals += yValList[i]
	if xlog is True and ylog is False:
		plt.semilogx(xVals, combinedVals, fmtSum)
	elif xlog is False and ylog is True:
		plt.semilogy(xVals, combinedVals, fmtSum)
	elif xlog is True and ylog is True:
		plt.loglog(xVals, combinedVals, fmtSum)
	else:
		plt.plot(xVals, combinedVals, fmtSum)

	
	if xrange is None:
		xrange = plt.xlim()
	plt.xlim(xrange[0], xrange[1])
	if (yrange is not None):
		plt.ylim(yrange[0], yrange[1])
	plt.ylabel(ylabel1)
	if title is not None:
		plt.title(title)
	# turn off x-tick labeling for this subplot
	ax1.set_xticklabels([])
	
	
	# II. Plot residuals
	ax2 = plt.axes([0.1, 0.1, 0.8, 0.25])
	
	# calculate residuals
	residVals = yDataVals - combinedVals
	# 1. Plot residuals
	if (xlog is True):
		plt.semilogx(xVals, residVals, 'o', mfc='None', mec='k', ms=ms)
	else:
		plt.plot(xVals, residVals, 'o', mfc='None', mec='k', ms=ms)
	plt.plot([xrange[0], xrange[1]], [0,0], 'k')

	# 2. Overplot residuals for fitted data (if latter supplied by user) as solid circles
	if (x_fitted is not None) and (y_fitted is not None):
		combinedVals_fitted = GetFuncSum(x_fitted, fitFuncList, fitParamList, mag=mag)
		residVals_fitted = y_fitted - combinedVals_fitted
		plt.plot(x_fitted, residVals_fitted, 'ko', ms=ms)

	plt.xlim(xrange[0], xrange[1])
	if (residrange is not None):
		plt.ylim(residrange[0], residrange[1])

	plt.xlabel("Radius [arcsec]")
	plt.ylabel(ylabel2)
	
	
	# III. Final cleanups:
	# somewhat kludgy way of suppressing the bottom-most y-axis tick label on the upper plot
	ylocs1 = ax1.get_yticks()
	print ylocs1
	new_labels = [str(yloc) for yloc in ylocs1]
	new_labels[-1] = ""
	ax1.set_yticklabels(new_labels)
	# same for the lower plot (turn off display of *last* (topmost) y-axis tick label)
	ylocs2 = ax2.get_yticks()
	new_labels2 = [str(yloc) for yloc in ylocs2]
	new_labels2[-1] = ""
	new_labels2[-2] = ""
	ax2.set_yticklabels(new_labels2)

	

def nicecont( imageData, xc=0, yc=0, width=None, height=None, levels=None, pix=1.0,
			axisLabel="pixels", title=None, imageExt=0, log=False ):
	"""Function which contour-plots an image.
		imageData = 2D Numpy array OR FITS image filename (image data is assumed to be
			in 0th header-data unit, unless imageExt is set to something else)
		xc, yc = optional center for axes (e.g., center of galaxy)
		width, height = width and height of subimage (centered on xc,yc) to be plotted;
			if height=None, then a square subimage of size width x width will be
			extracted
		levels = sequence (tuple, list, or Numpy array) of contour levels to be plotted
		pix = arcsec/pixel scale (for axis labeling)
		axisLabel = label for x and y axes
		title = title for plot
		imageExt = optional specification of a different header-data unit in input
			FITS image (if imageData points to a file)
		
	Example:
		>>> nicecont("image.fits", xc=202.4, yc=500.72, levels=N.arange(1.0, 20.0, 0.5))
	"""

	# handle case of user supplying a FITS filename
	if type(imageData) == str:
		hdulist = pyfits.open(imageData)
		imData = hdulist[imageExt].data
	else:
		imData = imageData
	
	if log is True:
		imData = N.log10(imData)
	
	ySize, xSize = imData.shape
	xPos = pix*(N.arange(1.0, xSize + 1.0) - xc)
	yPos = pix*(N.arange(1.0, ySize + 1.0) - yc)
	
	# extract centered subimage, if requested
	if width is not None:
		if height is None:
			height = width
		halfwidth = int(0.5*width)
		halfheight = int(0.5*height)
		x1 = xc - halfwidth
		x2 = xc + halfwidth
		y1 = yc - halfheight
		y2 = yc + halfheight
		xPos = xPos[x1:x2]
		yPos = yPos[y1:y2]
		imData = imData[y1:y2,x1:x2]
		
	
	if levels is not None:
		plt.contour(xPos, yPos, imData, levels, colors='k')
	else:
		plt.contour(xPos, yPos, imData, colors='k')
	plt.axes().set_aspect('equal')
	if axisLabel is not None:
		plt.xlabel(axisLabel)
		plt.ylabel(axisLabel)
	if title is not None:
		plt.title(title)



def HistPlot( xVals, nbins, xrange=None, offset=0.0, width=0.5, color=None ):
	if xrange is None:
		xrange = [min(xVals), max(xVals)]
	(counts, bin_edges) = N.histogram(xVals, nbins, range=xrange)
	midvals = N.zeros(nbins)
	for i in range(nbins):
		midvals[i] = 0.5*(bin_edges[i] + bin_edges[i + 1])
	
	plt.bar(midvals + offset, counts, width=width, color=color)



def DrawArrow( startCoord, endCoord, lineWidth=1, headWidth=8, headFrac=0.1,
                color="black", text="" ):
    """Draw an arrow from startCoords to endCoord (both should be 2-element
    tuples or lists of (x,y)), with lineWidth and headWidth in points (latter
    is width of arrowhead *base*).  headFrac specifies what fraction of the
    arrow the head occupies.
    text is optional text which will appear at the base of the arrow.
    """
    
    arrowProperties = dict(facecolor=color, width=lineWidth,
                        headwidth=headWidth, frac=headFrac)
    plt.annotate(text, xy=endCoord, xycoords="data", xytext=startCoord,
                arrowprops=arrowProperties)
    
    return None

def MultiArrows( xVals, yVals, dx, dy, **kwargs ):
	"""Plot arrows of same size and direction for a set of values
	
		xVals: list or array of x values
		yVals: list or array of y values
		dx: (scalar) offset in x for the arrow head
		dy: (scalar) offset in y for the arrow head
	
	Useful keyword arguments for DrawArrow:
		lineWidth=1 (points)
		headWidth=8 (points) -- width of arrow base
		headFrac=0.1 -- what fraction of total arrow is head
		color="black"
	"""
	nPts = len(xVals)
	if (nPts != len(yVals)):
		print("ERROR: xVals and yVals do not have the same length!")
		return None
	for i in range(nPts):
		startCoord = (xVals[i],yVals[i])
		endCoord = (xVals[i] + dx, yVals[i] + dy)
		DrawArrow(startCoord, endCoord, **kwargs)


# def MultiArrows2( xVals, yVals, dx, dy, **kwargs ):
# 	"""Plot arrows of same size and direction for a set of values
# 	
# 		xVals: list or array of x values
# 		yVals: list or array of y values
# 		dx: (scalar) offset in x for the arrow head
# 		dy: (scalar) offset in y for the arrow head
# 		
# 	Useful keyword parameters (passed on to matplotlib.axes.arrow and then
# 	to matplotlib.patches.FancyArrow):
# 		width -- defaults = 0.001
# 		head_width -- defaults to 3*width
# 		head_length -- defaults to 1.5*head_width
# 		overhang -- distance that the arrow is swept back (0 overhang means
#               triangular shape)
# 		length_includes_head: boolean (default=False) -- *True* if head is counted 
# 			in calculating the length.
# 		head_starts_at_zero: boolean (default=False) -- If *True*, the head
# 			starts being drawn at coordinate 0 instead of ending at coordinate 0.
# 		shape = ['full', 'left', 'right'] (default='full') -- arrow is full & symmetric,
# 			left-half-only, or right-half-only
# 		
# 	Other useful, more general keywords:
# 		color, edgecolor/ec, facecolor/fc
# 		
# 	Arrows are drawn from xVals[i],yVals[i] to xVals[i] + dx,yVals[i] + dy
# 	
# 	NOTE: this works, but seems to produce funny variation
# 	"""
# 	
# 	nPts = len(xVals)
# 	if (nPts != len(yVals)):
# 		print("ERROR: xVals and yVals do not have the same length!")
# 		return None
# 	for i in range(nPts):
# 		plt.arrow(xVals[i], yVals[i], dx, dy, **kwargs)
# 
