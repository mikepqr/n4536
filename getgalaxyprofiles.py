#!/usr/bin/env python
#
# Desired interface for profile extraction from galaxy image
# getprofile image-name radius =  Extract major-axis profile, using galaxy center, with length = 2*radius
# getprofile image-name  PA  = extract profile along specified *sky* PA
# getprofile image-name  PA --imagepa  =  extract profile along specfied image PA
#
# Useful/necessary files:
# galaxyinfo -- major-axis PA
# galaxycenter
# telescopePA
#
# pvect inputs
#    pvector.image = "filename.fits"
#    pvector.xc
#    pvector.yc
#    (pvector.width
#    (pvector.theta
#    (pvector.length
#    (pvector.vec_output = "output_filename.dat"
#    (pvector.out_type = "text"
#
#         image = "n3412rss.fit"  image containing vector to be plotted
#            x1 = 2170.14         x-coord of first point
#            y1 = 2431.79         y-coord of first point
#            x2 =                 x-coord of second point
#            y2 =                 y-coord of second point
#            xc = 1153.48         x-coord of center point
#            yc = 721.5           y-coord of center point
#        (width = 1)              width of strip
#        (theta = INDEF)          angle of vector (ccw from +x axis)
#       (length = INDEF)          length of vector in theta mode
#     (boundary = "constant")     type of boundary extension to use
#     (constant = 0.)             the constant for constant-valued boundary extension
#   (vec_output = "")             file or image name if output vector is desired
#     (out_type = "text")         type of output format (image|text)
#          (wx1 = 0.)             left user x-coord if not autoscaling
#          (wx2 = 0.)             right user x-coord if not autoscaling
#          (wy1 = 0.)             lower user y-coord if not autoscaling
#          (wy2 = 0.)             upper user y-coord if not autoscaling
#    (pointmode = no)             plot points instead of lines
#       (marker = "box")          point marker character
#     (szmarker = 0.005)          marker size (0 for list input)
#         (logx = no)             log scale x-axis
#         (logy = no)             log scale y-axis
#       (xlabel = "")             x-axis label
#       (ylabel = "")             y-axis label
#        (title = "imtitle")      title for plot
#          (vx1 = 0.)             left limit of device window (ndc coords)
#          (vx2 = 0.)             right limit of device window (ndc coords)
#          (vy1 = 0.)             bottom limit of device window (ndc coords)
#          (vy2 = 0.)             upper limit of device window (ndc coords)
#        (majrx = 5)              number of major divisions along x grid
#        (minrx = 5)              number of minor divisions along x grid
#        (majry = 5)              number of major divisions along y grid
#        (minry = 5)              number of minor divisions along y grid
#        (round = no)             round axes to nice values
#         (fill = yes)            fill device viewport regardless of aspect ratio?
#       (append = no)             append to existing plot
#



import sys, optparse, os, datetime, math
from pyraf import iraf
import angles
import profiles


TELESCOPE_PA_FILENAME = "telescope_pa.dat"
INFO_ROOTNAME = "galaxyinfo_"
HEADER_LINE1 = "# Profile generated via iraf task pvector (%s):\n"

class DummyOptions( object ):
	pass



def GetOffsetPosition( options, majorAxisPA, r, mode="parallel" ):
	"""
	Calculate pixel position for a position offset along the specified angle
	(options.imagePA) from the galaxy center (options.xCenter, options.yCenter),
	with the distance of the offset specified by r [in pixels].

	NOTE: although the code implies that options.imagePA = "major axis", this can be any
	position angle on the image, as long as it is relative to the image +y axis.
	"""

	x0, y0 = options.xCenter, options.yCenter


	# convert to angle relative to +x axis
	majorAxisPA_xref = majorAxisPA + 90.0
	if majorAxisPA_xref >= 360.0:
		majorAxisPA_xref -= 360.0
	print "majorAxisPA = %g, majorAxisPA_xref = %g" % (majorAxisPA, majorAxisPA_xref)
	angle_rad = math.radians(majorAxisPA_xref)
	if (mode == "parallel"):
		x_off = r * math.sin(angle_rad)
		y_off = r * math.cos(angle_rad)
	else:   # "perpendicular" offset
		x_off = r * math.cos(angle_rad)
		y_off = r * math.sin(angle_rad)
	print "x_off = %g, y_off = %g" % (x_off, y_off)
	return (x_off, y_off)



def GetDataLines( filename ):
    if os.path.isfile(filename): 
        lines = open(filename).readlines()
    else:
        return None
	dlines = [ line for line in lines if line[0] != "#" and len(line) > 1 ]
	return dlines



def GetGalaxyInfo(infoFileName):
	"""Extract coordinates for galaxy center (x.y) from infoFileName,
	given the image name.  Assumes that positionFilename has lines formatted
	like this:
	x_coord   y_coord   PA   telPA
	"""

	dlines = GetDataLines(infoFileName)

	if dlines == None:
		return None, None, 0., 0.

	if len(dlines) == 1:
		pp = dlines[0].split()
		if len(pp) == 2:
			return (float(pp[0]), float(pp[1]), 0., 0.)
		if len(pp) == 3:
			return (float(pp[0]), float(pp[1]), float(pp[2]), 0.)
		if len(pp) == 4:
			return (float(pp[0]), float(pp[1]), float(pp[2]), float(pp[3]))

	# only get here if we can't find image filename
	print "\n*** GetGalaxyInfo: can't find info for image file \"%s\"!\n" % infoFileName
	return None, None, 0., 0.


def GetTelescopePA():
	"""This version of GetTelescopePA reads the information in telescope_pa.dat,
	then stores the PA in a dictionary with the *filename* [rather than the filters,
	as in the sdssproc.py version of this function] as the key.
	"""

	try:
		dlines = [ line for line in open(TELESCOPE_PA_FILENAME) if line[0] != "#" ]
		paDict = {}
		for line in dlines:
			pp = line.split()
			fname = pp[0]
			pa = 0.5*(float(pp[2]) + float(pp[3]))
			paDict[fname] = pa
	except IOError:
		# file doesn't exist; assume that all images have standard orientation
		paDict = None
	return paDict



def AnnotateProfileFile( outputProfileFilename, imageName, xCenter, yCenter, width, theta, length ):
	"""Add a short header to a profile generated by pvect
	"""
	lines = open(outputProfileFilename).readlines()
	headerLines = []
	currentDateTime = datetime.datetime.now()
	newLine = HEADER_LINE1 % str(currentDateTime)
	headerLines.append(newLine)
	newLine = "# pvect %s xc=%f yc=%f width=%d theta=%f length=%d\n" % (imageName, xCenter,
				yCenter, width, theta, length)
	headerLines.append(newLine)
	outf = open(outputProfileFilename, 'w')
	for line in headerLines:
		outf.write(line)
	for line in lines:
		outf.write(line)
	outf.close()


def GetAndSaveProfiles( options, imageName, radius, outputRootName, xOffset=0.0, yOffset=0.0 ):

	# load IRAF packages
	iraf.plot(_doprint=0)

	# convert position angle to pvector format:
	posAng = options.imagePA + 90
	if (posAng >= 360.0):
		posAng -= 360.0

	# set up pvector
	iraf.pvector.theta = posAng
	iraf.pvector.length = 2*radius
	iraf.pvector.out_type = "text"

	xCenter = options.xCenter + xOffset
	yCenter = options.yCenter + yOffset

	outputFilenames = []
	for w in options.widths:
		outputName = "%s_w%g.dat" % (outputRootName, w)
		if (os.path.exists(outputName)):
			os.remove(outputName)
		iraf.pvector(imageName, xc=xCenter, yc=yCenter, width=w, vec_output=outputName)
		AnnotateProfileFile(outputName, imageName, xCenter, yCenter, w,	posAng, 2*radius)
		outputFilenames.append(outputName)

	return outputFilenames


def GetOneProfile( imageName, xCen, yCen, PA, width, radius, outputRootName="tempprofile",
					xOffset=0.0, yOffset=0.0, pixVal=1.0, fold=False ):
	"""Get pvector-derived profiles from an image.
		xCen = x position of profile center
		yCen = y position of profile center
		PA = position angle on image of profile

	Returns a tuple of (r, intensity)
	"""

	options = DummyOptions()
	options.widths = [width]
	options.imagePA = PA
	options.xCenter = xCen
	options.yCenter = yCen

	fnames = GetAndSaveProfiles(options, imageName, radius, outputRootName, xOffset, yOffset)
	fileName = fnames[0]
	if fold is False:
		x,y = profiles.ReadProfile(fileName, pix=pixVal)
	else:
		x,y = profiles.ReadAndFoldProfile(fileName, pix=pixVal)
	return x,y



def MakeOutputRoot( options, pa_str, shortName ):
	if options.outputName is None:
		outputRoot = "pvect_%s_%s" % (shortName, pa_str)
	else:
		outputRoot = "%s_%s" % (options.outputName, pa_str)
	return outputRoot



def main(argv=None):

	usageString = "%prog <image-filename>.fits <radius(pixels)> [options]\n"
	usageString += "\nAssumes the existence of galaxyinfo_<image-filename>.dat in the\n"
	usageString += "same directory. File should contain x0 y0 skyPA telPA on single line.\n"
	usageString += "x0, y0 compulsory. skyPA and telPA optional (assumed 0 if not present).\n"
	usageString += "telPA is output or IRAF north task. 0 = north up, east left, i.e. standard orientation\n"
	usageString += "x0, y0, skyPA and telPA may be overriden by command line\n"
	parser = optparse.OptionParser(usage=usageString, version="%prog ")

	parser.add_option("--output", "-o", type="str", dest="outputName", default=None,
						help="root name for output files (e.g. 'n1300rss'; default = 'pvect_<image-filename>_')")
	htext = "position angle of profile (on sky); default=taken from galaxyinfo_<image-filename>.dat"
	parser.add_option("--skyPA", type="float", dest="skyPA", default=None,
						help=htext)
	parser.add_option("--imagePA", type="float", dest="imagePA", default=None,
						help="position angle of profile (on image)")
	parser.add_option("--xc", type="float", dest="xCenter", default=None,
						help="x-coordinate of galaxy center (overrides galaxyinfo)")
	parser.add_option("--yc", type="float", dest="yCenter", default=None,
						help="y-coordinate of galaxy center (overrides galaxyinfo)")
	htext = "width of profile in pixels; can be used multiple times to build up list"
	htext += " (default is list of [1,3,5])"
	parser.add_option("--width", type="float", dest="widths", action="append",
						help=htext)
	parser.add_option("--perpendicular", action="store_true", dest="getPerpendicular", default=False,
						help="also extract perpendicular profile")
	parser.add_option("--parallel-offset", type="float", dest="parOffset", default=None,
						help="extract major-axis parallel profile at <parOffset> pixels along minor axis")
	parser.add_option("--perpendicular-offset", type="float", dest="perpOffset", default=None,
						help="extract minor-axis profile at <perpOffset> pixels along major axis")
	parser.add_option("--merge", action="store_true", dest="doMerge", default=False,
						help="merge w1, w3, w5 profiles [NOT IMPLEMENTED]")

	(options, args) = parser.parse_args(argv)
	# args[0] = name program was called with
	# args[1] = first actual argument, etc.

	# Regular mode (no input profile-specification files)
	if (len(args) < 3):
		print "You must supply an image filename and a radius!\n"
		return -1
	else:
		imageFilename = args[1]
		radius = int(args[2])

	galaxyName = os.path.splitext(imageFilename)[0]
	shortName = galaxyName

	infoFileName = INFO_ROOTNAME + galaxyName + ".dat"
	(x0, y0, galaxyPA, telPA) = GetGalaxyInfo(infoFileName)

	if (options.xCenter is None) and (options.yCenter is None):
		options.xCenter, options.yCenter = x0, y0

	# If no PA on command line
	if options.skyPA is None and options.imagePA is None:
		# user didn't specify a PA, so get sky PA from info file
		options.skyPA = galaxyPA

	if options.skyPA is not None:
		# using sky PA from file or user requested sky PA; convert to image PA
		options.imagePA = angles.SkyAngleToImageAngle(telPA, options.skyPA)

	# determine width or widths of profiles
	if options.widths is None:
		options.widths = [1.0, 3.0, 5.0]

	if options.skyPA is not None:
		paString = "pa%.1f" % options.skyPA
	else:
		paString = "imagepa%.1f" % options.imagePA
	outputRoot = MakeOutputRoot(options, paString, shortName)
	print "outputRoot = %s" % outputRoot

	if options.parOffset is not None:
		# calculate offset position
		majorAxisAngle = options.imagePA
		(xoff, yoff) = GetOffsetPosition(options, majorAxisAngle, options.parOffset)
		GetAndSaveProfiles(options, imageFilename, radius, outputRoot, xOffset=xoff, yOffset=yoff)
	elif options.perpOffset is not None:
		majorAxisAngle = options.imagePA
		options.imagePA += 90.0
		# calculate offset position
		(xoff, yoff) = GetOffsetPosition(options, majorAxisAngle, options.perpOffset, mode="perpendicular")
		GetAndSaveProfiles(options, imageFilename, radius, outputRoot, xOffset=xoff, yOffset=yoff)
	else:
		# standard, through-galaxy-center profile
		GetAndSaveProfiles(options, imageFilename, radius, outputRoot)

		if options.getPerpendicular:
			if options.skyPA is not None:
				paString = "pa%.1f" % (options.skyPA + 90.0)
			else:
				paString = "imagepa%.1f" % (options.imagePA + 90.0)
			outputRoot = MakeOutputRoot(options, paString, shortName)
			options.imagePA += 90.0
		GetAndSaveProfiles(options, imageFilename, radius, outputRoot)

	return 0


if __name__ == '__main__':

	main(sys.argv)
