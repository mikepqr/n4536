#!/usr/bin/env python

# Note that there is a unit-test file which accompanies this file;
# currently, it only tests SkyAngleToImageAngle and SphereFunc.


import math

RADIAN_TO_DEGREE = 180.0 / math.pi
DEGREE_TO_RADIAN = math.pi / 180.0


def RectifyPA( angle, maxAngle ):
	"""Convert angle to lie between 0 and maxAngle degrees."""
	if (angle < 0.0):
		return angle + maxAngle
	elif (angle >= maxAngle):
		return angle - maxAngle
	else:
		return angle



def SkyAngleToImageAngle( telescopePA, skyPA, flipped=False, symmetric=True ):
	"""Function which converts a position angle on the sky (degrees E of N)
	to a position angle on an image (degrees CCW from +y axis), given the
	"telescope PA", which is the PA on the sky of the image +y axis (e.g.,
	the output from the IRAF-STSDAS task "north" applied to an HST image,
	or the value of the "SPA" header in an SDSS image).

	If flipped=True, it corrects for images which are flipped about the y-axis
	(e.g., WIYN images).
	If symmetric=True, then values > 180 are converted to the equivalent
	0--180 values."""

	# First, correct for rotation, then for possible reflection about y-axis
	imagePA = skyPA - telescopePA
	if flipped:
		imagePA = 360 - imagePA

	# Now, correct for angles outside [0,360] and optionally for angles
	# outside [0,180]
	if imagePA < 0:
		imagePA += 360.0
	if imagePA >= 360.0:
		imagePA = imagePA - 360.0
	if symmetric:
		if imagePA >= 180.0:
			imagePA = imagePA - 180

	return imagePA



def ImageAngleToSkyAngle( telescopePA, imagePA, flipped=False, symmetric=True ):
	"""Function which converts a position angle on an image (assumed to be
	degrees CCW from +y axis) to position angle on the sky (degrees E of N),
	given the "telescope PA", which is the PA on the sky of the image +y axis
	(e.g., the output from the IRAF-STSDAS task "north" applied to an HST image,
	or the value of the "SPA" header in an SDSS image).

	If flipped=True, it corrects for images which are flipped about the y-axis
	(e.g., WIYN images).
	If symmetric=True, then values > 180 are converted to the equivalent
	0--180 values."""

	# First, correct for possible reflection about y-axis,
	# then for rotation
	if flipped:
		imagePA = 360 - imagePA
	skyPA = imagePA + telescopePA

	# Now, correct for angles outside [0,360] and optionally for angles
	# outside [0,180]
	if skyPA < 0:
		skyPA += 360.0
	elif skyPA >= 360.0:
		skyPA = skyPA - 360.0
	if symmetric:
		if skyPA >= 180.0:
			skyPA = skyPA - 180

	return skyPA



def SphereFunc(dec1, dec2):
	"""Function to calculate correction to RA, converting RA degrees to
	degrees on the sky, assuming a declination which is the average of
	the two specified declinations."""

	avgDec = (dec1 + dec2)/2.0
	scalingFactor = math.cos(avgDec * DEGREE_TO_RADIAN)
	#print "avgDec = %f, scaleFactor = %f" % (avgDec, scalingFactor)
	return scalingFactor


def DeltaRADec( coords1, coords2 ):
	"""Function to calculate distance between two coordinates (in degrees).
	Both coordinates are lists of (ra,dec) in degrees."""

	ra1 = coords1[0]
	dec1 = coords1[1]
	ra2 = coords2[0]
	dec2 = coords2[1]
	deltaRA = SphereFunc(dec1, dec2) * (ra1 - ra2)
	deltaDec = dec1 - dec2
	return math.sqrt( deltaRA*deltaRA + deltaDec*deltaDec )


def rpa( center, point ):
	"""Function to calculate radius and position angle, given two
	coordinates (e.g., galaxy center and some other point).
	Coordinates must be 2-element lists."""

	xP = point[0]
	yP = point[1]
	xdiff = xP - center[0]
	ydiff = yP - center[1]
	radius = math.sqrt(xdiff*xdiff + ydiff*ydiff)
	pa0 = math.asin(ydiff / radius) * RADIAN_TO_DEGREE
	if (xdiff < 0):
		pa =  90 - pa0
	else:
		pa = 270 + pa0

	return (radius, pa)


def ellipser( ellipticity, deltaPA ):
	"""Function to calculate length of radius vector to a point on an
	ellipse, given a user-specified angle relative to the ellipse major axis."""

	a = 1.0
	b = (1 - ellipticity)
	deltaPA_rad = deltaPA * DEGREE_TO_RADIAN
	sinepiece = a * math.sin(deltaPA_rad)
	cospiece = b * math.cos(deltaPA_rad)
	r = a*b / (math.sqrt( sinepiece*sinepiece + cospiece*cospiece ))
	return r


def ifrome( ellipticity, q = 0.2 ):
	"""Function to calculate the inclination of a galactic disk given
	an observed disk ellipticity (this assumes an intrinsic axis ratio
	q = c/a = 0.2, but this can be changed using the q keyword).
	Based on formula in Hubble (1926, ApJ 64: 321)."""

	denominator = 1.0 - q**2
	numerator = 2.0*ellipticity - ellipticity*ellipticity
	sini = math.sqrt(numerator / denominator)
	return (math.asin(sini) * RADIAN_TO_DEGREE)



# *** FUNCTIONS FOR PROJECTING AND DEPROJECTING ***

def deprojectv( v_obs, inclination ):
	"""Calculate deprojected velocity, given observed velocity v_obs
	and inclination in degrees."""

	return v_obs / math.sin( inclination * DEGREE_TO_RADIAN )



def projectpa( deltaPA, inclination ):
	"""Function to calculate a projected position angle, given an input
	(unprojected) position angle.  Position angles are relative to
	disk line-of-nodes."""

	deltaPA_rad = deltaPA * DEGREE_TO_RADIAN
	i_rad = inclination * DEGREE_TO_RADIAN
	deltaPA_proj = math.atan( math.tan(deltaPA_rad) * math.cos(i_rad) )
	return ( deltaPA_proj * RADIAN_TO_DEGREE )


def deprojectr( deltaPA, inclination, r ):
	"""Function to calculate a deprojected length, given an input
	observed position angle (*relative to disk line-of-nodes*, *not*
	straight position angle east of north!) and inclination, both in
	degrees, and an input observed (projected) length r.
	Returns the deprojected length."""

	deltaPA_rad = deltaPA * DEGREE_TO_RADIAN
	i_rad = inclination * DEGREE_TO_RADIAN
	cosi = math.cos(i_rad)
	sindp = math.sin(deltaPA_rad)
	cosdp = math.cos(deltaPA_rad)
	scale = math.sqrt( (sindp*sindp)/(cosi*cosi) + cosdp*cosdp )
	return ( scale * r )


def deprojectpa( deltaPA, inclination ):
	"""Function to calculate a deprojected position angle, given an input
	observed position angle (*relative to disk line-of-nodes*, *not*
	straight position angle east of north!) and an input inclination, both
	in degrees.  Returns the deprojected position angle, relative to disk
	line-of-nodes, in degrees."""

	deltaPA_rad = deltaPA * DEGREE_TO_RADIAN
	i_rad = inclination * DEGREE_TO_RADIAN
	deltaPA_deproj = math.atan( math.tan(deltaPA_rad) / math.cos(i_rad) )
	return ( deltaPA_deproj * RADIAN_TO_DEGREE )


def deprojectpa_abs( obsPA, diskPA, inclination, symmetric=True ):
	"""Function to calculate deprojected PA on the sky of a structure, as
	if the galaxy were rotated to face-on orientation.
		obsPA = observed (projected) PA of structure
		diskPA = major axis (line of nodes) of disk
		inclination = inclination of disk (0 = face-on).
	If symmetric=True, then values > 180 are converted to the equivalent
	0--180 values."""

	if (inclination == 0.0):
		dPA = obsPA
	else:
		deltaPA = obsPA - diskPA
		deltaPA_dp = deprojectpa(deltaPA, inclination)
		dPA = deltaPA_dp + diskPA
	if symmetric:
		dPA = RectifyPA(dPA, 180.0)
	return dPA


def minoraxis( structurePA, diskPA, inclination=0.0 ):
	"""Function to determine the PA of the minor axis of a given structure,
	optionally including the appropriate projection effects if the galaxy is
	inclined.  Returned value is constrained to lie in range 0 <= x < 180."""

	if (inclination == 0.0):
		minoraxis_proj = structurePA + 90.0
	else:
		structurePA_inplane = deprojectpa_abs( structurePA, diskPA, inclination )
		minorPA_inplane = RectifyPA(structurePA_inplane + 90.0, 180.0)
		minoraxis_proj = projectpa(minorPA_inplane - diskPA, inclination)
		minoraxis_proj += diskPA
		msg = "minoraxis: diskPA = %.1f, structPA = %.1f, inc = %.1f, " % (diskPA, structurePA, inclination)
		msg += "structurePA_inplane = %.1f, minorPA_inplane = %.1f" % (structurePA_inplane, minorPA_inplane)
		print msg

	return RectifyPA(minoraxis_proj, 180.0)


def tan2t_rad( a, b, u, i ):
	"""Function to calculate angle t_mm (angle of deprojected ellipse's semi-
	major or -minor axis, w.r.t. angle u).
		a,b = semi-major and -minor axes of *observed* ellipse
		u = PA of *observed* ellipse major axis w.r.t. line-of-nodes
		i = inclination angle of observed system
	*** inputs must be in radians! ***
	Formula from Seppo Laine (2003, private comm.)
	"""
	g = 1.0/math.cos(i)
	g2 = g*g
	a2 = a*a
	b2 = b*b
	numerator = a*b*(g2 - 1.0)*math.sin(2.0*u)
	denom = math.sin(u)**2 * (b2 - g2*a2) + math.cos(u)**2 * (g2*b2 - a2)
	t_mm = 0.5 * math.atan2(numerator,denom)
	return t_mm


def deprojectell_abraham( deltaPA, inclination, ell_obs ):
	"""Function to determine intrinsic (in-plane) ellipticity of an
	observed (projected) ellipse.
		deltaPA = PA of *observed* ellipse major axis w.r.t. line-of-nodes
		inclination = inclination assumed for deprojection (in degrees).
		ell_obs = *observed* ellipticity of (projected) ellipse
	Formulae from Abraham et al. (1999, MNRAS 308: 569, Eqn.2--3);
	appears to produce same results as deproject_laine().
	Tested July 2009.
	"""

	# Abraham et al. formulae use axis ratio q = b/a
	q_obs = 1 - ell_obs
	q_obs2 = q_obs * q_obs
	q_obs_inv2 = 1.0/q_obs2
	deltaPA_rad = deltaPA * DEGREE_TO_RADIAN
	i_rad = inclination * DEGREE_TO_RADIAN

	sini = math.sin(i_rad)
	sini_2 = sini * sini
	cosi = math.cos(i_rad)
	seci_2 = 1.0 / (cosi * cosi)
	sini_4 = sini_2 * sini_2
	sinPA = math.sin(deltaPA_rad)
	sinPA_2 = sinPA * sinPA
	cosPA = math.cos(deltaPA_rad)
	cosPA_2 = cosPA * cosPA

	# Abraham et al. Eqn. 3:
	X = seci_2 * ( 2.0*cosPA_2*sinPA_2*sini_4 + \
		q_obs2*(1.0 - sinPA_2*sini_2)**2 + \
		q_obs_inv2*(1.0 - cosPA_2*sini_2)**2 )

	# Abraham et al. Eqn. 2:
	# q0_2 is the intrinsic (deprojected) axis ratio, squared
	q0_2 = 0.5 * (X - math.sqrt(X*X - 4.0))

	return (1.0 - math.sqrt(q0_2))


def deprojectell_laine( a, ell, u_deg, i_deg, printAll=0 ):
	"""Function to deproject an ellipse.  Calculates new semi-major axis,
	ellipticity, and PA of ellipse major axis w.r.t. line of nodes.
		a = *observed* semi-major axis of (projected) ellipse
		ell = *observed* ellipticity of (projected) ellipse
		u_deg = PA of *observed* ellipse major axis w.r.t. line-of-nodes
		i_deg = inclination assumed for deprojection
	Returned value is deprojected ellipticity; if printAll=1, then the
	derived (deprojected) semi-major axis and position angle are also
	printed.
	Formulae and general algorithm from Seppo Laine (2003, private comm.);
	appears to produce same results as deproject_abraham().
	Tested July 2009.
	"""
	b = a * (1.0 - ell)
	i = i_deg * DEGREE_TO_RADIAN
	u = u_deg * DEGREE_TO_RADIAN
	g = 1.0/math.cos(i)

	# Get angle of deprojected major or minor axis w.r.t. u:
	t_mm = tan2t_rad(a, b, u, i)

	# Calculate semi-major and semi-minor axis lengths, first using t_mm:
	xp = a*math.cos(t_mm)*math.cos(u) + b*math.sin(t_mm)*math.sin(u)
	yp = g*a*math.cos(t_mm)*math.sin(u) - g*b*math.sin(t_mm)*math.cos(u)
	semiaxis_1 = math.sqrt(xp*xp + yp*yp)
	t_mm2 = t_mm + math.pi/2.0
	# and now using t_mm + pi/2:
	xp = a*math.cos(t_mm2)*math.cos(u) + b*math.sin(t_mm2)*math.sin(u)
	yp = g*a*math.cos(t_mm2)*math.sin(u) - g*b*math.sin(t_mm2)*math.cos(u)
	semiaxis_2 = math.sqrt(xp*xp + yp*yp)
# 	print "(t_mm, t_mm2) = (%f, %f) deg" % (t_mm*RADIAN_TO_DEGREE,
# 				t_mm2*RADIAN_TO_DEGREE)
	t_mm_list = [t_mm, t_mm2]

	# OK, now figure out *which* of those two axis lengths is the *major* axis:
	if (semiaxis_1 >= semiaxis_2):
		major_index = 0
		a_dp = semiaxis_1
		b_dp = semiaxis_2
	else:
		major_index = 1
		a_dp = semiaxis_2
		b_dp = semiaxis_1
	ell_dp = 1.0 - (b_dp/a_dp)
	if (printAll):
		print "Deprojected major axis =	%f,	angle =	%f deg,	ell	= %f" %	(a_dp,
			u_deg + 180 - t_mm_list[major_index]*RADIAN_TO_DEGREE, ell_dp)
	return ell_dp




def deprojectell_m95( barPA, diskPA, inclination, ell_obs ):
	"""Simple function to attempt deprojecting an ellipse using the crude
	approach of Martin (1995): deproject the semimajor and -minor axes
	independently, then use their deprojected ratio."""

	a_obs = 1.0
	b_obs = 1.0 - ell_obs

	deltaPA = math.fabs(barPA - diskPA)
	if deltaPA > 90:
		deltaPA = 180 - deltaPA
	a_dp = deprojectr(deltaPA, inclination, a_obs)

	if barPA >= 90:
		PAb = barPA - 90
	else:
		PAb = barPA + 90
	deltaPAb = math.fabs(PAb - diskPA)
	if deltaPAb > 90:
		deltaPAb = 180 - deltaPAb
	b_dp = deprojectr(deltaPAb, inclination, b_obs)

	if (b_dp < a_dp):
		return 1.0 - b_dp/a_dp
	else:
		return 1.0 - a_dp/b_dp



def deprojectOblateSpheroid( inclination, ell_obs ):
	"""Function to obtain the intrinsic, edge-on ellipticity of an oblate spheroid,
	given an observed inclination (in degrees) and ellipticity.  From Carollo &
	Danziger (1994), who are apparently quoting Mihalas & Binney (1981): "equations 5-21."
	"""

	i_rad = math.radians(inclination)
	sin2i = math.sin(i_rad)**2
	cot2i = 1.0/(math.tan(i_rad)**2)
	stuff = ((1.0 - ell_obs)**2)/sin2i - cot2i
	ell_intrinsic = 1.0 - math.sqrt(stuff)
	return ell_intrinsic


def projectOblateSpheroid( inclination, ell_in ):
	"""Function to obtain the observed ellipticity of an oblate spheroid having
	intrinsic (edge-on) ellipticity ell_in, assuming that it is observed with the
	specified inclination.  From Carollo & Danziger (1994), who are apparently
	quoting Mihalas & Binney (1981): "equations 5-21."
	"""

	i_rad = math.radians(inclination)
	sini2 = math.sin(i_rad)**2
	cosi2 = math.cos(i_rad)**2
	stuff = (1.0 - ell_in)**2 * sini2 + cosi2
	ell_obs = 1.0 - math.sqrt(stuff)
	return ell_obs



# Function to calculate apparent axis ratio of a thin
# elliptical ring with true axis ratio q0, observed with
# inclination i and orientation phi (where phi = angle
# between ellipse and line-of-nodes, in the disk plane).
# Also calculates the projected position angle of the
# ellipse's apparent major axis, relative to the line-of-nodes.
#    CURRENTLY SEEMS TO TREAT PHI AS 90 - PHI!!!
#    From Buta 1986 (ApJ Supp. 61: 609--630, eqns. 1 and 6,7).
def projectAxisRatio( q0, inclination, phi ):
	# Convert angles from degrees to radians:
	i_rad = inclination * DEGREE_TO_RADIAN
	phi_rad = phi * DEGREE_TO_RADIAN
	# Calculate squares and trig terms:
	q02 = q0 * q0
	sinphi = math.sin(phi_rad)
	sinphi2 = sinphi * sinphi
	cosphi = math.cos(phi_rad)
	cosphi2 = cosphi * cosphi
	cosi = math.cos(i_rad)
	cosi2 = cosi * cosi
	# Calculate j, l, k terms:
	j = q02 * sinphi2 * cosi2 + cosphi2 * cosi2
	k = (1.0 - q02) * sinphi * cosphi * cosi
	l = q02 * cosphi2 + sinphi2

	# Final calculation of projected axis ratio:
	bracketPiece = math.sqrt( (j - l)*(j - l) + 4.0*k*k )
	numerator = j + l - bracketPiece
	denominator = j + l + bracketPiece
	q_proj = math.sqrt( numerator/denominator )
	# Final calculation of projects position angle:
	projPA = 0.5 * math.atan(2*k / (l - j))

	return (q_proj, projPA * RADIAN_TO_DEGREE)


