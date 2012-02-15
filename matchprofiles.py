# Python code to determine best scaling and offset values to match to data
# vectors (e.g., two vectors of intensities with the same or overlapping radii).
# to match two profiles in the region rmin <= r <= rmax.  The first profile is
# assumed to be "correct

import math
import numpy as N
import scipy.optimize
import matplotlib.pyplot as p
#import spline2 as spline
import spline


def LeastsqFunc(p, i1, i2):
	"""Returns a NumPy array containing the differences between i1 (the
	calibrated fluxes, regridded to match i2 in radial spacing) and
	p[0]*(i2 - p[1]), where i2 is the array of uncalibrated fluxes, p[0] is
	the scale factor and p[1] is the offset (e.g., the residual sky level in
	i2)."""
	
	scaleFactor = p[0]
	offset = p[1]
	
	return i1 - scaleFactor*(i2 - offset)

	

def MatchProfiles( r1, i1, r2, i2, rmin, rmax, offset=False, magInput=False,
					zp1=None ):
	"""Procedure to determine best scaling and (optional) subtraction
	necessary to match two surface-brightness profiles in the region r_min <= r_max.
	The first profile [r1, i1] is assumed to be "correct" -- that is, it has no
	residual background-subtraction problems.  The second profile [r2, i2] is the
	"uncalibrated" one we want to match to the first by calculating the
	necessary scaling factor -- and optionally an additive offset to
	account for possible background-subtraction errors -- such that:
		i2_scaled = scaleFactor * (i2 + offsetValue)
	
	r1, i1 = radius and intensity vectors for the first ("correct"/reference)
		profile
	r2, i2 = radius and intensity vectors for the second profile
	rmin, rmax = radius values defining the overlap region (i.e., the region where
		the two profiles will be matched)
	
	If offset=True, then an additive offset (e.g., residual sky background
	correction) is calculated for the second profile
	
	If magInput=True, then profiles are assumed to be in magnitude units; they
	are internally converted to intensities, and the resulting scaleFactor
	is converted to a magnitude offset.
	
	zp1 = optional specification of zero point for r1,i1 profile; if supplied,
	the corresponding zero point for r2,i2 will be printed
	"""
	
	if (rmin < r1[0]) or (rmin < r2[0]):
		print("Error: rmin < minimum radius in r1 or r2!")
		return None
	if (rmax > r1[-1]) or (rmax > r2[-1]):
		print("Error: rmax > maximum radius in r1 or r2!")
		return None
	
	if (magInput is True):
		i1_flux = 10**( (30 - i1)/2.5 )
		i2_flux = 10**( (30 - i2)/2.5 )
	else:
		i1_flux = i1
		i2_flux = i2
	valid_indices = [ i for i in range(len(r2)) if r2[i] >= rmin and r2[i] <= rmax ]
	r2_imin = valid_indices[0]
	r2_imax = valid_indices[-1]
	uncalibrated_r = N.array(r2[r2_imin:r2_imax + 1])
	uncalibrated_flux = N.array(i2_flux[r2_imin:r2_imax + 1])
	npts = len(uncalibrated_r)
	
	spline_func = spline.Spline(N.array(r1), N.array(i1_flux))
	calibrated_flux_regrid = [ spline_func(uncalibrated_r[i]) for i in range(npts) ]
	calibrated_flux_regrid = N.array(calibrated_flux_regrid)
	
	initialScaling = N.mean(calibrated_flux_regrid / uncalibrated_flux)
	print("Simple scaling: scaled i2 = %g * i2" % initialScaling)
	zeropointOffset = -2.5*math.log10(initialScaling)
	print("(Add %f to i1 zero point to get effective zero point for raw i2)" % zeropointOffset)
	
	if offset is True:
		p0 = [initialScaling, 0.0]
		fluxData = (calibrated_flux_regrid, uncalibrated_flux)
		(pp, s) = scipy.optimize.leastsq(LeastsqFunc, p0, args=fluxData)
		print("Including offset: scaled i2 = %g * (i2 - %g)" % (pp[0], pp[1]))
		zeropointOffset = -2.5*math.log10(pp[0])
		print("(Add %f to i1 zero point to get effective zero point for raw i2)" % zeropointOffset)
		if (zp1 is not None):
			zp2 = zp1 + zeropointOffset
			print("ZP = %f for i2 profile" % zp2)
		return pp
	else:
		return initialScaling



def QuickMatchEfits( efit1, efit2, range=None, offset=False, plot=True, xlog=False, ylog=False ):
	"""Wrapper for QuickMatch(), allowing user to specify the two profiles as ellipse-fit
	objects (as produced by ellipse.Readellipse).  See QuickMatch() help for more info.
	"""
	
	r1 = efit1['sma']
	r2 = efit2['sma']
	i1 = efit1['intens']
	i2 = efit2['intens']
	result = QuickMatch(r1, i1, r2, i2, range=range, offset=offset, plot=plot,
						xlog=xlog, ylog=ylog)
	return result


	
def QuickMatch( r1, i1, r2, i2, range=None, offset=False, plot=True, xlog=False, ylog=False ):
	"""Function to permit quick, iterative matching of two profiles.
	The first profile, i1(r1), is the "reference" profile, and the goal is to
	scale the second profile, i2(r2), to match the reference.
	By default, we let the function determine the overlap region (maximized)
	and calculate the mean scaling necessary to match i2(r2) with i1(r1).
	The keyword parameter range can be used to specify a smaller matching
	range -- e.g. range=[45.0, 66.0].
	
	If offset = True, then an additive offset is calculated as well via
	least-squares optimization.
	
	If plot = True [the default], then i1(r1) is plotted, along with the
	scaled i2(r2).  The user can then zooom in, select a better overlap
	range, and re-run this function.  The plot is by default linear, but
	can be log in x or y or both by setting xlog and/or ylog to True.
	
	The output is the scale factor (and the offset, if requested), such that:
		i2_scaled = result * i2   {result = scalar}
		i2_scaled = result[0] * (i2 - result[1])   {result = 2-element list}
	"""
	
	if range is None:
		r1min = r1[0]

		r2min = r2[0]
		rmin = min([r1min, r2min])
		r1max = r1[0]
		r2max = r2[0]
		rmax = max([r1max, r2max])
	else:
		rmin = range[0]
		rmax = range[1]
	
	result = MatchProfiles(r1, i1, r2, i2, rmin, rmax, offset=offset)

	if plot is True:
		p.clf()
		if xlog is True:
			if ylog is True:
				p.loglog(r1, i1, 'k')
			else:
				p.semilogx(r1, i1, 'k')
		else:
			if ylog is True:
				p.semilogy(r1, i1, 'k')
			else:
				p.plot(r1, i1, 'k')
		i2_nd = N.array(i2)
		if offset is True:
			i2_corrected = result[0] * (i2_nd - result[1])
		else:
			i2_corrected = result * i2_nd
		p.plot(r2, i2_corrected, 'r')
	
	return result


