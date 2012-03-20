import ellipse
import copy
import pylab as p
import numpy as n
from scipy import optimize

# Constants
iracpix = 0.75
iraczp = 18.1639
nicmospix = 0.075
nicmoszp = 15.1304
nicmos_counts_per_irac_count = 16.3456
nicmostelPA=246.74
ell_transition_radius = 6.
flux_transition_radius = 15.
fitrange = [2., 18.]
replacement_ell = 0.5
d = "/Users/mike/research/2012/smbhpb/data/n4536/profilefit/"
el_bulge_bender_file = "el_n4536_bulge.bender"

# Read ellipse fits
el_nicmos_minus_disk = ellipse.ReadEllipse(d + "el_n4536_nicmos_minus_disk.ascii", pix = nicmospix, telPA = nicmostelPA, ZP = nicmoszp, dataFrame = False)
el_irac_minus_disk = ellipse.ReadEllipse(d + "el_n4536_irac_minus_disk.ascii", pix = iracpix, ZP = iraczp, dataFrame = False)
# Unused SINFONI efit
el_sinfoni_100comb = ellipse.ReadEllipse("/Users/mike/research/2012/smbhpb/data/n4536/sinfoni_collapsed/el_sinfoni_100comb.ascii", pix=0.05, telPA=125.0, dataFrame = False)

# Replace ellipticity with constant beyond transition radius
el_irac_minus_disk_replaced = copy.deepcopy(el_irac_minus_disk)
ii = el_irac_minus_disk['sma'] > ell_transition_radius
el_irac_minus_disk_replaced['ellip'][ii] = replacement_ell

# Replace IRAC flux with extrapolation beyond transition radius #2
# Assemble coordinates for fitting
ii = (el_irac_minus_disk['sma'] > ell_transition_radius) & (el_irac_minus_disk['sma'] < fitrange[1])
xx_irac = el_irac_minus_disk['sma'][ii]
yy_irac = el_irac_minus_disk['intens'][ii]
ii = (el_nicmos_minus_disk['sma'] > fitrange[0]) & (el_nicmos_minus_disk['sma'] < ell_transition_radius)
xx = n.append(xx_irac, el_nicmos_minus_disk['sma'][ii])
yy = n.append(yy_irac, nicmos_counts_per_irac_count * el_nicmos_minus_disk['intens'][ii])
# Fit
sersic = lambda pt, x: pt[0] * n.exp(-1. * ((x / pt[1])**(1./pt[2])))
# / y for traditional chi^2. /y**2 to weight larger radii more
err = lambda pt, x, y: (sersic(pt, x) - y) / y
p0 = [50., 10., 4.]
p1, success = optimize.leastsq(err, p0, args=(xx, yy))
# Save fit for plotting
xt = el_irac_minus_disk['sma'][el_irac_minus_disk['sma'] > fitrange[0]]
yt = sersic(p1, xt)
el_irac_minus_disk_replaced['intens'][el_irac_minus_disk_replaced['sma'] > flux_transition_radius] = yt[xt > flux_transition_radius]
el_irac_minus_disk_replaced['sb'][:] = -2.5 * n.log10(el_irac_minus_disk_replaced['intens']) + iraczp

# Now merge NICMOS and tidied IRAC efit and write in Bender format
el_n4536_bulge = ellipse.MergeEllipseFits(el_nicmos_minus_disk, el_irac_minus_disk_replaced, ell_transition_radius)
el_n4536_bulge_bender = ellipse.ConvertIrafToBender(el_n4536_bulge, irafColnames = [], dataFrame = False)
ellipse.WriteBenderEllipse(el_n4536_bulge_bender, d + el_bulge_bender_file)

# Plot results
# Ell/PA
p.figure(0)
p.clf()
(ax1, ax2) = ellipse.PlotEllPA([el_irac_minus_disk, el_nicmos_minus_disk, el_n4536_bulge], xlog=True)
p.suptitle("Free e/PA ellipse fit to IRAC1 bulge image (black), NICMOS bulge (red)")
ax1.axvline(ell_transition_radius, color = "Gray")
ax2.axvline(ell_transition_radius, color = "Gray")
ax2.axhline(0.5, color = "Gray")
p.show()
# Surface brightness
p.figure(1)
p.clf()
a, = p.semilogx(el_irac_minus_disk['sma'], el_irac_minus_disk['sb'])
b, = p.semilogx(el_n4536_bulge['sma'], el_n4536_bulge['sb'] + 1)
c, = p.semilogx(el_nicmos_minus_disk['sma'], el_nicmos_minus_disk['sb'])
d, = p.semilogx(xt, -2.5 * n.log10(yt) + iraczp)
p.legend([a, b, c, d], ["IRAC", "Stitched + offset", "NICMOS", "Fit"])
p.ylabel(r"$\mu_K$ (mag arcsec$^{-2}$)")
p.xlabel(r"$r$ (arcsec)")
p.axvline(ell_transition_radius, color = "Gray")
p.axvline(flux_transition_radius, color = "Gray")
p.suptitle("Free e/PA ellipse fit to IRAC1 bulge image (black), NICMOS bulge (red)")
p.xlim(0.07,100.)
p.ylim(8, 30)
p.ylim(tuple(reversed(p.ylim())))
p.show()
