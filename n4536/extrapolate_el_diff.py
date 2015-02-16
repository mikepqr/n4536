import ellipse
import numpy as n
import pylab as p
import copy

d = "/Users/mike/research/2012/smbhpb/data/n4536/profilefit/"
el_disk_file = "el_n4536_irac_minus_bulge.ascii"
el_disk_benderfile = "el_n4536_disk.bender"
el_inner_disk_file = "el_n4536_nicmos_minus_disk.ascii"
iracpix = 0.75
iraczp = 18.1639
nicmospix = 0.075
nicmoszp = 15.1304

el_disk = ellipse.ReadEllipse(d + el_disk_file, pix = iracpix, ZP = iraczp, dataFrame = False)
# Get rid of negative counts
neg = el_disk['intens'] < 0.
el_disk['intens'][neg] = 1e-4
el_disk['sb'][neg] = -2.5 * n.log10(el_disk['intens'][neg]) + iraczp
# Manually set a4 to zero
el_disk['a4'][:] = 0.

# Use NICMOS efit as dummy efit with correct radii
# Replace this with dummy SB and ell values
el_inner_disk = ellipse.ReadEllipse(d + el_inner_disk_file, dataFrame = False, pix = nicmospix, ZP = nicmoszp)
el_inner_disk['ellip'][:] = el_disk['ellip'][0]
el_inner_disk['pa'][:] = el_disk['pa'][0]
el_inner_disk['sb'][:] = -2.5 * n.log10(el_disk['intens'].max()) + iraczp

# Get radius at which maximum SB occurs in disk efit
rmax = el_disk['sma'][n.argmax(el_disk['intens'])]
# Merge el_disk and el_inner_disk at rmax
el_n4536_disk = ellipse.MergeEllipseFits(el_inner_disk, el_disk, rmax)

el_n4536_disk_bender = ellipse.ConvertIrafToBender(el_n4536_disk, dataFrame =  False, irafColnames = [])
# Initialize deps, a4, PA, dpa
el_n4536_disk_bender['dpa'][:] = 0.
el_n4536_disk_bender['deps'][:] = 0.
el_n4536_disk_bender['a4'][:] = 0.
ellipse.WriteBenderEllipse(el_n4536_disk_bender, d + el_disk_benderfile)

# Plot
p.clf()
a, = p.semilogx(el_disk['sma'], el_disk['sb'])
b, = p.semilogx(el_inner_disk['sma'], el_inner_disk['sb'])
c, = p.semilogx(el_n4536_disk['sma'], el_n4536_disk['sb'])
p.ylim(tuple(reversed(p.ylim())))
p.show()
