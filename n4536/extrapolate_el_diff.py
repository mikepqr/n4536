import ellipse
import pylab as p
import numpy as n
import os

d = "/Users/mike/research/2012/smbhpb/data/n4536/profilefit/"
el_minus_bulge_file = "el_n4536_irac_minus_bulge.ascii"
el_minus_disk_file = "el_n4536_irac_minus_disk.ascii"
el_minus_bulge_benderfile = os.path.splitext(el_minus_bulge_file)[0] + ".bender"
el_minus_disk_benderfile = os.path.splitext(el_minus_disk_file)[0] + ".bender"

iracpix = 0.75
iraczp = 18.1639

el_minus_bulge = ellipse.ReadEllipse(d + el_minus_bulge_file, pix = iracpix, ZP = iraczp)
el_minus_disk = ellipse.ReadEllipse(d + el_minus_disk_file, pix = iracpix, ZP = iraczp)
p.loglog(el_minus_bulge.a, el_minus_bulge.i, color = "LightSteelBlue")
p.loglog(el_minus_disk.a, el_minus_disk.i, color = "SandyBrown")

# Extrapolate inwards assuming constant sb
imax = n.argmax(el_minus_bulge.i)
el_minus_bulge.i[0:imax] = el_minus_bulge.i[imax]
el_minus_bulge.sb[0:imax] = el_minus_bulge.sb[imax]

# Get rid of negative counts
neg = el_minus_bulge.i < 0.
el_minus_bulge.i[neg] = 1e-4
el_minus_bulge.sb[neg] = -2.5 * n.log10(el_minus_bulge.i[neg]) + iraczp

p.loglog(el_minus_bulge.a[0:imax], el_minus_bulge.i[0:imax], "o", color = "LightSteelBlue")
p.show()

el_minus_bulge_bender = ellipse.ConvertIrafToBender(el_minus_bulge, dataFrame =  False)
ellipse.WriteBenderEllipse(el_minus_bulge_bender, d + el_minus_bulge_benderfile)
