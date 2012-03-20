import ellipse
import numpy as n
import pylab as p
import copy

d = "/Users/mike/research/2012/smbhpb/data/n4536/profilefit/"
el_disk_file = "el_n4536_irac_minus_bulge.ascii"
el_disk_benderfile = "el_n4536_disk.bender"
el_bulge_benderfile = "el_n4536_bulge.bender"
iracpix = 0.75
iraczp = 18.1639

el_disk = ellipse.ReadEllipse(d + el_disk_file, pix = iracpix, ZP = iraczp, dataFrame = False)
# Get rid of negative counts
neg = el_disk['intens'] < 0.
el_disk['intens'][neg] = 1e-4
el_disk['sb'][neg] = -2.5 * n.log10(el_disk['intens'][neg]) + iraczp

# Use bulge as dummy efit objects with correct radii
# Replace bulge with dummy SB and ell values
el_inner_disk = ellipse.ReadBenderEllipse(d + el_bulge_benderfile)
nPts = len(el_inner_disk['a'])
el_inner_disk['eps'][:] = el_disk['ellip'][0]
el_inner_disk['pa'] = n.zeros(nPts)
el_inner_disk['sb'] = -2.5 * n.log10(el_disk['intens'].max()) + iraczp
# TODO: deal with deps, dpa, a4

el_inner_disk = ellipse.ConvertBenderToIraf(el_inner_disk, dataFrame = False)

# Get radius at which maximum SB occurs in disk efit
rmax = el_disk['sma'][n.argmax(el_disk['intens'])]

# TODO: Merge el_disk and el_inner_disk at rmax
# Inspect around rmax to check it's not too jumpy
# Read IN bender ellipse files and plot to check no basic errors

el_disk_bender = ellipse.ConvertIrafToBender(el_disk, dataFrame =  False, irafColnames = [])
ellipse.WriteBenderEllipse(el_disk_bender, d + el_disk_benderfile)
