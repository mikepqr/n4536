import pylab as p
import ellipse, profiles, plotutils

d = "/Users/mike/research/2012/smbhpb/data/n4536/profilefit/"
el_bulge_bender = ellipse.ReadBenderEllipse(d + "el_n4536_bulge.bender")
el_disk_bender = ellipse.ReadBenderEllipse(d + "el_n4536_disk.bender")
el_bulge = ellipse.ConvertBenderToIraf(el_bulge_bender)
el_disk = ellipse.ConvertBenderToIraf(el_disk_bender)
r_concat, mu_concat = profiles.ReadProfile(d + "n4536_nicmosirac1_pa120.0_w1.dat")

p.figure(0)
p.clf()
ellipse.PlotEllPA(el_bulge, xlog=True)
p.show()

p.figure(1)
p.clf()
b, = p.semilogx(el_disk['sma'], el_disk['sb'])
a, = p.semilogx(el_bulge['sma'], el_bulge['sb'])
c, = p.semilogx(r_concat, mu_concat)
d, = p.semilogx(el_disk['sma'], plotutils.addmag(el_disk['sb'], el_bulge['sb']))
p.ylim(tuple(reversed(p.ylim())))
p.show()
