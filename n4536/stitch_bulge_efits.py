import ellipse
import pylab as p

iracpix = 0.75
iraczp = 18.1639
nicmospix = 0.075
nicmoszp = 15.1304
nicmostelPA=246.74

el_nicmos_minus_disk = ellipse.ReadEllipse("/Users/mike/research/2012/smbhpb/data/n4536/profilefit/el_n4536_nicmos_minus_disk.ascii", pix = nicmospix, telPA = nicmostelPA, ZP = nicmoszp)

el_irac_minus_disk = ellipse.ReadEllipse("/Users/mike/research/2012/smbhpb/data/n4536/profilefit/el_n4536_irac_minus_disk.ascii", pix = iracpix, ZP = iraczp)

p.figure(0)
p.clf()
(ax1, ax2) = ellipse.PlotEllPA([el_irac_minus_disk, el_nicmos_minus_disk], xlog=True)
p.suptitle("Free e/PA ellipse fit to IRAC1 bulge image (black), NICMOS bulge (red)")
p.show()

p.figure(1)
p.clf()
a, = p.semilogx(el_irac_minus_disk.a, el_irac_minus_disk.sb)
b, = p.semilogx(el_nicmos_minus_disk.a, el_nicmos_minus_disk.sb)
p.legend([a, b], ["IRAC", "NICMOS"])
p.ylim(tuple(reversed(p.ylim())))
p.ylabel(r"$\mu_K$ (mag arcsec$^{-2}$)")
p.xlabel(r"$r$ (arcsec)")
p.suptitle("Free e/PA ellipse fit to IRAC1 bulge image (black), NICMOS bulge (red)")
p.show()

