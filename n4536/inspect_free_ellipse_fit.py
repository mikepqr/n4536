import ellipse
import os
import pylab as p

d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/")

el_irac1_fixed = ellipse.ReadEllipse(d + "el_fixed0epa/el_n4536_irac1_fixed0epa.ascii", pix=0.75)
el_irac1_free = ellipse.ReadEllipse(d + "el_free/el_n4536_irac1_free.ascii", pix=0.75)

p.figure(0)
p.clf()
(ax1, ax2) = ellipse.PlotEllPA(el_irac1_free)
ax1.axhline(110., color = "Gray")
ax1.axhline(120., color = "Gray")
ax1.axhline(130., color = "Gray")
ax2.axhline(0.5, color = "Gray")
ax2.axhline(0.6, color = "Gray")
p.suptitle("Free e/PA ellipse fit to IRAC1 image")
p.show()

p.figure(1)
p.clf()
a, = p.loglog(el_irac1_free.a, el_irac1_free.i)
b, = p.loglog(el_irac1_fixed.a, el_irac1_fixed.i)
p.legend([a, b], ["Free", "Fixed"])
p.suptitle("Free e/PA and fixed e/PA ellipse fits to IRAC1 image\nEllipse fit center fixed for both")
p.show()
