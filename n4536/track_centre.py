import ellipse
import pylab as p
import os

elfile = os.path.expanduser("~/research/2012/smbhpb/data/n4536/el_free/el_n4536_irac1_free.ascii")

n4536_irac1_free = ellipse.ReadEllipse(elfile, pix=0.75)

p.figure(0)
p.clf()
p.plot(n4536_irac1_free.x0, n4536_irac1_free.y0, "o")
p.show()
