import profiles
import os, subprocess
import pylab as p
import numpy as n

l =  1000
xc = 992.18
yc = 967.53
pa1 = 110.0
pa2 = 125.0
pix_irac = 0.75
d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/")
exe = os.path.expanduser("~/research/2012/smbhpb/code/getgalaxyprofiles.py")
output_dir = d + "profilefit/"
# getgalaxyprofiles doesn't work unless input/output dir is working dir
os.chdir(output_dir)

output_root = "n4536_conceptual"
suffix = "1_Sersic"

# # Write profiles at different PAs
# cmd = " ".join([exe, output_root + suffix + ".fits", str(l), \
#         "--output=" + output_root + suffix, "--width=1", \
#         "--skyPA=" + str(pa1), "--xc="+str(xc), "--yc="+str(yc)])
# print cmd
# subprocess.call(cmd, shell=True)
# cmd = " ".join([exe, output_root + suffix + ".fits", str(l), \
#         "--output=" + output_root + suffix, "--width=1", \
#         "--skyPA=" + str(pa2), "--xc="+str(xc), "--yc="+str(yc)])
# print cmd
# subprocess.call(cmd, shell=True)

pvect_suffix = "_pa" + str(pa1) + "_w1.dat"
r_b1, c_b1 = profiles.ReadAndFoldProfile(output_dir + \
        output_root + suffix + pvect_suffix, pix =  pix_irac)
pvect_suffix = "_pa" + str(pa2) + "_w1.dat"
r_b2, c_b2 = profiles.ReadAndFoldProfile(output_dir + \
        output_root + suffix + pvect_suffix, pix =  pix_irac)

p.figure(0)
p.clf()
p.suptitle("Cuts through conceputal model at bulge PA (%g) and disk PA (%g)" % (pa1, pa2))
a, = p.loglog(r_b1, c_b1, color = "Red")
b, = p.loglog(r_b2, c_b2, color = "Blue")
p.ylabel("IRAC counts")
p.xlabel("Radius [arcsec]")
p.legend([a,b], [str(pa1), str(pa2)])
p.show()
