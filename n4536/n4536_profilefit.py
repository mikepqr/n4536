import os, subprocess
import numpy as n
import pylab as p
import profiles, imfit, plotutils
import concatenate_nicmos_sings

d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/profilefit_tinybulge/")
d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/profilefit/")
exe = os.path.expanduser("~/research/2012/smbhpb/code/imfit/profilefit")
datfile = d + "n4536_nicmosirac1_pa120.0_w1.dat"
cfgfile = d + "n4536_nicmosirac1_pa120.0_w1.config"
outfile = d + "n4536_nicmosirac1_pa120.0_w1.profilefit"
opt1 = "--usemask"
opt2 = "--save-params"

subprocess.call([exe, datfile, cfgfile, opt1, opt2, outfile])

r_dat, mu_dat = profiles.ReadProfile(datfile)

funcs, params = imfit.ReadConfigFile(outfile)

plotutils.PlotFit(r_dat, mu_dat, funcs, params, xlog=True, yrange=[25,9])
p.axvline(concatenate_nicmos_sings.inner, color = "Gray")
p.axvline(concatenate_nicmos_sings.outer, color = "Gray")
p.ylim(-1,1)
p.axes([0.1, 0.35, 0.8, 0.55])
p.title("profilefit of concatenated NICMOS/IRAC1 major axis profile")
p.show()
