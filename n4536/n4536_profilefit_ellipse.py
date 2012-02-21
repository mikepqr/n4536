import os, subprocess
import numpy as n
import pylab as p
import profiles, imfit, plotutils
import match_n4536_elfixed0epa

d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/profilefit_ellipse/")
exe = os.path.expanduser("~/research/2012/smbhpb/code/imfit/profilefit")
datfile = d + "n4536_nicmosirac1_el_fixed0epa.dat"
cfgfile = d + "n4536_nicmosirac1_el_fixed0epa.config"
outfile = d + "n4536_nicmosirac1_el_fixed0epa.profilefit"
opt1 = "--usemask"
opt2 = "--save-params"

subprocess.call([exe, datfile, cfgfile, opt1, opt2, outfile])

r_dat, mu_dat = profiles.ReadProfile(datfile)

funcs, params = imfit.ReadConfigFile(outfile)

plotutils.PlotFit(r_dat, mu_dat, funcs, params, xlog=True, yrange=[25,9])
p.axvline(match_n4536_elfixed0epa.inner, color = "Gray")
p.axvline(match_n4536_elfixed0epa.outer, color = "Gray")
p.ylim(-1,1)
p.axes([0.1, 0.35, 0.8, 0.55])
p.title("profilefit of concatenated NICMOS/IRAC1 fixed ellipse fit profile")
p.show()
