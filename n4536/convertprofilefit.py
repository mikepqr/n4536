import imfit
import os, copy

pix = 0.75
zp = 18.1639

def convertprofilefit(d = "~/research/2012/smbhpb/data/n4536/profilefit/", \
        root = "n4536_nicmosirac1_pa120.0_w1"):

    d = os.path.expanduser(d)
    infile = d + root + ".profilefit"

    funcs, params = imfit.ReadConfigFile(infile)

    params_converted = copy.deepcopy(params)
    # Sersic surface brightness to counts
    params_converted[0][2] = 10**(0.4 * (zp - params[0][2]))
    # Exponential surface brightness to counts
    params_converted[1][1] = 10**(0.4 * (zp - params[1][1]))
    # Sersic scale to arcseconds
    params_converted[0][3] = params[0][3]/pix
    # Exponential surface brightness to counts
    params_converted[1][2] = params[1][2]/pix

    print params
    print params_converted
