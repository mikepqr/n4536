import imfit
import os, copy

iracpix = 0.75
iraczp = 18.1639
nicmospix = 0.075
nicmoszp = 15.1304

def convertprofilefit(d = "~/research/2012/smbhpb/data/n4536/profilefit/", \
        root = "n4536_nicmosirac1_pa120.0_w1"):

    d = os.path.expanduser(d)
    infile = d + root + ".profilefit"

    funcs, params = imfit.ReadConfigFile(infile)

    params_irac = copy.deepcopy(params)
    # Sersic surface brightness to counts
    params_irac[0][2] = 10**(0.4 * (iraczp - params[0][2]))
    # Exponential surface brightness to counts
    params_irac[1][1] = 10**(0.4 * (iraczp - params[1][1]))
    # Sersic scale to arcseconds
    params_irac[0][3] = params[0][3]/iracpix
    # Exponential surface brightness to counts
    params_irac[1][2] = params[1][2]/iracpix

    params_nicmos = copy.deepcopy(params)
    # Sersic surface brightness to counts
    params_nicmos[0][2] = 10**(0.4 * (nicmoszp - params[0][2]))
    # Exponential surface brightness to counts
    params_nicmos[1][1] = 10**(0.4 * (nicmoszp - params[1][1]))
    # Sersic scale to arcseconds
    params_nicmos[0][3] = params[0][3]/nicmospix
    # Exponential surface brightness to counts
    params_nicmos[1][2] = params[1][2]/nicmospix

    print params
    print params_irac
    print params_nicmos
