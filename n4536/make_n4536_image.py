import os, subprocess
import numpy as n
import pylab as p
import pyfits
import profiles

def make_n4536_image(d = "~/research/2012/smbhpb/data/n4536/profilefit/", \
        root = "n4536_nicmosirac1_pa120.0_w1"):

    d = os.path.expanduser(d)
    exe = os.path.expanduser("~/research/2012/smbhpb/code/imfit/makeimage")
    datfile = d + root + ".dat"
    cfgfile = d + root + ".profilefitconverted"
    refimage = d + "../irac1/ngc4536_v7.phot.1.fit"
    output_root = d + "n4536_conceptual"
    psfimg = d + "psf_moffat_irac.fits"

    # Make image
    cmd = " ".join([exe, cfgfile, "--output", output_root + ".fits", \
            "--refimage", refimage, "--output-functions", output_root, \
            "--psf", psfimg])
    subprocess.call(cmd, shell=True)

    # Subtract bulge from IRAC image and save
    irac = pyfits.open(refimage)
    conceptual = pyfits.open(output_root + ".fits")
    bulge = pyfits.open(output_root + "1_Sersic.fits")
    disk = pyfits.open(output_root + "2_Exponential.fits")
    irac[0].data = irac[0].data - bulge[0].data
    irac.writeto(d + "n4536_irac_minus_bulge.fits", clobber = True)
