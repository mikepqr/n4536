import os, subprocess
import numpy as n
import pylab as p
import pyfits
import profiles

def make_n4536_image_irac(d = "~/research/2012/smbhpb/data/n4536/profilefit/", \
        root = "n4536_nicmosirac1_pa120.0_w1"):

    d = os.path.expanduser(d)
    exe = os.path.expanduser("~/research/2012/smbhpb/code/imfit/makeimage")
    datfile = d + root + ".dat"
    cfgfile = d + root + ".profilefit_irac"
    refimage = d + "../irac1/ngc4536_v7.phot.1.fit"
    output_root = d + "n4536_irac_conceptual"
    psfimg = d + "psf_moffat_irac.fits"

    # Make image
    cmd = " ".join([exe, cfgfile, "--output", output_root + ".fits", \
            "--refimage", refimage, "--output-functions", output_root, \
            "--psf", psfimg])
    subprocess.call(cmd, shell=True)

    # Subtract bulge from IRAC image and save
    irac = pyfits.open(refimage)
    bulge = pyfits.open(output_root + "1_Sersic.fits")
    irac[0].data -= bulge[0].data
    irac.writeto(d + "n4536_irac_minus_bulge.fits", clobber = True)

    # Subtract disk from IRAC image and save
    irac = pyfits.open(refimage)
    disk = pyfits.open(output_root + "2_Exponential.fits")
    irac[0].data -= disk[0].data
    irac.writeto(d + "n4536_irac_minus_disk.fits", clobber = True)

def make_n4536_image_nicmos(d = "~/research/2012/smbhpb/data/n4536/profilefit/", \
        root = "n4536_nicmosirac1_pa120.0_w1"):

    d = os.path.expanduser(d)
    exe = os.path.expanduser("~/research/2012/smbhpb/code/imfit/makeimage")
    datfile = d + root + ".dat"
    cfgfile = d + root + ".profilefit_nicmos"
    refimage = d + "../nicmos/n4536_nicmos.fits"
    output_root = d + "n4536_nicmos_conceptual"
    psfimg = d + "psf_tinytim_nicmos.fits"

    # Make image
    cmd = " ".join([exe, cfgfile, "--output", output_root + ".fits", \
            "--refimage", refimage, "--output-functions", output_root, \
            "--psf", psfimg])
    subprocess.call(cmd, shell=True)

    # Subtract bulge from NICMOS image and save
    nicmos = pyfits.open(refimage)
    bulge = pyfits.open(output_root + "1_Sersic.fits")
    nicmos[0].data = nicmos[0].data - bulge[0].data
    nicmos.writeto(d + "n4536_nicmos_minus_bulge.fits", clobber = True)

    # Subtract disk from nicmos image and save
    nicmos = pyfits.open(refimage)
    disk = pyfits.open(output_root + "2_Exponential.fits")
    nicmos[0].data -= disk[0].data
    nicmos.writeto(d + "n4536_nicmos_minus_disk.fit", clobber = True)
