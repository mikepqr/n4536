import os, subprocess
import numpy as n
import pylab as p
import profiles

zp_irac = 18.1639
irac_counts_per_nicmos_count = 16.3456
pix_irac = 0.75
pix_nicmos = 0.075
l =  1000
xc = 992.18
yc = 967.53
pa = 120.0

def compare_model(d = "~/research/2012/smbhpb/data/n4536/profilefit/"):
    # Take major axis cuts of total model, bulge and disk

    d = os.path.expanduser(d)
    exe = os.path.expanduser("~/research/2012/smbhpb/code/getgalaxyprofiles.py")
    # getgalaxyprofiles doesn't work unless input/output dir is working dir
    os.chdir(d)
    output_root = "n4536_conceptual"
    suffix = ["", "1_Sersic", "2_Exponential"]
    pvect_suffix = "_pa" + str(pa) + "_w1.dat"
    for i in [0, 1, 2]:
        cmd = " ".join([exe, output_root + suffix[i] + ".fits", str(l), \
                "--output=" + output_root + suffix[i], "--width=1", \
                "--skyPA=" + str(pa), "--xc="+str(xc), "--yc="+str(yc)])
        print cmd
        subprocess.call(cmd, shell=True)

    # Read major axis cuts of total model, bulge and disk
    r_conceptual, c_conceptual = profiles.ReadAndFoldProfile(d + \
            output_root + suffix[0] + pvect_suffix, pix =  pix_irac)
    r_bulge, c_bulge = profiles.ReadAndFoldProfile(d + \
            output_root + suffix[1] + pvect_suffix, pix =  pix_irac)
    r_disk, c_disk = profiles.ReadAndFoldProfile(d + \
            output_root + suffix[2] + pvect_suffix, pix =  pix_irac)

    # Read in major axis profiles and rescale to IRAC counts where necessary of images
    r_concat, mu_concat = profiles.ReadProfile(d + "../cuts/n4536_nicmosirac1" + pvect_suffix, pix=1)
    c_concat = 10**(0.4 * (zp_irac - mu_concat))
    r_nicmos, c_nicmos = profiles.ReadAndFoldProfile(d + "../cuts/n4536_nicmos" + pvect_suffix, pix=pix_nicmos)
    c_nicmos *= irac_counts_per_nicmos_count
    r_irac, c_irac = profiles.ReadAndFoldProfile(d + "../cuts/n4536_irac1" + pvect_suffix, pix=pix_irac)


    cmd = " ".join([exe, "n4536_irac_minus_bulge.fits", str(l), \
            "--output=" + "n4536_irac_minus_bulge", "--width=1", \
            "--skyPA=" + str(pa), "--xc="+str(xc), "--yc="+str(yc)])
    print cmd
    subprocess.call(cmd, shell=True)
    r_irac_minus_bulge, c_irac_minus_bulge = profiles.ReadAndFoldProfile(\
            d + "n4536_irac_minus_bulge" + pvect_suffix, pix = pix_irac)

    cmd = " ".join([exe, "n4536_irac_minus_disk.fits", str(l), \
            "--output=" + "n4536_irac_minus_disk", "--width=1", \
            "--skyPA=" + str(pa), "--xc="+str(xc), "--yc="+str(yc)])
    print cmd
    subprocess.call(cmd, shell=True)
    r_irac_minus_disk, c_irac_minus_disk = profiles.ReadAndFoldProfile(\
            d + "n4536_irac_minus_disk" + pvect_suffix, pix = pix_irac)

    p.figure(0)
    p.clf()
    p.suptitle("Cuts through image and conceputal model at sky PA = %g degrees" % pa)
    a, = p.loglog(r_irac, c_irac, "o", mfc='None', color="Gray")
    b, = p.loglog(r_nicmos, c_nicmos, "o", color="Gray")
    c, = p.loglog(r_disk, c_disk, color = "Blue")
    d, = p.loglog(r_bulge, c_bulge, color = "Red")
    e, = p.loglog(r_conceptual, c_conceptual, color = "Black")
    f, = p.loglog(r_irac_minus_bulge, c_irac_minus_bulge, "o", color = "Blue")
    g, = p.loglog(r_irac_minus_disk, c_irac_minus_disk, "o", color = "Red")
    p.ylabel("IRAC counts")
    p.xlabel("Radius [arcsec]")
    p.legend([a,b,c,d,e,f,g], ["IRAC", "NICMOS", "Disk model", "Bulge model", "Total model", "IRAC - bulge model", "IRAC - disk model"])
    p.show()
