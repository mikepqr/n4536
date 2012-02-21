# rescale then concatenate nicmos and irac1 major axis cuts in units of 2MASS K mag/arcsec2

import profiles, os
import pylab as p
import numpy as n

nicmoszp = 15.1304
irac1zp = 18.1639
merge_radius = 5.
inner = 0.5
outer = 90.

def concatenate_nicmos_sings():
    d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/cuts/")
    h = """# Major axis (PA = 120) cut of n4536. SBs in 2MASS K mag/arcsec2. Formed by
# concatenating NICMOS and IRAC1 cuts at r = %g arcsec. Masked r < %g, r > %g""" % (merge_radius, inner, outer)

    r_nicmos, c_nicmos = profiles.ReadAndFoldProfile(d + "n4536_nicmos_pa120.0_w1.dat", pix = 0.075)
    r_irac1, c_irac1 = profiles.ReadAndFoldProfile(d + "n4536_irac1_pa120.0_w1.dat", pix = 0.75)

    mu_nicmos = -2.5 * n.log10(c_nicmos) + nicmoszp
    mu_irac1 = -2.5 * n.log10(c_irac1) + irac1zp
    r, mu = profiles.MergeTwoProfiles(r_nicmos, mu_nicmos, r_irac1, mu_irac1, merge_radius)

    # unmasked is a boolean array. ~ unary operator inverts it somehow. nonzero 
    # returns tuple containing list of indices where True.
    unmasked = (r > inner) & (r < outer) & n.isfinite(mu)
    masked = ~unmasked
    maskIndices = masked.nonzero()[0]

    profiles.WriteProfile(r, mu, d + "n4536_nicmosirac1_pa120.0_w1.dat", mask=True, \
            header=h, maskIndices=maskIndices)

    p.figure(0)
    p.clf()
    a, = p.semilogx(r_nicmos, mu_nicmos)
    b, = p.semilogx(r_irac1, mu_irac1)
    c, = p.semilogx(r, mu - 1., color = "Gray")
    d, = p.semilogx(r[unmasked], mu[unmasked] - 1., ".", color = "Red")

    p.axvline(merge_radius, color = "Gray")
    p.axvline(inner, color = "Gray")
    p.axvline(outer, color = "Gray")

    p.ylabel(r"$\mu_K$ [mag arcsec$^{-2}$]")
    p.title("Major axis cuts of NICMOS and IRAC1 images")
    p.legend([a, b, c], ["NICMOS", "IRAC1", "Merged + arbitrary offset"])
    p.ylim(tuple(reversed(p.ylim())))
    p.show()

if __name__ == "__main__":
    concatenate_nicmos_sings()
