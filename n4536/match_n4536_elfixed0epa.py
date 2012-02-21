import matchprofiles, ellipse, math, os, profiles
import pylab as p
import numpy as n

merge_radius = 5.
inner = 0.5
outer = 90.
tmasskzp = 19.868
irac_range = [8, 30]
nicmos_range = [4, 15]

def match_n4536_elfixed0epa():
    d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/el_fixed0epa/")

    el_nicmos = ellipse.ReadEllipse(d + "el_n4536_nicmos_fixed0epa.ascii", pix=0.075)
    el_nicmos_lin = ellipse.ReadEllipse(d + "el_n4536_nicmos_fixed0epa_lin0.5.ascii", pix=0.075)
    el_2massk = ellipse.ReadEllipse(d + "el_n4536_2massk_fixed0epa.ascii", pix=1.0)
    el_irac1 = ellipse.ReadEllipse(d + "el_n4536_irac1_fixed0epa.ascii", pix=0.75)
    el_irac1_lin = ellipse.ReadEllipse(d + "el_n4536_irac1_fixed0epa_lin0.5.ascii", pix=0.75)

    print "\n2MASS K sb in 2MASS K: mu [mag/arcsec2] = -2.5 log C + %g\n" % tmasskzp

    # Match IRAC1 to 2MASS
    irac1 = matchprofiles.QuickMatchEfits(el_2massk, el_irac1, range = irac_range, offset = False, plot = False)
    iraczp = -2.5 * math.log10(irac1) + 19.868
    p.figure(0)
    p.clf()
    a, = p.loglog(el_2massk.a, el_2massk.i)
    b, = p.loglog(el_irac1.a, irac1 * el_irac1.i)
    p.axvline(irac_range[0])
    p.axvline(irac_range[1])
    p.title("Match IRAC1 to 2MASS K fixed ellipse fit profile")
    p.legend([a, b], ["2MASS K", "%g * IRAC1" % irac1])
    p.show()

    print "\nIRAC1 sb in 2MASS K: mu [mag/arcsec2] = -2.5 log (%g * C) + 19.868" % irac1
    print "                                      = -2.5 log C + %g" % iraczp
    print "               IRAC1 C [counts/pixel] = 10**[0.4 * (%g - mu)]\n" % iraczp

    # Match NICMOS to IRAC1
    nicmos = matchprofiles.QuickMatchEfits(el_irac1_lin, el_nicmos_lin, range = nicmos_range, offset = True, plot = False)
    nicmoszp = -2.5 * math.log10(irac1 * nicmos[0]) + 19.868
    p.figure(1)
    p.clf()
    a, = p.loglog(el_irac1.a, el_irac1.i)
    b, = p.loglog(el_nicmos.a, nicmos[0] * (el_nicmos.i - nicmos[1]))
    p.axvline(nicmos_range[0])
    p.axvline(nicmos_range[1])
    p.title("Match NICMOS to IRAC1 fixed ellipse fit profile")
    p.legend([a, b], ["IRAC1", "%g * (NICMOS - %g)" % tuple(nicmos)])
    p.show()

    print "\nNICMOS sb in 2MASS K: mu [mag/arcsec2] = -2.5 log (%g * (%g * (C - %g))) + 19.868" % (irac1, nicmos[0], nicmos[1])
    print "                                       = -2.5 log C' + %g" % nicmoszp
    print "                                         (where C' = C + %g)" % (-1. * nicmos[1])
    print "                             NICMOS C' = 10**[0.4 (%g - mu)]\n" % nicmoszp

    # Concatenate and save matched NICMOS and IRAC1 profiles

    mu_nicmos = -2.5 * n.log10(el_nicmos.i - nicmos[1]) + nicmoszp
    mu_irac1 = -2.5 * n.log10(el_irac1.i) + iraczp
    r, mu = profiles.MergeTwoProfiles(el_nicmos.a, mu_nicmos, el_irac1.a, mu_irac1, merge_radius)
    unmasked = (r > inner) & (r < outer) & n.isfinite(mu)
    masked = ~unmasked
    maskIndices = masked.nonzero()[0]

    profiles.WriteProfile(r, mu, d + "n4536_nicmosirac1_el_fixed0epa.dat", mask=True, \
            maskIndices=maskIndices)
    p.figure(2)
    p.clf()
    a, = p.semilogx(el_nicmos.a, mu_nicmos)
    b, = p.semilogx(el_irac1.a, mu_irac1)
    c, = p.semilogx(r, mu - 1., color = "Gray")
    d, = p.semilogx(r[unmasked], mu[unmasked], ".", color = "Red")

    p.axvline(merge_radius, color = "Gray")
    p.axvline(inner, color = "Gray")
    p.axvline(outer, color = "Gray")

    p.ylabel(r"$\mu_K$ [mag arcsec$^{-2}$]")
    p.title("Fixed ellipse fits to IRAC and NICMOS")
    p.legend([a, b, c], ["NICMOS", "IRAC1", "Merged + arbitrary offset"])
    p.ylim(tuple(reversed(p.ylim())))
    p.show()
