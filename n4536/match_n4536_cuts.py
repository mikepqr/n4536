import matchprofiles, profiles, os

d = os.path.expanduser("~/research/2012/smbhpb/data/n4536/cuts/")

r_nicmos, c_nicmos = profiles.ReadAndFoldProfile(d + "n4536_nicmos_pa120.0_w1.dat", pix=0.075)
r_irac1, c_irac1 = profiles.ReadAndFoldProfile(d + "n4536_irac1_pa120.0_w1.dat", pix=0.75)
r_2massk, c_2massk = profiles.ReadAndFoldProfile(d + "n4536_2massk_pa120.0_w1.dat", pix=1.)

irac1 = matchprofiles.QuickMatch(r_2massk, c_2massk, r_irac1, c_irac1, range=[8,30], offset=False)
nicmos = matchprofiles.QuickMatch(r_irac1, c_irac1, r_nicmos, c_nicmos, range=[3,6], offset=False)

# Radial cuts lack S/N and do not extend as far in radius.
