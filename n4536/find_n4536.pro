pro find_n4536, nblob = nblob
    fits_read, "~/research/2012/smbhpb/data/n4536/irac1/ngc4536_v7.phot.1.fit", img, h
    img[where(~finite(img))] = 0.
    find_galaxy, img, majoraxis, eps, axisangle, xpeak, ypeak, xmid, ymid, nblob = nblob
end
