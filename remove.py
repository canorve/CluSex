#! /usr/bin/env python

import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path
import scipy


def main():


    if len(sys.argv[1:]) < 3:
        print ('Missing arguments')
        print ("Usage:\n %s [Mask] [SexCatalog] [Number] [Scale optional]" % sys.argv[0])
        print ("Example:\n %s mask.fits sex.cat 1 1 " % sys.argv[0])
        sys.exit()

    image = sys.argv[1]
    sexcat = sys.argv[2]
    number = sys.argv[3]

    if len(sys.argv[1:])  == 4:

        scale  = sys.argv[4]
    else:
        scale = 1

    number = int(number)


#    scale =2  #change this


    (NCol, NRow) = GetAxis(image)


    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(sexcat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    Rkron = scale * ai * kr

    Bim = (1 - e) * Rkron


    for idx, val in enumerate(n):

        if n[idx] == number:

            print ("removing elipse \n")
            (sxmin, sxmax, symin, symax) = GetSize(xx[idx], yy[idx], Rkron[idx], theta[idx], e[idx], NCol, NRow)


            if sxmin > 1:
                sxmin=np.int(np.round(sxmin-1))
            if sxmax < NCol:
                sxmax=np.int(np.round(sxmax+1))
            if symin > 1:
                symin=np.int(np.round(symin-1))
            if symax < NRow:
                symax=np.int(np.round(symax+1))


            SetZero(image, number, sxmin, sxmax, symin, symax)  # remove main object from mask


def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow



def GetSize(x, y, R, theta, ell, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    # k Check
    q = (1 - ell)
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

# getting size

    xmin = x - np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    xmax = x + np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    ymin = y - np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    ymax = y + np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    mask = xmin < 1
    if mask.any():
        xmin = 1

    mask = xmax > ncol
    if mask.any():
        xmax = ncol

    mask = ymin < 1
    if mask.any():
        ymin = 1

    mask = ymax > nrow
    if mask.any():
        ymax = nrow

    return (xmin, xmax, ymin, ymax)



def SetZero(image, num, xlo, xhi, ylo, yhi):
    "remove object from mask"
    # put all pixels
    # indicated in file to zero
    # in the image provided

#    my ($pixx,$pixy,$pixs,$errno,$flag,$subx,$suby);
#    my ($nam,$num,$temp);
#    flag=0



#############  ATENTION  ###############

#  Correct when xlo,xhi, ylo or yhi reach the top or bottom of the image

#########################################

    hdu = fits.open(image)
    dat = hdu[0].data

    mask = dat[ylo - 1:yhi, xlo - 1:xhi] == num

    if mask.any():
        dat[ylo - 1:yhi, xlo - 1:xhi][mask] = 0

    hdu[0].data = dat
    hdu.writeto(image, overwrite=True)
    hdu.close()

    return True

#end of program
if __name__ == '__main__':
    main()
