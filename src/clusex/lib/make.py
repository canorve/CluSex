#! /usr/bin/env python

import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
import matplotlib.colors as colors

import argparse
import subprocess as sp


from clusex.lib.check import CheckSatReg2
from clusex.lib.check import GetAxis 
from clusex.lib.check import CheckFlag 
from clusex.lib.check import GetWCS 
from clusex.lib.check import GetCounts
from clusex.lib.ds9 import ds9kron




def MakeMask(maskimage, catfile, scale, offset, regfile):
    "Create a mask image using ellipses for every Object of catfile. Now includes offset"
# k Check

    checkflag = 0
    flagsat = 4  # flag value when object is saturated (or close to)

    regflag = 0  # flag for saturaded regions

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, sxmin, sxmax, symin, symax, sxsmin, sxsmax, sysmin, sysmax = np.genfromtxt(
        catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    print("Creating Masks ... \n")

    Rkron = scale * ai * kr + offset

    mask = Rkron < 1
    if mask.any():
        Rkron[mask] = 1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    count1 = count2 = 0

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = CheckFlag(flg[idx], flagsat)
        # check if object is inside of a saturaded box region indicated by
        # user in ds9
        #regflag = CheckSatReg(xx[idx], yy[idx], Rkron[idx], theta[idx], e[idx],regfile)
        regflag = CheckSatReg2(xx[idx], yy[idx],regfile)



        if (checkflag == False) and (regflag == False):

            count1 += 1
            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rkron[idx], theta[idx], e[
                idx], sxsmin[idx], sxsmax[idx], sysmin[idx], sysmax[idx])

        elif(checkflag == True or regflag == True):
            count2 += 1
            

    print ("Number of ellipse masks created: {}  ".format(count1))
    print ("Number of objects rejected because they are saturated: {} \n".format(count2))

    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeKron(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

# Check

    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist < dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


def MakeSatBox(maskimage, region, val, ncol, nrow):
    "Create a mask for saturated regions"
    "Regions must be in DS9 box regions format"

# k Check

#	fileflag=1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    with open(region) as f_in:

        next(f_in)
#        next(f_in)
#        next(f_in)

        # All lines including the blank ones
        lines = (line.rstrip() for line in f_in)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            (box, info) = line.split('(')

            if(box == "box"):

                (xpos, ypos, xlong, ylong, trash) = info.split(',')

                xpos = float(xpos)
                ypos = float(ypos)
                xlong = float(xlong)
                ylong = float(ylong)

                xlo = (xpos - xlong / 2)
                xhi = (xpos + xlong / 2)

                ylo = (ypos - ylong / 2)
                yhi = (ypos + ylong / 2)

                xlo = int(xlo)
                xhi = int(xhi)

                ylo = int(ylo)
                yhi = int(yhi)

                if (xlo < 1):

                    xlo = 1

                if (xhi > ncol):

                    xhi = ncol

                if (ylo < 1):

                    ylo = 1

                if (yhi > nrow):

                    yhi = nrow

                img[ylo - 1:yhi, xlo - 1:xhi] = val

    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True




def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"
# k Check
    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex))
    hdu.writeto(newfits, overwrite=True)

    return True


def CatArSort(SexCat,scale,offset,SexArSort,NCol,NRow):
    # k Check

    # sort the sextractor
    # catalog by magnitude,
    # get sizes for objects
    # and write it in a new file

    # The sextractor catalog must contain the following parameters:
    #   1 NUMBER                 Running object number
    #   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
    #   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
    #   4 X_IMAGE                Object position along x                                    [pixel]
    #   5 Y_IMAGE                Object position along y                                    [pixel]
    #   6 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   7 KRON_RADIUS            Kron apertures in units of A or B
    #   8 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
    #   9 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
    #  10 A_IMAGE                Profile RMS along major axis                               [pixel]
    #  11 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
    #  12 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
    #  13 BACKGROUND             Background at centroid position                            [count]
    #  14 CLASS_STAR             S/G classifier output
    #  15 FLAGS                  Extraction flags


    print("Sorting and getting sizes for objects \n")

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(
        SexCat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    Rkron = scale * ai * kr + offset

    Rwsky = scale * ai * kr + 10  + 20


#   considering to use only  KronScale instead of SkyScale
#    Rwsky = parvar.KronScale * ai * kr + parvar.Offset + parvar.SkyWidth

    Bim = (1 - e) * Rkron

    Area = np.pi * Rkron * Bim *(-1)

    (sxmin, sxmax, symin, symax) = GetSize(xx, yy, Rkron, theta, e, NCol, NRow)

    (sxsmin, sxsmax, sysmin, sysmax) = GetSize(xx, yy, Rwsky, theta, e, NCol,NRow)


    f_out = open(SexArSort, "w")

    index = Area.argsort()
    for i in index:

        line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(n[i], alpha[i], delta[i], xx[i], yy[i], mg[i], kr[i], fluxrad[i], ia[i], ai[i], e[i], theta[i], bkgd[i], idx[i], flg[i], int(
            np.round(sxmin[i])), int(np.round(sxmax[i])), int(np.round(symin[i])), int(np.round(symax[i])), int(np.round(sxsmin[i])), int(np.round(sxsmax[i])), int(np.round(sysmin[i])), int(np.round(sysmax[i])))

        f_out.write(line)

    f_out.close()

    return len(n)


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
        xmin[mask] = 1

    mask = xmax > ncol
    if mask.any():
        xmax[mask] = ncol

    mask = ymin < 1
    if mask.any():
        ymin[mask] = 1

    mask = ymax > nrow
    if mask.any():
        ymax[mask] = nrow

    xmin = np.rint(xmin)
    ymin = np.rint(ymin)
    xmax = np.rint(xmax)
    ymax = np.rint(ymax)

    return (xmin.astype(int), xmax.astype(int), ymin.astype(int), ymax.astype(int))


def EraseObjectMask(MaskFile,tempMask,obj,fill = 0):

    hdumask = fits.open(MaskFile)
    data = hdumask[0].data

    mask = data == obj

    data[mask] = fill # removing object from mask
    
    hdumask[0].data = data

    hdumask.writeto(tempMask, overwrite=True)
    
    hdumask.close()
    
def EraseObjectMask2(maskimg,obj):

    tempmask=maskimg.copy()

    mask = tempmask == obj

    tempmask[mask] = 0 # removing object from mask
   
    return tempmask


def MakeStamps(image, catalog, maskimage, stretch, skyoff, dpi, 
                cmap, scale, offset, bright, contrast, frac, fracmax, galclass):


    hdu = fits.open(image)
    img = hdu[0].data
    hdu.close()


    hdu = fits.open(maskimage)
    mask = hdu[0].data
    hdu.close()


    NCol, NRow = GetAxis(image)


    wcs = GetWCS(image)
    #counts = GetCounts(image)


    STRETCH_CONST = stretch #for stamps sizes

    flagsat=4      ## flag value when object is saturated (or close to)

    #creates a object image: same image but counts = 0 for regions where 
    #there is no objets
    objimage = MakeObjImg(img,mask) 



    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Class,Flg=np.genfromtxt(catalog,delimiter="",unpack=True)


    if not os.path.exists("stamps"):
        print("Creating directory for stamps ... ")
        os.makedirs("stamps")



    Rkron = scale * Ai * Kr + offset

    maskron = Rkron == 0 
    if maskron.any():
        Rkron[maskron] = 1


    Rkron = Rkron * STRETCH_CONST    

    q = (1 - E)
    bim = q * Rkron


    (xmin, xmax, ymin, ymax) = GetSize(X, Y, Rkron, Theta, E, NCol, NRow) #size of every object


    line="Creating image stamps for every object of the catalog"
    print (line)


    num = 0
    num_stamps = 0
    for idx, item in enumerate(N):


        check = CheckFlag(Flg[idx],flagsat)  ## check if object doesn't has saturated regions


        if  ((check == False) and (Class[idx] < galclass)):


            objimg = objimage.copy() # a new objimg for every new stamp 


            objmask = N[idx] == mask
            objimg[objmask] = 0 # removing main object from object image 
     
            objmask2 = (mask != N) & (mask != 0)

            # removes all the objects except the main target.
            stamp = img - objimg 

            imgstmp = "obj-" + str(round(N[idx])) + ".png"

            yy = int(X[idx]) - xmin[idx]  #interchange because numpy arrays 

            xx = int(Y[idx]) - ymin[idx] 

            #quick routine:
            #ShowImg(stamp[ymin[idx]-1:ymax[idx]-1,xmin[idx]-1:xmax[idx]-1], 
            #            xx, yy, wcs, imgstmp, dpival = dpi, sky = Bkgd[idx], 
            #            cmap = cmap, bri = bright, con = contrast, 
            #            frac = frac, fracmax = fracmax)

            ShowImg(stamp[ymin[idx]-1:ymax[idx]-1,xmin[idx]-1:xmax[idx]-1], 
                        wcs, imgstmp, dpival = dpi, sky = Bkgd[idx], 
                        cmap = cmap, bri = bright, con = contrast, 
                        frac = frac, fracmax = fracmax)



            #slow routine

            #GetPng(stamp[ymin[idx]-1:ymax[idx]-1,xmin[idx]-1:xmax[idx]-1], 
            #        counts, wcs, dpi=dpi, cmap = cmap, namepng = imgstmp, 
            #        bri = bright, con = contrast)


            #move to folder change this for an function of os library
            runcmd = "mv  {}  stamps/{}".format(imgstmp,imgstmp)
            errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                               stderr=sp.PIPE, universal_newlines=True)

            num_stamps +=1
        else:

            num +=1

    print ("{} objects rejected because they have saturated pixels or are classified as stars  \n".format(num))
    print ("{} stamps created  \n".format(num_stamps))


def MakeObjImg(image,mask):
    """make object image """
    objimg = image.copy()

    zero_pix = mask == 0 

    objimg[zero_pix]  = 0

    return objimg


def ShowImg(img: np.array ,xc: int, yc: int, wcs, namepng="obj.png", 
            dpival=100, sky=1, cmap='viridis', bri = 0, con = 1, 
            frac = 1, fracmax = 1):
    """This routine shows the image"""

        
    #hdu = fits.open(cubeimg)
    #data = (hdu[1].data.copy()).astype(float)
    #hdu.close()

    data =  (img.copy()).astype(float)

    root_ext = os.path.splitext(namepng)

    objname = root_ext[0]

    #flatten image

    flatdatimg = data.flatten()  

    flatdatimg.sort()

    datimgpatch = flatdatimg#[modbot:modtop]

    datmin = np.min(datimgpatch)
    datmax = np.max(datimgpatch)


    data = data.clip(datmax/1e4, datmax)

    datmin = datmax/1e4


    middle = (datmax -datmin)/2

   #brightness auto-adjust according to the contrast value 

    Autobri = middle*(con -1) + modmin*(1-con) 

    #user can re-adjust according to the contrast value
    newdata = con*(data - middle) + middle + Autobri + bri*(datmax-middle)

    mask = data < 0 
    data[mask] = 1 # avoids problems in log
     
    fig, ax1 = plt.subplots(figsize=(8, 8)) 
    #fig, ax1 = plt.subplots(figsize=(700/dpival, 700/dpival),dpi=dpival) 
    #plt.figure(figsize = (800 / my_dpi, 800 / my_dpi), dpi = my_dpi)


    fig.subplots_adjust(left=0.08, right=0.94, bottom=0.04, top=0.94)

    ax1 = fig.add_subplot(projection=wcs)

    #flatdata = data.flatten()  

    #flatdata.sort()

    #tot=len(flatdata)

    #top=round(.9*tot)
    #bot=round(.1*tot)

    #imgpatch=flatdata#[bot:top]


    #galmin = np.min(imgpatch)
    #galmin = sky 
    #galmax = np.max(imgpatch)

    #median=np.median(imgpatch)

    #galmin = (frac)*galmin 
    #galmax = fracmax*galmax


    #if (galmin > galmax): #to prevent from failing
    #    galmin, galmax = galmax, galmin

    #my old routine
    #ax1.imshow(con*data+bri, origin = 'lower', norm
    #            = colors.LogNorm(vmin = galmin, vmax = galmax), 
    #            cmap = cmap, interpolation='nearest')


    im = ax1.imshow(newdata, origin ='lower', interpolation='nearest', norm 
                    = colors.LogNorm(vmin=datmin, vmax=datmax), cmap = cmap)


    ax1.set_xlabel('Right Ascension')
    ax1.set_ylabel('Declination')
    ax1.grid(color='black', ls='solid', alpha=0.1)

    ax1.set_title(objname)
    plt.savefig(namepng,dpi=dpival)
     
    plt.close()


def GetPng(data, counts, wcs, dpi=200, cmap='gray_r',namepng="obj.png", bri = 33, con = 0.98):
    "Converts image into a PNG image with axis coordinates, inverted colormap, log/zmax style"

    #Slow routine, it was avoided

    #filename = get_pkg_data_filename(Image)

    #hdu = fits.open(filename)[0]

    n, bins, patches = plt.hist(counts, bins=512, range=(min(counts), max(counts)), color='black')
    n_max = np.argsort(n)[::-1]
    n_idx = n_max[0:20]
    lim_min = min(bins[n_idx])
    lim_max = max(bins[n_idx])

    #bri = 33 # brightness, source: docs.opencv.org/3.4/d3/dc1/tutorial_basic_linear_transform.html
    #con = 0.98 # contrast, > 0

    plt.subplot(projection=wcs)
    plt.imshow(con*data+bri, cmap=cmap, norm=LogNorm(lim_min, lim_max))

    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.grid(color='black', ls='solid', alpha=0.1)
    plt.title(namepng)
    plt.tight_layout()
    #plt.savefig('%s.png' % (Image), bbox_inches='tight', dpi=dpi)
    plt.savefig(namepng, bbox_inches='tight', dpi=dpi)


