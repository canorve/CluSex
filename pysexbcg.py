#! /usr/bin/env python

import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path
import scipy

# This program creates a catalog of Sextractor with
# a combination of two runs of Sextractor with
# different configuration parameters.

# pybcg runs sextractor for the third time to extract Bright cluster Galaxies

# ATENTION
#  fusionar las librerias con DGCG


def main():

    if len(sys.argv[1:]) != 2:
        print ('Missing arguments')
        print ("Usage:\n %s [ConfigFile] [ImageFile]" % sys.argv[0])
        print ("Example:\n %s Config.txt image.fits " % sys.argv[0])
        sys.exit()

    ConfigFile= sys.argv[1]
    image= sys.argv[2]

# default parameters
# init parameters
    outhot  = "hot.sex"
    outcold = "cold.sex"
    sexfile="default.sex"
#    image="A1413-cD.fits"
    output="hc.cat"
    output2="hotcold.cat"
    scale=1

    run1=run2=run3=0

    create = 0

    regoutfile = "hotcold.reg"

    maskfile="mask.fits"

    SexArSort="sortar.cat"

    outbcg = "bcg.sex"


    dn1 = dn2 = dn3 = 1
    dm1 = dm2 = dm3 = 1
    at1 = at2 = at3 = 1
    dt1 = dt2 = dt3 = 1
    da1 = da2 = da3 = 1
    bs1 = bs2 = bs3 = 1
    bf1 = bf2 = bf3 = 1



    if not os.path.exists(ConfigFile):
        print ('%s: filename does not exist!' %sys.argv[1])
        sys.exit()

    count=0


    with open(ConfigFile) as f_in:

        lines = (line.rstrip() for line in f_in) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) #remove comments
        lines = (line for line in lines if line) # Non-blank lines



        for line2 in lines:

            (param,val)=line2.split()


# param 1
            if param == "FirstRun":

                (run1) = val.split()[0]
                run1=int(run1)



            if param == "DEBLEND_NTHRESH1":
                (dn1) = val.split()[0]
                dn1=int(dn1)

            if param == "DEBLEND_MINCONT1":
                (dm1) = val.split()[0]
                dm1=float(dm1)


            if param == "ANALYSIS_THRESH1":

                (at1) = val.split()[0]
                at1=float(at1)


            if param == "DETECT_THRESH1":

                (dt1) = val.split()[0]
                dt1= float(dt1)

            if param == "DETECT_MINAREA1":

                (da1) = val.split()[0]
                da1=int(da1)

            if param == "BACK_SIZE1":

                (bs1) = val.split()[0]
                bs1= float(bs1)

            if param == "BACK_FILTERSIZE1":

                (bf1) = val.split()[0]
                bf1=int(bf1)


#  param2

            if param == "SecondRun":

                (run2) = val.split()[0]
                run2=int(run2)


            if param == "DEBLEND_NTHRESH2":

                (dn2) = val.split()[0]
                dn2=float(dn2)


            if param == "DEBLEND_MINCONT2":

                (dm2) = val.split()[0]
                dm2=float(dm2)


            if param == "ANALYSIS_THRESH2":

                (at2) = val.split()[0]
                at2=float(at2)


            if param == "DETECT_THRESH2":

                (dt2) = val.split()[0]
                dt2=float(dt2)

            if param == "DETECT_MINAREA2":

                (da2) = val.split()[0]
                da2=int(da2)

            if param == "SatDs9":
                (satfileout) = val.split()[0]
                satfileout=str(satfileout)

            if param == "SatScale":
                (satscale) = val.split()[0]
                satscale=float(satscale)

            if param == "SatOffset":
                (satoffset) = val.split()[0]
                satoffset=int(satoffset)

            if param == "BACK_SIZE2":

                (bs2) = val.split()[0]
                bs2= float(bs2)

            if param == "BACK_FILTERSIZE2":

                (bf2) = val.split()[0]
                bf2=int(bf2)

            if param == "MakeMask":

                (create) = val.split()[0]
                create=int(create)

            if param == "Scale":

                (scale) = val.split()[0]
                scale=float(scale)

            if param == "OutCatalog":
                (output2) = val.split()[0]
                output2=str(output2)

            if param == "RegDs9":
                (regoutfile) = val.split()[0]
                regoutfile=str(regoutfile)

# bcg.sex
# param 3

            if param == "BGCRun":

                (run3) = val.split()[0]
                run3=int(run3)


            if param == "BCG":

                (bcg) = val.split()[0]
                bcg=int(bcg)

            if param == "DEBLEND_NTHRESH3":

                (dn3) = val.split()[0]
                dn3=float(dn3)



            if param == "DEBLEND_MINCONT3":

                (dm3) = val.split()[0]
                dm3=float(dm3)


            if param == "ANALYSIS_THRESH3":

                (at3) = val.split()[0]
                at3=float(at3)


            if param == "DETECT_THRESH3":

                (dt3) = val.split()[0]
                dt3=float(dt3)

            if param == "DETECT_MINAREA3":

                (da3) = val.split()[0]
                da3=int(da3)


            if param == "BACK_SIZE3":

                (bs3) = val.split()[0]
                bs3= float(bs3)

            if param == "BACK_FILTERSIZE3":

                (bf3) = val.split()[0]
                bf3=int(bf3)

######

    flaghead = False




    f_out = open(outhot, "w")

    with open(sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

#every line of the file:
        for line2 in lines:

        #print line2
            (params)=line2.split()

            if params[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(dn1)

            if params[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(dm1)

            if params[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(at1)

            if params[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(dt1)

            if params[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(da1)

            if params[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME hot.cat"

            if params[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(bs1)

            if params[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(bf1)

######  check output catalog
            if params[0] ==  "CATALOG_TYPE":
                if params[1] !=  "ASCII":
	                   flaghead = True
######




            line2 = line2+"\n"
            f_out.write(line2)


    f_out.close()

#######
# cold.sex

    f_out = open(outcold, "w")

    with open(sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

#every line of the file:
        for line2 in lines:

            (params)=line2.split()

            if params[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(dn2)

            if params[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(dm2)

            if params[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(at2)

            if params[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(dt2)

            if params[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(da2)

            if params[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME cold.cat"

            if params[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(bs2)

            if params[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(bf2)


            f_out.write(line2+"\n")


    f_out.close()


######

    f_out = open(outbcg, "w")

    with open(sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

#every line of the file:
        for line2 in lines:


            (params)=line2.split()

            if params[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(dn3)

            if params[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(dm3)

            if params[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(at3)

            if params[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(dt3)

            if params[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(da3)

            if params[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME bcg.cat"

            if params[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(bs3)

            if params[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(bf3)


            line2 = line2+"\n"
            f_out.write(line2)


    f_out.close()



    if flaghead:
        print ("CATALOG_TYPE must be ASCII in default.sex. Ending program.. \n")
        sys.exit()


    (NCol, NRow) = GetAxis(image)

    print("Running hot.sex, cold.sex and bcg.sex  \n")


######  Running Sextractor  #######

    if (run1 == 1):

        runcmd="sextractor -c {} {} ".format(outhot,image)
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT


    if (run2 == 1):

        runcmd="sextractor -c {} {} ".format(outcold,image)
        err2 = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT


    if (run3 == 1):

        runcmd="sextractor -c {} {} ".format(outbcg,image)
        err3 = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT

####################################

    if (run1 == 1 and run3 ==1):
        print ("adding the bright cluster galaxies from bcg.cat to hot.cat catalogs ....\n");

        Total = CatArSort("bcg.cat",scale,"bcgsort.cat",NCol,NRow)

        addbcgcat("hot.cat","bcgsort.cat","hotbcg.cat",bcg,scale)

        if (run2 == 1):
            print ("joining hotbcg.cat and cold.cat catalogs ....\n");
            joinsexcat("hotbcg.cat","cold.cat",output,scale)

    elif(run2 == 1 and run3 ==1):
        print ("adding the bright cluster galaxies from bcg.cat to cold.cat catalogs ....\n");

        Total = CatArSort("bcg.cat",scale,"bcgsort.cat",NCol,NRow)

        addbcgcat("cold.cat","bcgsort.cat","coldbcg.cat",bcg,scale)


###############

    if (run1 == 1 and run2 ==1 and run3 == 1 ):
        print ("creating {0} for ds9 ....\n".format(satfileout))
        ds9satbox(satfileout,output,satscale,satoffset) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat(output, output2,satfileout)


    elif(run1 == 1 and run3 == 1):
        print ("creating {0} for ds9 ....\n".format(satfileout))
        ds9satbox(satfileout,"hotbcg.cat",satscale,satoffset) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("hotbcg.cat", output2,satfileout)


    elif(run2 ==1 and run3 == 1):
        print ("creating {0} for ds9 ....\n".format(satfileout))
        ds9satbox(satfileout,"coldbcg.cat",satscale,satoffset) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("coldbcg.cat", output2,satfileout)


    elif(run1 ==1):
        print ("creating {0} for ds9 ....\n".format(satfileout))
        ds9satbox(satfileout,"hot.cat",satscale,satoffset) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("hot.cat", output2,satfileout)



    elif(run2 ==1):
        print ("creating {0} for ds9 ....\n".format(satfileout))
        ds9satbox(satfileout,"cold.cat",satscale,satoffset) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("cold.cat", output2,satfileout)


    elif(run3 ==1):
        print ("creating {0} for ds9 ....\n".format(satfileout))
        ds9satbox(satfileout,"bcg.cat",satscale,satoffset) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("bcg.cat", output2,satfileout)


    print ("{0} is the output catalog  ....\n".format(output2))


    print ("Creating ds9 check region files....\n")


#    ds9kron("bcg.cat","bcg.reg",scale)

#    ds9kron("hot.cat","hot.reg",scale)

#    ds9kron("cold.cat","cold.reg",scale)


#####
    if (run1 == 1 or run2 ==1 or run3 == 1 ):
        ds9kron(output2,regoutfile,scale)
#####

    if (create == 1):

        if (run1 == 1 or run2 ==1 or run3 == 1 ):

            print ("Creating masks....\n")

#        print(output2,scale,SexArSort,NCol,NRow)
            Total = CatArSort(output2,scale,SexArSort,NCol,NRow)

#        ParVar.Total = catfil.CatSort(ParVar)


##### segmentation mask

            MakeImage(maskfile, NCol, NRow)

            MakeMask(maskfile, SexArSort, scale,0,satfileout)  # offset set to 0
            MakeSatBox(maskfile, satfileout, Total + 1, NCol, NRow)

            print ("Running ds9 ....\n")
            runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} {} ".format(image,regoutfile,satfileout,maskfile)
            err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT

        else:
            print("can't create mask because sextractor was not used \n")


    else:

        if (run1 == 1 or run2 ==1 or run3 == 1 ):

            print ("Running ds9 ....\n")
            runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} ".format(image,regoutfile,satfileout)
            err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT

        else:
            print("can't create mask because sextractor was not used \n")


#############################################################################
######################### End of program  ######################################


#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/


##############################################################################
########################## Functions ###########################################

def MakeMask(maskimage, catfile, scale, offset, regfile):
    "Create a mask image using ellipses for every Object of catfile. Now includes offset"
# k Check

    checkflag = 0
    flagsat = 4  # flag value when object is saturated (or close to)
    maxflag = 128  # max value for flag

    regflag = 0  # flag for saturaded regions

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, sxmin, sxmax, symin, symax, sxsmin, sxsmax, sysmin, sysmax = np.genfromtxt(
        catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    print("Creating Masks for sky \n")

    Rkron = scale * ai * kr + offset

    mask = Rkron < 1
    if mask.any():
        Rkron[mask] = 1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = CheckFlag(flg[idx], flagsat, maxflag)
        # check if object is inside of a saturaded box region indicated by
        # user in ds9
        regflag = CheckSatReg(xx[idx], yy[idx], Rkron[idx], theta[idx], e[idx], regfile)

        if (checkflag == False) and (regflag == False):

            print ("Creating ellipse mask for object {}  ".format(n[idx]))
            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rkron[idx], theta[idx], e[
                idx], sxsmin[idx], sxsmax[idx], sysmin[idx], sysmax[idx])

        elif(checkflag == True or regflag == True):

            print ("Skipping object {}, one or more pixels are saturated ".format(n[idx]))

    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeKron(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

# Check

    xmin = np.int(xmin)
    xmax = np.int(xmax)
    ymin = np.int(ymin)
    ymax = np.int(ymax)

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


def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow


def CatArSort(SexCat,scale,SexArSort,NCol,NRow):
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

    Rkron = scale * ai * kr

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

        line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(n[i], alpha[i], delta[i], xx[i], yy[i], mg[i], kr[i], fluxrad[i], ia[i], ai[i], e[i], theta[i], bkgd[i], idx[i], flg[i], np.int(
            np.round(sxmin[i])), np.int(np.round(sxmax[i])), np.int(np.round(symin[i])), np.int(np.round(symax[i])), np.int(np.round(sxsmin[i])), np.int(np.round(sxsmax[i])), np.int(np.round(sysmin[i])), np.int(np.round(sysmax[i])))

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

    return (xmin, xmax, ymin, ymax)


def joinsexcat (maincat,secondcat,output,KronScale):
    "merges two Sextractor catalogs"

    f_out = open(output, "w")

    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(maincat,delimiter="",unpack=True)

    AR         = 1 - E
    RKron      = KronScale * Ai * Kr



    maskron = RKron <= 0
    RKron[maskron]=1

    maskar = AR <= 0.005

    AR[maskar]=0.005

    for idx, item in enumerate(N):

        line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

        f_out.write(line)


    total =len(N)

    NewN=total + 1


#second cat
    N2,Alpha2,Delta2,X2,Y2,Mg2,Kr2,Fluxr2,Isoa2,Ai2,E2,Theta2,Bkgd2,Idx2,Flg2=np.genfromtxt(secondcat,delimiter="",unpack=True)

    AR2         = 1 - E2
    RKron2      = KronScale * Ai2 * Kr2


    for idx2, item2 in enumerate(N2):

        flag =False
        for idx, item in enumerate(N):

#            flag=CheckKron(X2[idx2],Y2[idx2],X[idx],Y[idx],RKron[idx],Theta[idx],AR[idx])
            flag=CheckKron(X[idx],Y[idx],X2[idx2],Y2[idx2],RKron2[idx2],Theta2[idx2],AR2[idx2])

            if flag:   # boolean value
                break

        if not flag:

            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(NewN, Alpha2[idx2], Delta2[idx], X2[idx2], Y2[idx2], Mg2[idx2], Kr2[idx2], Fluxr2[idx2], Isoa2[idx2], Ai2[idx2], E2[idx2], Theta2[idx2], Bkgd2[idx2], Idx2[idx2], Flg2[idx2])

            f_out.write(line)

            NewN+=1

    f_out.close()


def addbcgcat (maincat,secondcat,output,bcg,KronScale):
    "add Bright cluster galaxies to second catalog two Sextractor catalogs"


    flagsat=4      ## flag value when object is saturated (or close to)
    maxflag=128    ## max value for flag
    checkflag=0


    f_out = open(output, "w")


    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(maincat,delimiter="",unpack=True)

    AR         = 1 - E
    RKron      = KronScale * Ai * Kr


    maskron = RKron <= 0
    RKron[maskron]=1

    maskar = AR <= 0.005

    AR[maskar]=0.005


#second cat bcg

    N2,Alpha2,Delta2,X2,Y2,Mg2,Kr2,Fluxr2,Isoa2,Ai2,E2,Theta2,Bkgd2,Idx2,Flg2,xmin2,xmax2,ymin2,ymax2,xsmin2,xsmax2,ysmin2,ysmax2=np.genfromtxt(secondcat,delimiter="",unpack=True)

    idx2=0
    while idx2 < bcg:

        flag =False

        distmin=10

        checkflag=CheckFlag(Flg2[idx2],flagsat,maxflag)

        if checkflag == False:

            for idx, item in enumerate(N):

                dist = np.sqrt((X2[idx2] - X[idx])**2 + (Y2[idx2] - Y[idx])**2)

                if dist < distmin:
                    distmin=dist
                    flag = True
                    sidx=idx

            if flag:

                print("BGC encontrado")
                Alpha[sidx] = Alpha2[idx2]
                Delta[sidx] = Delta2[idx2]
                X[sidx] = X2[idx2]
                Y[sidx] = Y2[idx2]
                Mg[sidx] = Mg2[idx2]
                Kr[sidx] = Kr2[idx2]
                Fluxr[sidx] = Fluxr2[idx2]
                Isoa[sidx] = Isoa2[idx2]
                Ai[sidx] = Ai2[idx2]
                E[sidx] = E2[idx2]
                Theta[sidx]  = Theta2[idx2]
                Bkgd[sidx] = Bkgd2[idx2]
                Idx[sidx] = Idx2[idx2]
                Flg[sidx] = Flg2[idx2]

        else:
            bcg+=1

        idx2 +=1


## writing results

    for idx, item in enumerate(N):

        line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


        f_out.write(line)


    total =len(N)

    NewN=total + 1


    f_out.close()



def CheckKron (xpos,ypos,x,y,R,theta,q):
    "check if position is inside of the Kron Ellipse saturaded region returns True if object center is in Ellipse region"


    bim = q * R

    theta = theta * np.pi /180  ## Rads!!!


    flag =False


    dx = xpos - x
    dy = ypos - y


    landa = np.arctan2( dy,dx )

    if landa < 0:
        landa=landa + 2 * np.pi


    landa = landa - theta


    angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

    xell =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
    yell =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

    dell = np.sqrt ( (xell - x)**2 + (yell - y)**2 )
    dist = np.sqrt ( dx**2 + dy**2 )


    if  dist < dell:
	      flag=True


    return flag



def CheckFlag(val,check,max):
   "Check for flag contained in $val, returns 1 if found "

   flag = False
   mod = 1

   while (mod != 0):


       res = int(val/max)

       if (max == check and res == 1 ):

           flag=True


       mod = val % max

       val = mod
       max = max/2



   return flag



def ds9satbox (satfileout,output,satscale,satoffset):
    "Creates a file for ds9 which selects bad saturated regions"

    scaleflag=1
    offsetflag=1
    regfileflag=1
    magflag=1
    clasflag=1

    flagsat=4      ## flag value when object is saturated (or close to)
    maxflag=128    ## max value for flag
    check=0
    regflag = 0    ## flag for saturaded regions




    f_out = open(satfileout, "w")

    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(output,delimiter="",unpack=True)


    line="image \n"
    f_out.write(line)



    for idx, item in enumerate(N):


        if (E[idx] >= 0.3):
             E[idx]-=0.2

        bi=Ai[idx]*(1-E[idx])

        Theta[idx] = Theta[idx] * np.pi /180  #rads!!!


        Rkronx = satscale * 2 * Ai[idx] * Kr[idx]  + satoffset
        Rkrony = satscale * 2 * bi * Kr[idx]  + satoffset



        if Rkronx == 0:
            Rkronx = 1

        if Rkrony == 0:
            Rkrony = 1

        check=CheckFlag(Flg[idx],flagsat,maxflag)  ## check if object has saturated regions

        if (check):

            line="box({0},{1},{2},{3},0) # color=red move=0 \n".format(X[idx],Y[idx],Rkronx,Rkrony)
            f_out.write(line)

            line2="point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={{ {2} }} \n".format(X[idx],Y[idx],N[idx])
            f_out.write(line2)


    f_out.close()




def putflagsat(sexfile,sexfile2,regfile):
    "Put flags on objects which are inside saturated regions"


    f_out= open(sexfile2, "w")

    scale = 1
    offset=0


    flagsat=4      ## flag value when object is saturated (or close to)
    maxflag=128    ## max value for flag


    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(sexfile,delimiter="",unpack=True)

    for idx, item in enumerate(N):

        Rkron = scale * Ai[idx] * Kr[idx] + offset

        if Rkron == 0:

            Rkron = 1

        q = (1 - E)
        bim = q * Rkron


        check=CheckFlag(Flg[idx],flagsat,maxflag)  ## check if object doesn't has saturated regions
        regflag=CheckSatReg(X[idx],Y[idx],Rkron,Theta[idx],E[idx],regfile)    ## check if object is inside of a saturaded box region indicated by user in ds9

        if  (check == False ) and ( regflag == True) :

            Flg[idx] = Flg[idx] + 4


            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

            f_out.write(line)


        else:

            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


            f_out.write(line)



    f_out.close()



def ds9kron(sexfile,regfile,scale):
    "Creates ds9 region file to check output catalog "


    f_out= open(regfile, "w")

#    scale = 1
    offset=0


    flagsat=4      ## flag value when object is saturated (or close to)
    maxflag=128    ## max value for flag



#print OUT "image \n";


    line="image \n"
    f_out.write(line)


    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Star,Flg=np.genfromtxt(sexfile,delimiter="",unpack=True)

    for idx, item in enumerate(N):

        Rkron = scale * Ai[idx] * Kr[idx] + offset

#        print (Rkron)

        if Rkron == 0:

            Rkron = 1

        q = (1 - E)
        bim = q * Rkron


        check=CheckFlag(Flg[idx],flagsat,maxflag)  ## check if object doesn't has saturated regions
#        regflag=CheckSatReg(X[idx],Y[idx],Rkron,Theta[idx],E[idx],regfile)    ## check if object is inside of a saturaded box region indicated by user in ds9

        if  (check == False ) :

            line="ellipse({0},{1},{2},{3},{4}) # color=blue move=0 \n".format(X[idx],Y[idx],Rkron,bim[idx],Theta[idx])

            f_out.write(line)


            line2="point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={2} {3} {4} \n".format(X[idx],Y[idx],"{",int(N[idx]),"}")

            f_out.write(line2)



        else:

            print ("Skipping object {} one or more pixels are saturated ".format(N[idx]))


    #        f_out.write(line)



    f_out.close()





def CheckSatReg(x,y,R,theta,ell,filein):
   "Check if object is inside of saturated region. returns True if at least one pixel is inside"
## check if object is inside of
## saturaded region as indicated by ds9 box region
## returns 1 if object center is in saturaded region


   q = (1 - ell)

   bim = q * R

   theta = theta * np.pi /180  ## Rads!!!

   flag = False
   fileflag =1



   with open(filein) as f_in:

       lines = (line.rstrip() for line in f_in) # All lines including the blank ones
       lines = (line.split('#', 1)[0] for line in lines) # remove comments
       lines = (line.rstrip() for line in lines)   # remove lines containing only comments
       lines = (line for line in lines if line) # Non-blank lines

       for line in lines:


           if (line != "image"):

               (box,info)=line.split('(')

               if(box == "box"):

                   (xpos,ypos,xlong,ylong,trash)=info.split(',')

                   xpos=float(xpos)
                   ypos=float(ypos)
                   xlong=float(xlong)
                   ylong=float(ylong)


                   xlo = xpos - xlong/2
                   xhi = xpos + xlong/2

                   ylo = ypos - ylong/2
                   yhi = ypos + ylong/2

                   dx = xpos - x
                   dy = ypos - y

                   landa=np.arctan2( dy,dx )

                   if landa < 0:
                       landa=landa + 2 * np.pi


                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xell =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yell =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)


                   if ( (xell > xlo and xell < xhi) and (yell > ylo and yell < yhi)  ):

                       flag=True
                       break


   return flag

#end of program
if __name__ == '__main__':
    main()
