#! /usr/bin/env python

import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path
import scipy

from pathlib import Path
import shutil
import argparse


from clusex.lib.init import Params

from clusex.lib.read import readcon 
from clusex.lib.writesex import wsex
from clusex.lib.writesex import runsex

from clusex.lib.join import joinsexcat 


from clusex.lib.satbox  import Ds9SatBox
from clusex.lib.satbox  import putflagsat  
from clusex.lib.satbox  import GetAxis

from clusex.lib.ds9 import ds9kron


from clusex.lib.mask import MakeMask 
from clusex.lib.mask import MakeImage
from clusex.lib.mask import MakeSatBox 
from clusex.lib.mask import CatArSort

# This program creates a catalog of Sextractor with
# a combination of two runs of Sextractor with
# different configuration parameters


def main():

#########################################
##### lectura de archivo #######
#########################################


    parser = argparse.ArgumentParser(description="CluSex: combines two Sextractor catalogs among other stuff")

    # required arguments
    parser.add_argument("ConfigFile",help="CluSex configuration file ")
    parser.add_argument("image",help="FITS image file of the galaxy cluster ")


    args = parser.parse_args()

    confile = args.ConfigFile 
    image = args.image

#########################################
##### parametros iniciales #######
#########################################


    #copying files to actual folder

    filsex = Path(__file__).parent / "../../def/default.sex"
    filnnw = Path(__file__).parent / "../../def/default.nnw"
    filcon = Path(__file__).parent / "../../def/default.conv"
    filpar = Path(__file__).parent / "../../def/sex.param"

    to_des = Path('.')

    shutil.copy2(filsex, to_des)  
    shutil.copy2(filnnw, to_des)  
    shutil.copy2(filcon, to_des) 
    shutil.copy2(filpar, to_des)


# init parameters
# default values

    params=Params()

#########################################
##### lectura del archivo de parametros #######
#########################################

    readcon(params,confile)

##############################################################
##### lectura y modificacion de los archivos de Sextractor ########
#########################################################

    wsex(params)

##############################################################
##### ejecucion de los archivos de Sextractor ########
#########################################################



    runsex(params,image)


##############################################################
##### llamado a la funcion radcor, rad correccion################### 
#########################################################




##############################################################
##### une los dos archivos hot.cat y cold.cat ####### 
#########################################################


    if (params.run1 == 1 and params.run2 == 1):
        print ("joining hot.cat and cold.cat catalogs ....\n");
        joinsexcat("hot.cat","cold.cat",params.output,params.scale,params.scale2)

    else:
        print("Can not join catalogs because sextractor was not used \n")



##############################################################
##### busca y selecciona las estrellas saturadas ####### 
##### modifica las banderas para los objetos cercanos a regiones 
###### de saturacion ####
#########################################################
    


    if (params.run1 == 1 and params.run2 == 1):
        print ("creating {0} for ds9 ....\n".format(params.satfileout))
        Ds9SatBox(image,params.satfileout,params.output,params.satscale,params.satoffset,params.satlevel,params.minsatsize,params.satq) # crea archivo  de salida de reg

        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat(params.output,params.output2,params.satfileout)

    elif(params.run1 ==1):
        print ("creating {0} for ds9 ....\n".format(params.satfileout))
        Ds9SatBox(image,satfileout,"hot.cat",params.satscale,params.satoffset,params.satlevel,params.minsatsize,params.satq) # crea archivo  de salida de reg
        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("hot.cat","hot2.cat",params.satfileout)
        # renaming hot2.cat catalog
        runcmd="mv hot2.cat hot.cat"
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT


    elif(params.run2 == 1):
        print ("creating {0} for ds9 ....\n".format(params.satfileout))
        Ds9SatBox(image,params.satfileout,"cold.cat",params.satscale,params.satoffset,params.satlevel,params.minsatsize,params.satq) # crea archivo  de salida de reg
        print ("recomputing flags on objects which are inside saturated regions  ....\n")
        putflagsat("cold.cat","cold2.cat",params.satfileout)
        # renaming cold2.cat catalog
        runcmd="mv cold2.cat cold.cat"
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT




##################################################################
###### crea los archivos de salida para ds9 ###################### 
##################################################################


    if (params.run1 == 1 and params.run2 == 1):
        print ("{0} is the output catalog  ....\n".format(params.output2))
        print ("Creating ds9 check region file....\n")
        ds9kron(params.output2,params.regoutfile,params.scale)
    elif(run1 ==1):
        print ("{0} is the output catalog  ....\n".format("hot.cat"))
        print ("Creating ds9 check region file....\n")
        ds9kron("hot.cat",params.regoutfile,params.scale)
    elif(run2==1):
        print ("{0} is the output catalog  ....\n".format("cold.cat"))
        print ("Creating ds9 check region file....\n")
        ds9kron("cold.cat",params.regoutfile,params.scale)





#####
#####           creating mask

#########################################################################
########### crea la mascara a partir del archivo  ###################### 
#####################################################################



    if (create == 1):

        print ("Creating masks....\n")

        (NCol, NRow) = GetAxis(image)

#        print(output2,scale,SexArSort,NCol,NRow)

        if (run1 == 1 and run2 == 1):
            Total = CatArSort(params.output2,params.scale,params.SexArSort,NCol,NRow)
        elif(run1 ==1):
            Total = CatArSort("hot.cat",params.scale,params.SexArSort,NCol,NRow)
        elif(run2==1):
            Total = CatArSort("cold.cat",params.scale,params.SexArSort,NCol,NRow)

#        ParVar.Total = catfil.CatSort(ParVar)

##### segmentation mask

        MakeImage(params.maskfile, NCol, NRow)

        MakeMask(params.maskfile, params.SexArSort, params.scale,0,params.satfileout)  # offset set to 0
        MakeSatBox(params.maskfile, params.satfileout, Total + 1, NCol, NRow)

        print ("Running ds9 ....\n")
        runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} {} ".format(image,params.regoutfile,params.satfileout,params.maskfile)
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT

    else:

        if (run1 == 1 or run2 == 1):
            print ("Running ds9 ....\n")
            runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} ".format(image,params.regoutfile,params.satfileout)
            err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT
        else:
            print ("Ds9 can not run because sextractor was not used ")




    #lastmod

#########################################################################
########### anadir calculo de cielo para cada objeto del catalogo  ###################### 
#######################################################################




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
        checkflag = CheckFlag(flg[idx], flagsat)
        # check if object is inside of a saturaded box region indicated by
        # user in ds9
#        regflag = CheckSatReg(xx[idx], yy[idx], Rkron[idx], theta[idx], e[idx],regfile)
        regflag = CheckSatReg2(xx[idx], yy[idx],regfile)



        if (checkflag == False) and (regflag == False):

            print ("Creating ellipse mask for object {}  ".format(n[idx]))
            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rkron[idx], theta[idx], e[
                idx], sxsmin[idx], sxsmax[idx], sysmin[idx], sysmax[idx])

        elif(checkflag == True or regflag == True):

            print ("Skipping object {}, one or more pixels are saturated \n".format(n[idx]))

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






def CheckSatReg(x,y,R,theta,ell,filein):
   "Check if object is inside of saturated region. returns True if at least one pixel is inside"
## check if object is inside of
## saturaded region as indicated by ds9 box region
## returns True if object center is in saturaded region

   q = (1 - ell)

   bim = q * R

   theta = theta * np.pi /180  ## Rads!!!

   flag = False
#   fileflag =1



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

                   distcen = np.sqrt(dx**2 + dy**2)

##

                   dxlb = xlo - x
                   dylb = ylo - y

                   distleftbot = np.sqrt(dxlb**2 + dylb**2)

                   dxlt = xlo - x
                   dylt = yhi - y

                   distleftop = np.sqrt(dxlt**2 + dylt**2)

                   dxrb = xhi - x
                   dyrb = ylo - y

                   distrightbot = np.sqrt(dxrb**2 + dyrb**2)

                   dxrt = xhi - x
                   dyrt = yhi - y

                   distrightop = np.sqrt(dxrt**2 + dyrt**2)

#####
                   landa=np.arctan2( dy,dx )

                   if landa < 0:
                       landa=landa + 2 * np.pi

                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xell =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yell =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

                   dxe = xell - x
                   dye = yell - y
#####
                   landa=np.arctan2( dylb,dxlb )

                   if landa < 0:
                       landa=landa + 2 * np.pi

                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xellb =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yellb =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

                   dxelb = xellb - x
                   dyelb = yellb - y
#####
                   landa=np.arctan2( dylt,dxlt )

                   if landa < 0:
                       landa=landa + 2 * np.pi


                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xellt =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yellt =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

                   dxelt = xellt - x
                   dyelt = yellt - y
#####
                   landa=np.arctan2( dyrb,dxrb )

                   if landa < 0:
                       landa=landa + 2 * np.pi


                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xellrb =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yellrb =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

                   dxerb = xellrb - x
                   dyerb = yellrb - y
#####
                   landa=np.arctan2( dyrt,dxrt )

                   if landa < 0:
                       landa=landa + 2 * np.pi


                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xellrt =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yellrt =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

                   dxert = xellrt - x
                   dyert = yellrt - y
#####

####################

                   distell = np.sqrt(dxe**2 + dye**2)

                   distellb = np.sqrt(dxelb**2 + dyelb**2)
                   distellt = np.sqrt(dxelt**2 + dyelt**2)
                   distellrb = np.sqrt(dxerb**2 + dyerb**2)
                   distellrt = np.sqrt(dxert**2 + dyert**2)


                   if ( (xell > xlo and xell < xhi) and (yell > ylo and yell < yhi)  ):
                       flag=True
                       break


                   if ( (xellb > xlo and xellb < xhi) and (yellb > ylo and yellb < yhi)  ):
                       flag=True
                       break
                   if ( (xellt > xlo and xellt < xhi) and (yellt > ylo and yellt < yhi)  ):
                       flag=True
                       break
                   if ( (xellrb > xlo and xellrb < xhi) and (yellrb > ylo and yellrb < yhi)  ):
                       flag=True
                       break
                   if ( (xellrt > xlo and xellrt < xhi) and (yellrt > ylo and yellrt < yhi)  ):
                       flag=True
                       break


                   if (distcen < distell):
                       flag=True
                       break
                   if (distleftbot < distellb):
                       flag=True
                       break
                   if (distleftop < distellt):
                       flag=True
                       break
                   if (distrightbot < distellrb):
                       flag=True
                       break
                   if (distrightop < distellrt):
                       flag=True
                       break

   return flag




#end of program
if __name__ == '__main__':
    main()
