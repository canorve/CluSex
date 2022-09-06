#! /usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
#from astropy.utils.data import get_pkg_data_filename

def CheckFlag(val,check):
   "Check for flag contained in $val, returns 1 if found "

   flag = False
   mod = 1
   max=128


   while (mod != 0):

       res = int(val/max)

       if (max == check and res == 1 ):

           flag=True

       mod = val % max

       val = mod
       max = max/2


   return flag


def CheckSatReg2(x,y,filein):
   "Check if object is inside of saturated region. returns True if at least one pixel is inside"
## check if object is inside of
## saturaded region as indicated by ds9 box region
## returns True if object center is in saturaded region

   flag = False

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

                   if ( (x > xlo and x < xhi) and (y > ylo and y < yhi) ):
                       flag=True
                       break


   return flag


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


def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow

def GetWCS(Image):
    # k Check
    "Get World Cordinate System info"
    hdu = fits.open(Image)[0]
    wcs = WCS(hdu.header)

    return wcs 

def GetCounts(Image):
    # k Check
    "Get Counts from Image"

    hdu = fits.open(Image)[0]
    counts = hdu.data.ravel()

    return counts




