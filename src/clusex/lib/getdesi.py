#!/usr/bin/env python3

import os
import numpy as np
from urllib.request import urlretrieve
from PIL import Image
from numpy import array, uint8
import matplotlib.pylab as pl
from pandas.io.parsers import read_csv

from PIL import UnidentifiedImageError
from urllib import request, error
from io import BytesIO
from os import system
import socket
import argparse
import sys

TIMEOUT = 20  # seconds

def desiListParser():
    """
    Function for argument parsing
    """

    parser = argparse.ArgumentParser(description="desilist: Download png images from DESI (with SDSS optional) given a file containing: name, ra, dec colums")

    parser.add_argument("file", help="file containing: name, ra, and dec columns")
    parser.add_argument("-s","--scale",type=float, default=0.262, 
            help="Desi pixel scale default=0.262 (value of DESI. Value of SDSS is 0.396). Changing this value will zoom in or zoom out the image")
    parser.add_argument("-w","--width",type=int, default=512 , 
            help="width of the image. Default = 512")
    parser.add_argument("-he","--height",type=int, default=512 , 
            help="height of the image. Default = 512")

    parser.add_argument("-sh","--show", action='store_true', help="shows the image ")
    parser.add_argument("-sdss","--SDSS", action='store_true', help="download image from Sloan digital Sky Survey")


    args = parser.parse_args()

    desiList(args.file, args.scale, args.width, args.height,  args.show, args.SDSS)


def desiList(file, scale=0.262, width=512, height=512, show=False, sdss=False):

    (direc,ext) = file.split(".")

    # creating directories
    if not os.path.exists(direc):
        os.makedirs(direc)

    name, ra, dec = np.genfromtxt(file, delimiter="", unpack=True,comments="#",dtype="U")

    name = name.astype(str)
    ra = ra.astype(float)
    dec = dec.astype(float)

    for idx, item in enumerate(name):
        namfil=direc + "/"
        namefil=name[idx] + ".jpg"
        namfil= namfil + namefil
        simg(ra=ra[idx], dec=dec[idx], show=show, savename=namfil, scale=scale, width=width, height=height, sdss=sdss)


def simg(ra=37.228,dec=0.37,
    scale=0.396, width=512, height=512,
    savename=None, show=True, sdss=False):
    '''
    Get jpg image, view it,
    and remove temporary jpg file.
    '''
    jpg = DownImage(ra, dec, scale=scale,width=width,height=height,
    savename=savename, sdss=sdss)

    if show == True:
        jpg.show()



class DownImage(object):
    '''
    Class for an DESI/SDSS png image.

      RA, DEC - J2000, degrees
      SCALE - plate scale in arsec per pixel
      WIDTH, HEIGHT - size of image in pixels
      SAVENAME - if none provided, defaults to 'object.jpg'
      DR - integer value for SDSS data release.
    '''

    def __init__(self, ra, dec,
                 scale=0.3515625, width=512, height=512,
                 savename=None, sdss=False):
        self.ra = ra
        self.dec = dec
        self.scale = scale
        self.width = width
        self.height = height
        self.DR=14
        if savename==None:
            savename = 'object.jpg'
        self.savename = savename
        if sdss:
            self.downloadSDSS()
        else:
            self.downloadDESI()
    def __del__(self):
        pass
        #system('rm '+self.savename)

    def downloadDESI(self):
        url = 'http://legacysurvey.org/viewer/jpeg-cutout?'
        url += 'ra=%0.5f&dec=%0.5f&'%(self.ra, self.dec)
        url += 'width=%i&height=%i&'%(self.width, self.height)
        url += 'layer=dr8&'
        url += 'pixscale=%0.5f&'%self.scale
        url += 'bands=grz'
        print(url)

        try:
            # Descarga con timeout; lanza HTTPError (4xx/5xx) o URLError (DNS, conexión, etc.)
            with request.urlopen(url, timeout=TIMEOUT) as r:
                data = r.read()

            # Verifica que realmente es una imagen abrible por PIL
            img = Image.open(BytesIO(data))
            img.load()  # fuerza la decodificación para detectar corrupciones temprano
            img.save(self.savename, format="jpeg")
            print(f"Download image: {self.savename}")

        except (error.HTTPError, error.URLError, socket.timeout) as e:
            # Error de red o del servidor → reporta y termina
            print(f"ERROR: fail to download ({self.savename}, ra={self.ra:.5f}, dec={self.dec:.5f}). Cause: {e}")
            sys.exit(1)
        except (UnidentifiedImageError, OSError) as e:
            # Respuesta no válida o imagen corrupta → reporta y termina
            print(f"ERROR: respond is not a valid image ({self.savename}, ra={self.ra:.5f}, dec={self.dec:.5f}). Causa: {e}")
            sys.exit(1)



    def downloadSDSS(self):
        url = 'http://skyservice.pha.jhu.edu/dr%i/ImgCutout/getjpeg.aspx?'%self.DR
        url += 'ra=%0.5f&dec=%0.5f&'%(self.ra, self.dec)
        url += 'scale=%0.5f&'%self.scale
        url += 'width=%i&height=%i'%(self.width, self.height)
        print(url)

        try:
            # Descarga con timeout; lanza HTTPError (4xx/5xx) o URLError (DNS, conexión, etc.)
            with request.urlopen(url, timeout=TIMEOUT) as r:
                data = r.read()

            # Verifica que realmente es una imagen abrible por PIL
            img = Image.open(BytesIO(data))
            img.load()  # fuerza la decodificación para detectar corrupciones temprano
            img.save(self.savename, format="jpeg")
            print(f"Download image: {self.savename}")

        except (error.HTTPError, error.URLError, socket.timeout) as e:
            # Error de red o del servidor → reporta y termina
            print(f"ERROR: fail to download ({self.savename}, ra={self.ra:.5f}, dec={self.dec:.5f}). Cause: {e}")
            sys.exit(1)
        except (UnidentifiedImageError, OSError) as e:
            # Respuesta no válida o imagen corrupta → reporta y termina
            print(f"ERROR: respond is not a valid image ({self.savename}, ra={self.ra:.5f}, dec={self.dec:.5f}). Causa: {e}")
            sys.exit(1)



    def show(self):
        pl.imshow(self.data)
        pl.ion()
        pl.show()


if __name__ == '__main__':
    desiListParser()




