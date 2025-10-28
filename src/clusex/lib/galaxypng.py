#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 16:13:00 2025

@author: samnay
"""

import sys
import os, shutil
import pandas as pd
from urllib.request import urlretrieve
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import argparse
from PIL import UnidentifiedImageError
from urllib import request, error
from io import BytesIO
import socket

TIMEOUT = 20  # seconds


def DesiFromSexCat(sexcat="hotcold.cat", output_dir="desi_images", scale=1, offset=0, desi_pixscale=0.262, image_plate=0.68, max_mag=18, starclass=0.6, sdss=False):

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        print(f"Deleted: {output_dir}")
    else:
        print("Directory not found. Creating one")

    DR = 14 # sdss data release

    os.makedirs(output_dir, exist_ok=True)

    data = pd.read_csv(sexcat,sep=" ", header=None)


    #columnas X_new y Y_new para RA y DEC
    numbers = data.iloc[:,0]
    ra_list = data.iloc[:,1]
    dec_list = data.iloc[:,2]
    X = data.iloc[:,3]
    Y = data.iloc[:,4]
    mags = data.iloc[:,5]
    krons = data.iloc[:,6]
    Ais = data.iloc[:,9]
    Es = data.iloc[:,10]
    Thetas = data.iloc[:,11]

    classtar= data.iloc[:,13]
    flags = data.iloc[:,14]

    maxflag = 4 # avoid saturation

    # Filter
    mask  = (mags < max_mag) & (classtar < starclass) & (flags < maxflag)

    numbers = numbers[mask]
    ra_list = ra_list[mask]
    dec_list = dec_list[mask]
    mags = mags[mask] 
    krons = krons[mask] 
    Ais = Ais[mask] 
    Es = Es[mask] 
    Thetas = Thetas[mask] 


    classtar = classtar[mask] 
    flags = flags[mask]
    X = X[mask]
    Y = Y[mask]


    Rkrons = scale*Ais*krons + offset



    (xmin, xmax, ymin, ymax) = GetSize(X, Y, Rkrons, Thetas, Es) #size of every object

    #size in pix
    sizex = xmax - xmin
    sizey = ymax - ymin

    # Ajustes de tamaño 
    #pixscale = 0.262 # "/pix DESI 
    #locos_plate = 0.68 #"/pix  LOCOS 

    sizex = sizex * image_plate / desi_pixscale
    sizey = sizey * image_plate / desi_pixscale

    sizex = sizex.round().astype(int)
    sizey = sizey.round().astype(int)


    Tot = len(numbers)


    print(f"Number of images to download: {Tot}")

    #Recorrer cada galaxia
    for idx, (number, ra, dec, mag, width, height, clas, flag) in enumerate(zip(numbers,ra_list, dec_list, mags, sizex, sizey, classtar, flags)):


        if sdss:
            url = 'http://skyservice.pha.jhu.edu/dr%i/ImgCutout/getjpeg.aspx?'%DR
            url += 'ra=%0.5f&dec=%0.5f&'%(ra, dec)
            url += 'scale=%0.5f&'%desi_pixscale
            url += 'width=%i&height=%i'%(width, height)
        else:

            url = f'http://legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&width={width}&height={height}&layer=dr8&desi_pixscale={desi_pixscale}&bands=grz'


        try:
            # Descarga con timeout; lanza HTTPError (4xx/5xx) o URLError (DNS, conexión, etc.)
            with request.urlopen(url, timeout=TIMEOUT) as r:
                data = r.read()

            # Verifica que realmente es una imagen abrible por PIL
            img = Image.open(BytesIO(data))
            img.load()  # fuerza la decodificación para detectar corrupciones temprano

            output_png = os.path.join(output_dir, f"imagen_{number:04d}.png")
            img.save(output_png, format="PNG")
            print(f"Download image: {output_png}")

        except (error.HTTPError, error.URLError, socket.timeout) as e:
            # Error de red o del servidor → reporta y termina
            print(f"ERROR: fail to download (idx={idx}, ra={ra:.5f}, dec={dec:.5f}). Cause: {e}")
            sys.exit(1)
        except (UnidentifiedImageError, OSError) as e:
            # Respuesta no válida o imagen corrupta → reporta y termina
            print(f"ERROR: respond is not a valid image (idx={idx}, ra={ra:.5f}, dec={dec:.5f}). Causa: {e}")
            sys.exit(1)




    print(f"\n¡Process completed! Images live at '{output_dir}'.")




def GetSize(x, y, R, theta, ell):
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

    mask = ymin < 1
    if mask.any():
        ymin[mask] = 1

    xmin = np.rint(xmin)
    ymin = np.rint(ymin)
    xmax = np.rint(xmax)
    ymax = np.rint(ymax)

    return (xmin.astype(int), xmax.astype(int), ymin.astype(int), ymax.astype(int))



