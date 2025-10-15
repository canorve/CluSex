#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 16:13:00 2025

@author: samnay
"""

import os, shutil
import pandas as pd
from urllib.request import urlretrieve
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import argparse


def downloadDesi(sexcat="hotcold.cat", output_dir="desi_images", scale=1, offset=0, desi_pixscale=0.262, image_plate=0.68, max_mag=18, starclass=0.6):

    output_dir = "desi_images"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        print(f"Deleted: {output_dir}")
    else:
        print("Directory not found. Creating one")



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

        url = f'http://legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&width={width}&height={height}&layer=dr8&desi_pixscale={desi_pixscale}&bands=grz'


        temp_jpg = os.path.join(output_dir, f"temp_{idx:04d}.jpg")
        urlretrieve(url, temp_jpg)

        img = Image.open(temp_jpg)

        output_png = os.path.join(output_dir, f"imagen_{number:04d}.png")
        img.save(output_png, format="PNG")

        
        os.remove(temp_jpg)

        print(f"Download image: {output_png}")

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


if __name__ == '__main__':
    getDesiParser()


