#! /usr/bin/env python3

# Ellipse masking demo on your sample image.
# Edit the parameters (xc, yc, a, b, theta_deg) to your ellipse.

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse
import os
import shutil
import pandas as pd




def maskImage(sexcat, paths, scale, offset, desi_pixscale, 
                image_plate, out_dir):


    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        print(f"Deleted: {out_dir}")
    else:
        print("Directory not found. Creating one")

    os.makedirs(out_dir, exist_ok=True)


    df2 = pd.read_csv(sexcat,sep=" ", header=None)

    #scale = 1.5
    #offset = 3



    df1 = pd.read_csv(paths, sep=r"\s+", header=None, names=["path"], engine="python")

    # ID desde el nombre del archivo
    df1["ID"] = df1["path"].apply(extrae_id)


    df2.rename(columns={0: "ID"}, inplace=True)

    # 3) Unir (left join) usando el ID. Mantiene el número de filas de df1.
    merged = pd.merge(df1, df2, on="ID", how="left", sort=False)

    # 4) Reordenar columnas: primero path, clase, ID y luego las columnas originales del segundo archivo
    otras_cols = [c for c in merged.columns if c not in ["path",  "ID"]]
    merged = merged[["path", "ID"] + otras_cols]



    #columnas X_new y Y_new para RA y DEC
    img_paths =  merged.iloc[:,0]
    numbers = merged.iloc[:,1]
    X = merged.iloc[:,4]
    Y = merged.iloc[:,5]
    krons = merged.iloc[:,7]
    Ais = merged.iloc[:,10]
    Es = merged.iloc[:,11]
    Thetas = merged.iloc[:,12]

    Thetas = Thetas*(-1) # for conversion from Sextractor

    Rkrons = scale*Ais*krons + offset

    Rkrons = Rkrons * image_plate / desi_pixscale #modifying scale from original to DESI

    Rkronsb = Rkrons*(1-Es)

    ################

    # Load example image placed at /mnt/data by the system message
    #img_path = "desi_images/imagen_0006.png"


    for idx, (img_path, number, a, b, theta_deg) in enumerate(zip(img_paths, numbers, Rkrons, Rkronsb, Thetas)):

        outimage = applyMask(img_path, number , a, b, theta_deg, out_dir)



def applyMask(img_path, number , a, b, theta_deg, outdir):
   
    # Load example image placed at /mnt/data by the system message
    #img_path = "desi_images/imagen_0006.png"
    im = Image.open(img_path).convert("RGB")
    arr = np.array(im)

    h, w, _ = arr.shape
    
    # ---- Ellipse parameters (in pixels) ---
    xc, yc = w/2, h/2

    theta = np.deg2rad(theta_deg)
    
    # Create coordinate grid
    y, x = np.mgrid[0:h, 0:w]
    x_shift = x - xc
    y_shift = y - yc
   
    # Rotate coordinates (clockwise rotation of the ellipse frame)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    xr =  x_shift * cos_t + y_shift * sin_t
    yr = -x_shift * sin_t + y_shift * cos_t
    
    # Ellipse implicit equation: (xr/a)^2 + (yr/b)^2 <= 1
    mask_inside = (xr/a)**2 + (yr/b)**2 <= 1.0
    
    # Apply mask: keep inside ellipse, set outside to black
    clean = arr.copy()
    clean[~mask_inside] = 0
    
    # Save result
    #out_path = "imagen_0006_masked.png"

    out_path = os.path.join(outdir, f"imagen_{number:04d}.png")

    Image.fromarray(clean).save(out_path)
    
    return(out_path)






def extrae_id(path: str) -> int:
    """
    Extrae el ID numérico del nombre de archivo tipo imagen_0004.png -> 4.
    Lanza ValueError si no encuentra un entero.
    """
    m = re.search(r'(\d+)\.png$', path.strip())
    if not m:
        raise ValueError(f"No se pudo extraer ID de: {path}")
    return int(m.group(1).lstrip("0") or "0")




if __name__ == '__main__':
    maskArgParser()


