#!/usr/bin/env python3

import re
import argparse
import pandas as pd

def extrae_id(path: str) -> int:
    """
    Extract the ID from namefile:  imagen_0004.png -> 4.
    returns ValueError if it is not found.
    """
    m = re.search(r'(\d+)\.png$', path.strip())
    if not m:
        raise ValueError(f"ID not found: {path}")
    return int(m.group(1).lstrip("0") or "0")



def joinclass(archivo1, archivo2, salida):
    # 1) read the first file 
    df1 = pd.read_csv(archivo1, sep=",", engine="python")


    # ID desde el nombre del archivo
    # ID read from the name of image file 
    df1["ID"] = df1["path"].apply(extrae_id)

    # 2) read the second file 
    df2 = pd.read_csv(archivo2, sep=r"\s+", header=None)

    # make sure the first columns is the ID 
    df2.rename(columns={0: "ID"}, inplace=True)

    # 3) join using ID. It keeps the rows of df1 
    merged = pd.merge(df1, df2, on="ID", how="left", sort=False)

    # 4) sort columns: first path, class, ID and other columns 
    otras_cols = [c for c in merged.columns if c not in ["path", 'ID' ,'class']]
    merged = merged[["path", 'ID' ,'class'] + otras_cols]

    # 5) Save output as CSV file 
    sep = "," if salida.endswith((".csv", ".txt")) else ","
    merged.to_csv(salida, index=False, sep=sep)

    # brief verification messages 
    print(f"rows file1: {len(df1)}")
    print(f"row output file: {len(merged)}")
    faltantes = merged[merged.filter(regex=r"^\d+$").columns].isna().all(axis=1).sum() if len(df2.columns) > 1 else 0
    if faltantes:
        print(f"Advertencia: {faltantes} filas sin coincidencia en archivo2 (valores NaN).")

#if __name__ == "__main__":

#    joinparser()
