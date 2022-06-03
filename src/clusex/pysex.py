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


from clusex.lib.check  import GetAxis

from clusex.lib.ds9 import ds9kron


from clusex.lib.mask import MakeMask 
from clusex.lib.mask import MakeImage
from clusex.lib.mask import MakeSatBox 
from clusex.lib.mask import CatArSort


from clusex.lib.radcor import RadMod


from clusex.lib.satbox import SatBox 

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



    runsex(image,params)


##############################################################
##### llamado a la funcion radcor, rad correccion################### 
#########################################################

    print("correcting radius of galaxies")


    #incorporate this two parameters to input file
    tol = 0.5 #tolerance for radius differences between the two catalogs (proportion) 
    red = 0.3  # reduction factor

    RadMod("hot.cat","cold.cat","hot2.cat",tol)
    RadMod("cold.cat","hot.cat","cold2.cat",tol,red)


    os.rename("hot2.cat","hot.cat")
    os.rename("cold2.cat","cold.cat")





##############################################################
##### une los dos archivos hot.cat y cold.cat ####### 
#########################################################


    if (params.run1 == 1 and params.run2 == 1):
        print ("joining hot.cat and cold.cat catalogs ....\n");
        joinsexcat("hot.cat","cold.cat",params.output,params.scale,params.scale2)

    else:
        print("Unable to join catalogs because sextractor was not used \n")



##############################################################
##### busca y selecciona las estrellas saturadas ####### 
##### modifica las banderas para los objetos cercanos a regiones 
###### de saturacion ####
#########################################################
    

    SatBox(image,params)
    




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



    if (params.create == 1):

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

        if (params.run1 == 1 or params.run2 == 1):
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




#end of program
if __name__ == '__main__':
    main()