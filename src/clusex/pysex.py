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
    ##### Reading argument parser #######
    #########################################


    parser = argparse.ArgumentParser(description="CluSex: combines two Sextractor catalogs among other stuff")

    # required arguments
    parser.add_argument("ConfigFile",help="CluSex configuration file ")
    #parser.add_argument("image",help="FITS image file of the galaxy cluster ")  # incorporated   


    args = parser.parse_args()

    confile = args.ConfigFile 
    #image = args.image

    #########################################
    ##### initial parameters ########
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



    params=Params()

    #########################################
    ##### Reading paramter file #######
    #########################################

    readcon(params,confile)

    ##############################################################
    ##### Reading and writing sextractor files ########
    #########################################################

    wsex(params)

    ##############################################################
    ##### Running Sextractor ########
    #########################################################



    runsex(params)


    #################################################################
    ##### Correcting radius for large low surface brightness galaxies ################### 
    #################################################################

    print("correcting radius of galaxies")



    if (params.run1 == 1 and params.run2 == 1):

        #incorporate this two parameters to input file
        #tol = 0.5 #tolerance for radius differences between the two catalogs (proportion) 
        #red = 0.3  # reduction factor

        if not(params.flagminrad):
            params.minrad = params.seeing


        RadMod("hot.cat","cold.cat","hot2.cat",tol=params.tol,
                red=params.red,minrad=params.minrad,scalecor=params.scalecor)
        RadMod("cold.cat","hot.cat","cold2.cat",tol=params.tol,
                red=params.red,minrad=params.minrad,scalecor=params.scalecor)


        os.rename("hot2.cat","hot.cat")
        os.rename("cold2.cat","cold.cat")

    else:

        print("Unable to correct radius because it needs to run sextractor twice. run = 1 run2 = 1 \n")

    #####################################################################
    ##### joins the two output sextractor files into one single file ############## 
    #####################################################################


    if (params.run1 == 1 and params.run2 == 1):
        print ("joining hot.cat and cold.cat catalogs ....\n");
        joinsexcat("hot.cat","cold.cat",params.output,params.joinscale)

    else:
        print("Unable to join catalogs because sextractor was not used \n")



    ###########################################################################
    ##### search for saturated stars and figure out the saturated area  ####### 
    ##### It also modify flags of objects near to saturared stars ################ 
    ##########################################################################
        

    SatBox(params)
    

    ##################################################################
    ###### creates the output files for DS9 ###################### 
    ##################################################################


    if (params.run1 == 1 and params.run2 == 1):
        print ("{0} is the output catalog  ....\n".format(params.output2))
        print ("Creating ds9 check region file....\n")
        ds9kron(params.output2,params.regoutfile,params.scale,params.offset)
    elif(params.run1 ==1):
        print ("{0} is the output catalog  ....\n".format("hot.cat"))
        print ("Creating ds9 check region file....\n")
        ds9kron("hot.cat",params.regoutfile,params.scale,params.offset)
    elif(params.run2==1):
        print ("{0} is the output catalog  ....\n".format("cold.cat"))
        print ("Creating ds9 check region file....\n")
        ds9kron("cold.cat",params.regoutfile,params.scale,params.offset)




    #####################################################################
    ########### Make masks from output ################################## 
    #####################################################################



    if (params.create == 1):

        print ("Creating masks....\n")

        (NCol, NRow) = GetAxis(params.image)

#        print(output2,scale,SexArSort,NCol,NRow)

        if (run1 == 1 and run2 == 1):
            Total = CatArSort(params.output2,params.scale,params.offset,params.SexArSort,NCol,NRow)
        elif(run1 ==1):
            Total = CatArSort("hot.cat",params.scale,params.offset,params.SexArSort,NCol,NRow)
        elif(run2==1):
            Total = CatArSort("cold.cat",params.scale,params.offset,params.SexArSort,NCol,NRow)


        ##### segmentation mask

        MakeImage(params.maskfile, NCol, NRow)

        MakeMask(params.maskfile, params.SexArSort, params.scale,0,params.satfileout)  # offset set to 0
        MakeSatBox(params.maskfile, params.satfileout, Total + 1, NCol, NRow)

        print ("Running ds9 ....\n")
        runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} {} ".format(params.image,params.regoutfile,params.satfileout,params.maskfile)
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT

    else:

        if (params.run1 == 1 or params.run2 == 1):
            print ("Running ds9 ....\n")
            runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} ".format(params.image,params.regoutfile,params.satfileout)
            err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT
        else:
            print ("Ds9 can not run because sextractor was not used ")





    #########################################################################
    ########### Computes sky for every object and produces a new catalog ####
    #######################################################################

    #lastmod



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
