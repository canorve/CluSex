#! /usr/bin/env python

#    CluSex is a set of routines that aids Sextractor 
#    to perform on images of cluster galaxies
#    Copyright (C) 2022  Christopher Añorve 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path

from pathlib import Path
import shutil
import argparse


from clusex.lib.init import Params
from clusex.lib.init import printWelcome

from clusex.lib.read import readcon 
from clusex.lib.writesex import wsex
from clusex.lib.writesex import runsex

from clusex.lib.join import joinsexcat 



from clusex.lib.ds9 import ds9kron


from clusex.lib.make import MakeMask 
from clusex.lib.make import MakeImage
from clusex.lib.make import MakeSatBox 
from clusex.lib.make import CatArSort


from clusex.lib.radcor import RadMod


from clusex.lib.satbox import SatBox 


def main():

    ######################################
    ##### Reading argument parser ########
    ######################################

    printWelcome()


    parser = argparse.ArgumentParser(description="CluSex: combines two Sextractor catalogs among other stuff")

    # required arguments
    parser.add_argument("ConfigFile",help="CluSex configuration file ")
    #parser.add_argument("image",help="FITS image file of the galaxy cluster ")  # incorporated   


    args = parser.parse_args()

    confile = args.ConfigFile 
    #image = args.image



    ######################################
    ##### initial parameters #############
    ######################################


    #copying files to actual folder

    filsex = Path(__file__).parent / "config/default.sex"
    filnnw = Path(__file__).parent / "config/default.nnw"
    filcon = Path(__file__).parent / "config/default.conv"
    filpar = Path(__file__).parent / "config/sex.param"

    to_des = Path('.')

    shutil.copy2(filsex, to_des)  
    shutil.copy2(filnnw, to_des)  
    shutil.copy2(filcon, to_des) 
    shutil.copy2(filpar, to_des)



    params=Params()

    #########################################
    ##### Reading paramter file #############
    #########################################

    readcon(params,confile)

    ####################################################
    ##### Reading and writing sextractor files #########
    ####################################################

    wsex(params)

    ##################################
    ##### Running Sextractor #########
    ##################################



    runsex(params)


    ######################################################################
    ### Correcting radius for large low-surface brightness galaxies  ##### 
    ######################################################################

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
    ##### joins the two output sextractor files into one single file ####
    #####################################################################


    if (params.run1 == 1 and params.run2 == 1):
        print ("joining hot.cat and cold.cat catalogs ....\n");
        joinsexcat("hot.cat","cold.cat",params.output,params.joinscale)

    else:
        print("Unable to join catalogs because sextractor was not used \n")



    ###########################################################################
    ##### search for saturated stars and figure out the saturated area  ####### 
    ##### It also modify flags of objects near to saturated stars ############# 
    ##########################################################################
        

    SatBox(params)
    

    ##################################################################
    ###### creates the output files for DS9 ###################### 
    ##################################################################

    if params.flagDs9 == 1:


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


        print ("Running ds9 ...\n")
        runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} ".format(params.image,params.regoutfile,params.satfileout)
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  


    else:

        print("CluSex won't run ds9 because DisplayDs9 is set to 0 ")


    #####################################################################
    ###### Makemask has been moved to an independient routine ###########
    #####################################################################


    print("CluSex has finished.")






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
