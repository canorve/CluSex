#! /usr/bin/env python

import numpy as np
import argparse
import os
from astropy.io import fits
import subprocess as sp

from clusex.lib.join import joinsexcat 
from clusex.lib.join import putFlagSat 


from clusex.lib.init import printWelcome

from clusex.lib.check import GetAxis
from clusex.lib.make import CatArSort
from clusex.lib.make import MakeImage 
from clusex.lib.make import MakeMask
from clusex.lib.make import MakeSatBox 
from clusex.lib.make import MakeStamps
from clusex.lib.make import MakeObjImg


from clusex.lib.ds9 import ds9kron

from clusex.lib.make import EraseObjectMask2
from clusex.lib.make import EraseObjectMask


from clusex.lib.sky import SkyCal 



#console scripts



def joincat():
    """ joins two sextractor catalogs"""

    printWelcome()

    parser = argparse.ArgumentParser(description="Joincat: quickly combines two Sextractor catalogs")

    # required arguments
    parser.add_argument("FirstCatalog",help="First sextractor catalog")
    parser.add_argument("SecondCatalog",help="Second sextractor catalog")


    parser.add_argument("-s","--joinscale", type=float, help="factor that multiplies the radius of the First catalog objects. Objects outside of this catalog will be added. Default = 1",default=1)

    parser.add_argument("-o","--output", type=str, help="output catalog ",default='out.cat')

    parser.add_argument("-sf","--SatFile", type=str, help="Saturation DS9 reg file")

    parser.add_argument("-i","--include",action="store_true", help="Include all the galaxies from the second catalog that were not in the first catalog ")


    parser.add_argument("-r","--reduction", type=float, help="Reduction factor for objetc radius for -i option only",default=0.1)

    parser.add_argument("-m","--minrad", type=float, help="minrad for object radius for -i option only ",default=5)



    args = parser.parse_args()

    firstsex = args.FirstCatalog
    secondsex = args.SecondCatalog
    joinscale = args.joinscale
    output = args.output
    satfile = args.SatFile
    incFlag = args.include
    minrad= args.minrad
    red = args.reduction


    ##
    line="joining {} with {} using a scale of {}".format(firstsex,secondsex,joinscale)
    print(line)


    joinsexcat(firstsex,secondsex,output,joinscale,incFlag,red=red,minrad=minrad)


    if satfile:

        print("recomputing saturation flags for catalog")
        putFlagSat(output,"temp.cat",satfile)

        os.rename("temp.cat",output)

    else:

        print("no saturation DS9 reg file")

    line="joincat finished. output file {} created".format(output)
    print(line)


def makemask():
    """makes a mask from sextractor catalog"""


    printWelcome()

    parser = argparse.ArgumentParser(description="MakeMask: Creates mask from sextractor catalog")

    # required arguments
    parser.add_argument("SexCatalog",help="sextractor catalog")
    parser.add_argument("Image",help="Fits image of the objects")

    #optional arguments
    parser.add_argument("-s","--scale", type=float, help="factor that multiplies the radius of the catalog objects. Default = 1",default=1)
    
    parser.add_argument("-off","--offset", type=float, help="factor that it is added to the scale times radius of the catalog objects. Default = 0",default=0)

    parser.add_argument("-o","--outmask", type=str, help="name of the output mask ",default='mask.fits')

    parser.add_argument("-sf","--SatFile", type=str, help="Saturation DS9 reg file",default='ds9sat.reg')

    #options without arguments

    parser.add_argument("-n","--nodisplay",action="store_true", help="doesn't run DS9")


    args = parser.parse_args()

    sexcatalog = args.SexCatalog
    scale = args.scale
    offset = args.offset
    output = args.outmask
    satfile = args.SatFile
    image = args.Image
    flagds9 = args.nodisplay


    sexarsort="sexarea.cat"
    regoutfile="mask.reg"
    

    print ("Creating mask....\n")

    (NCol, NRow) = GetAxis(image)

    #check if exits satfile

    if not(os.path.exists(satfile)):
        with open(satfile, 'x') as f:
            f.write('box(1,1,0,0,0)') #to avoid empty file
            f.close()

    Total = CatArSort(sexcatalog,scale,offset,sexarsort,NCol,NRow)
   
    ##### segmentation mask

    MakeImage(output, NCol, NRow)

    MakeMask(output, sexarsort, scale, offset, satfile)  # offset set to 0
    MakeSatBox(output, satfile, Total + 1, NCol, NRow)

    #calling ds9kron to create ds9 reg objects
    ds9kron(sexcatalog,regoutfile,scale,offset)


    if not(flagds9): 
        print ("Running ds9 ...\n")
        runcmd="ds9 -tile column -cmap grey -invert -log -zmax -regions shape box {} -regions {} -regions {} {} ".format(image,regoutfile,satfile,output)
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  

    print('done') 


def remellmask():
    """Removes an ellipse from mask image"""


    printWelcome()

    parser = argparse.ArgumentParser(description="remellmask: Removes an ellipse from mask image")

    # required arguments
    parser.add_argument("mask",help="Fits mask image of the objects")

    parser.add_argument("number", type=int, help="ellipse number to be removed (or changed check option -f)")
    
    #optional arguments
    parser.add_argument("-f","--fill", type=float, help="Value number to be filled within the ellipse mask. 0 = object removed. Default = 0",default=0)

    parser.add_argument("-o","--outmask", type=str, help="name of the new mask ",default='objmask.fits')


    #options without arguments



    args = parser.parse_args()

    mask = args.mask
    number = args.number
    fill = args.fill
    output = args.outmask


    sexarsort="sexarea.cat"
    regoutfile="mask.reg"
    

    line="removing ellipse object:{} from mask: {}".format(number,mask)

    print (line)

    EraseObjectMask(mask,output,number,fill)

    print('new mask image:',output) 

    print('done')


def makestamps():
    """Make image stamps for each object"""


    printWelcome()

    parser = argparse.ArgumentParser(description="Make image stamps for every object in the catalog")

    # required arguments
    parser.add_argument("image",help="Fits image of the objects")
    parser.add_argument("catalog",help="sextractor catalog of the objects")
    parser.add_argument("mask",help="Fits image mask (created with makemask) ")

    #optional arguments
    parser.add_argument("-sr","--stretch", type=float, help="stretch factor to enlarge the stamps",default=5)

    parser.add_argument("-so","--skyoff", type=float, help="sky offset to be added the value of the mean sky ",default=1)


    parser.add_argument("-dp","--dpi", type=int, help="dots per inch resolution for image stamps",default=100)

    parser.add_argument("-cm","--cmap", type=str, help="color map",default='gray_r')


    parser.add_argument("-s","--scale", type=float, help="factor that multiplies the radius of the catalog objects. Default = 1",default=1)
    
    parser.add_argument("-off","--offset", type=float, help="factor that it is added to the scale times radius of the catalog objects. Default = 0",default=0)


    parser.add_argument("-br","--bright", type=float, help="brightness of the image. Default = 33",default=33)

    parser.add_argument("-co","--contrast", type=float, help="contrast of the image. Default = 0.98",default=0.98)


    parser.add_argument("-gc","--galclass", type=float, help="galaxy/star sextractor classification limit. Sextractor Classification  1 = Star. 0 = galaxy",default=1)

    args = parser.parse_args()

    image= args.image
    catalog= args.catalog
    mask = args.mask
    stretch = args.stretch
    skyoff = args.skyoff
    dpi= args.dpi
    cmap = args.cmap

    scale = args.scale
    offset  = args.offset
    bright = args.bright
    contrast = args.contrast
    galclass= args.galclass



    line="Creating image stamps for every object of the catalog"

    print (line)

    MakeStamps(image, catalog, mask, stretch, skyoff, dpi, cmap, 
                scale, offset, bright, contrast, galclass)


    print('done')




def sex2ds9():
    """creates a ds9 reg file from sextractor catalog"""

    printWelcome()


    parser = argparse.ArgumentParser(description="sex2ds9: Creates a DS9 reg file from sextractor catalog")

    # required arguments
    parser.add_argument("SexCatalog",help="sextractor catalog")

    #optional arguments
    parser.add_argument("-s","--scale", type=float, help="factor that multiplies the radius of the catalog objects. Default = 1",default=1)
    
    parser.add_argument("-off","--offset", type=float, help="factor that it is added to the scale times radius of the catalog objects. Default = 0",default=0)


    parser.add_argument("-o","--outreg", type=str, help="name of the output DS9 reg file ",default='ds9.reg')


    args = parser.parse_args()

    sexcatalog = args.SexCatalog
    scale = args.scale
    offset = args.offset
    regoutfile = args.outreg
 
    print("Creating DS9 reg file: ",regoutfile)

    ds9kron(sexcatalog,regoutfile,scale,offset)

    print("done.")
    

def makeobjimg():
    """Make objects image """

    printWelcome()

    parser = argparse.ArgumentParser(description="Make object image")

    # required arguments
    parser.add_argument("image",help="Fits image of the objects")
    parser.add_argument("mask",help="Fits image mask (created with makemask) ")

    parser.add_argument("-o","--outimage", type=str, help="name of the output obj image file ",default='objects.fits')

    args = parser.parse_args()

    image= args.image
    maskimage = args.mask
    outimage = args.outimage
 


    hdu = fits.open(image)
    img = hdu[0].data

    hdumask = fits.open(maskimage)
    mask = hdumask[0].data
    hdumask.close()


    objimage = MakeObjImg(img, mask)


    hdu[0].data = objimage 
    hdu.writeto(outimage, overwrite=True)
    hdu.close()


    print('new objects image:',outimage) 

    print("done")


def compsky():
    """Computes the sky for every object of a sextractor catalog"""

    printWelcome()

    parser = argparse.ArgumentParser(description="skysex: computes the sky for every object of a sextractor catalog")

    # required arguments
    parser.add_argument("SexCatalog",help="sextractor catalog")
    parser.add_argument("Image",help="Fits image of the objects")

    parser.add_argument("MaskFile",help="Fits mask image. File created with makemask")


    parser.add_argument("-s","--scaleRadius", type=float, help="factor that multiplies the radius of the catalog objects. For grad sky it is the Inital radius. For rand sky is the minimum radius of the box around main object. Default = 1",default=1)
   
    parser.add_argument("-w","--width", type=int, help="width of the ring for the grad method. ",default=20)

    parser.add_argument("-b","--box", type=int, help="size of the box for the random method. ",default=20)


    parser.add_argument("-nb","--numBox", type=int, help="number of boxes for the random method. ",default=20)

    parser.add_argument("-sm","--scaleRadMax", type=float, help="factor that multiplies the radius of the catalog objects. For rand sky is the maximum radius of the box around main object. Default = 20",default=20)
 

    parser.add_argument("-m","--method", type=int, help="method used for compute sky. Grad sky =1, rand sky = 2. Default = 1  ",default=1)


    parser.add_argument("-o","--outcat", type=str, help="name of the output catalog ",default='out.cat')




    args = parser.parse_args()

    sexcatalog = args.SexCatalog
    image = args.Image
    mask = args.MaskFile
    scale = args.scaleRadius
    width = args.width

    box = args.box
    numbox = args.numBox

    scaleMax = args.scaleRadMax

    method = args.method

    output = args.outcat


    ##end input




    ### reading sex catalog:
    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(sexcatalog,delimiter="",unpack=True)

    AR         = 1 - E
    RKron      = scale * Ai * Kr
    RKronMax      = scaleMax * Ai * Kr

    Theta = Theta + 90 #it must be the same as GALFIT 

    maskron = RKron <= 0
    RKron[maskron] = 1

    maskar = AR <= 0.005

    AR[maskar] = 0.005


    

    #######################################
    ############## SKY ####################
    #######################################

    f_out = open(output, "w")

    #  gradient sky method:
    if method == 1:


        # computing sky with the gradient method
        print("Computing sky with the gradient method")

        ImageFile = image
        MaskFile = mask

        tempMask = "tempmask.fits"
        #width = params.skywidth

        hdumask = fits.open(MaskFile)
        maskimg = hdumask[0].data
        hdumask.close()

        hdu = fits.open(ImageFile)
        datimg = hdu[0].data
        hdu.close()


       
        for idx, item in enumerate(N):

            xx = X[idx]
            yy = Y[idx] 


            thetadeg = Theta[idx]
            q = AR[idx]
            Rinit = RKron[idx]


            print("computing sky for object ",N[idx])

            line="using Rinit = {:.2f} width = {}".format(Rinit,width)
            print(line)

            line="using thetadeg = {:.2f} q = {:.2f}".format(thetadeg,q)
            print(line)
        
            line="using x = {} y  = {}".format(xx,yy)
            print(line)

            #EraseObjectMask(MaskFile,tempMask,N[idx])
            tempmask=EraseObjectMask2(maskimg,N[idx]) #routine more efficient 


            #mean,std, median,rad = SkyCal().GetEllipSky(ImageFile,tempMask,xx,yy,
            #                                            thetadeg,q,Rinit,width,
            #                                            "ring.fits","ringmask.fits")
            mean,std, median,rad = SkyCal().GetEllipSky(datimg,tempmask,xx,yy,
                                                        thetadeg,q,Rinit,width,
                                                        "ring.fits","ringmask.fits")


            line="Total sky:  mean = {:.2f}; std={:.2f}; median = {:.2f} ".format(mean,std,median)
            print(line)

            #saving for output
            Bkgd[idx] = mean  
            #galpar.gradskystd = std
            #galpar.gradskymed = median


            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11:.2f} {12:.2f} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


            f_out.write(line)



    #  random sky method:
    elif method == 2:

        # computing sky  using random boxes across the image
        print("computing sky with the random box method")

        ImageFile = image
        MaskFile = mask

        hdu = fits.open(ImageFile)
        datimg = hdu[0].data
        hdu.close()


        hdumask = fits.open(MaskFile)
        maskimg = hdumask[0].data
        hdumask.close()





        for idx, item in enumerate(N):

            xx = X[idx]
            yy = Y[idx] 


            thetadeg = Theta[idx]
            q = AR[idx]
            Rinit = RKron[idx]
            Rmax = RKronMax[idx]


            ###

            #box = params.skybox
            #num = params.skynum
            #box = args.box
            #numbox = args.numBox

            #scaleMax = args.scaleRadMax

            print("computing sky for object ",N[idx])
            line="using Rad = {:.2f}, box size= {}, number of boxes = {}".format(Rinit,box,numbox)
            print(line)


           # mean,std, median = SkyCal().RandBox(ImageFile,MaskFile,xx,yy,
           #                                     thetadeg,q,Rinit,box,numbox,Rmax)
            mean,std, median = SkyCal().RandBox(datimg,maskimg,xx,yy,
                                                thetadeg,q,Rinit,box,numbox,Rmax)

            line="Total sky:  mean = {:.2f}; std = {:.2f}; median = {:.2f}".format(mean,std,median)
            print(line)

            #saving for output
            Bkgd[idx] = mean  
            #galpar.randskystd = std
            #galpar.randskymed = median
            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11:.2f} {12:.2f} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


            f_out.write(line)



    else: 

        # computing sky  using random boxes across the image
        print("method not found. use '-m 1' for grad sky and '-m 2' for random sky")




    f_out.close()


    print("computation of sky has finished")

    #######################################
    ############## SKY End ################
    #######################################



