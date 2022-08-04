
from ellipsect.lib.libs import *

from ellipsect import *
import argparse

def compsky():
    """Computes the sky for every object of a sextractor catalog"""


    parser = argparse.ArgumentParser(description="skysex: computes the sky for every object of a sextractor catalog")

    # required arguments
    parser.add_argument("SexCatalog",help="sextractor catalog")
    parser.add_argument("Image",help="Fits image of the objects")

    parser.add_argument("MaskFile",help="Fits mask image. File created with makemask")


    parser.add_argument("-s","--scaleRadius", type=float, help="factor that multiplies the radius of the catalog objects. For grad sky it is the Inital radius. For rand sky is the minimum radius of the box around main object. Default = 1",default=1)
   
    parser.add_argument("-w","--width", type=int, help="width of the ring for the grad method. ",default=20)

    parser.add_argument("-b","--box", type=int, help="size of the box for the random method. ",default=20)


    parser.add_argument("-nb","--numBox", type=int, help="number of boxes for the random method. ",default=20)

    parser.add_argument("-sm","--scaleRadMax", type=float, help="factor that multiplies the radius of the catalog objects. For rand sky is the maximum radius of the box around main object. Default = 10",default=10)
 

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



    #galpar.ang = angle same as GALFIT

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


       
        for idx, item in enumerate(N):

            xx = X[idx]
            yy = Y[idx] 


            thetadeg = Theta[idx]
            q = AR[idx]
            Rinit = RKron[idx]


            print("computing sky for object ",N[idx])

            line="using Rinit = {:.2f} width = {}".format(Rinit,width)
            print(line)

            line="using thetadeg = {:.2f} q = {}".format(thetadeg,q)
            print(line)
        
            line="using xx = {} yy  = {}".format(xx,yy)
            print(line)

            EraseObjectMask(MaskFile,tempMask,N[idx])


            mean,std, median,rad = SkyCal().GetEllipSky(ImageFile,MaskFile,xx,yy,
                                                        thetadeg,q,Rinit,width,
                                                        "ring.fits","ringmask.fits")

            line="Total sky:  mean = {:.2f}; std={:.2f}; median = {:.2f} ".format(mean,std,median)
            print(line)

            #saving for output
            Bkgd[idx] = mean  
            #galpar.gradskystd = std
            #galpar.gradskymed = median


            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


            f_out.write(line)



    #  random sky method:
    if params.flagrandboxsky:

        # computing sky  using random boxes across the image
        print("computing sky with the random box method")

        ImageFile = image
        MaskFile = mask

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


            mean,std, median = SkyCal().RandBox(ImageFile,MaskFile,xx,yy,
                                                thetadeg,q,Rinit,box,numbox,Rmax)
        #

            line="Total sky:  mean = {:.2f}; std = {:.2f}; median = {:.2f}".format(mean,std,median)
            print(line)

            #saving for output
            Bkgd[idx] = mean  
            #galpar.randskystd = std
            #galpar.randskymed = median
            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


            f_out.write(line)






    f_out.close()

    #######################################
    ############## SKY End ################
    #######################################

def EraseObjectMask(MaskFile,tempMask,obj):

    hdumask = fits.open(MaskFile)
    data = hdumask[0].data

    mask = data == obj

    data[mask] = 0 # removing object from mask
    
    hdumask[0].data = data

    hdumask.writeto(tempMask, overwrite=True)
    
    hdumask.close()
    


class SkyCal:
    "This class compute the sky using two methods: random boxes and sky gradient"

    def RandBox(self,ImageFile,MaskFile,xx,yy,thetadeg,q,Rinit,box,num,Rmax,outliers=True):
        "random box method to compute sky"

        self.xx = xx 
        self.yy = yy

        self.thetadeg=90 + thetadeg

        self.q = q
        self.Rinit = Rinit
       
        self.box = box 
        self.num = num

        self.outliers = outliers

        ###

        hdumask = fits.open(MaskFile)
        self.maskimg = hdumask[0].data
        hdumask.close()


        hdu=fits.open(ImageFile)
        self.img = hdu[0].data
        hdu.close()

        ####
        
        (self.nrow,self.ncol)=self.img.shape

        xmin,xmax,ymin,ymax,Rend = self.GetXYRBorder()

        if Rmax != 0:
            if Rmax <= Rend:
                Rend = Rmax
            else:
                print("input skyRadmax is greater than image size")
            print("using Rmax = {:0.2f} ".format(Rend))

        if Rinit > Rend:

            Rinit= Rend/3 #is this number ok?
            print("Rinit is greater than image size")
            print("using Rinit = {:0.2f} ".format(Rinit))
            print("change this value with --radinit ")



        #mean, std, median  = self.GetRandBoxSky( Rinit, Rend )
        mean, var, median, N  = self.GetRandBoxSky2( Rinit, Rend )


        meansky = np.mean(mean)
        medsky = np.median(median)
        Npix = N.sum()
        stdsky = np.sqrt(var.sum()/Npix )


        #rmstd = stats.sem(mean)

        return meansky, stdsky,  medsky


    def GetXYRBorder(self):
        "this subroutine get the coordinates of the border"

        q =  self.q

        theta = self.thetadeg * (np.pi / 180)  # rads!!

        thetax=np.sqrt((np.cos(theta))**2 + (q**2)*(np.sin(theta))**2 )
        thetay=np.sqrt((q**2)*(np.cos(theta))**2 + (np.sin(theta))**2 )


        if (self.thetadeg >-45 and self.thetadeg <= 45):

            xmax=self.ncol
            xmin =1
            R1 = (xmax - self.xx)/thetax
            R2 = (self.xx - xmin)/thetax

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

        elif (self.thetadeg >45 and self.thetadeg <= 135):
            ymax=self.nrow
            ymin =1
            R1 = (ymax - self.yy)/thetay
            R2 = (self.yy - ymin)/thetay

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

        elif (self.thetadeg >135 and self.thetadeg <= 225):
            xmax=1
            xmin =self.ncol

            R1 = (xmax - self.xx)/thetax
            R2 = (self.xx - xmin)/thetax

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

     
        elif (self.thetadeg >225 and self.thetadeg <= 315):
            ymax=1
            ymin =self.nrow
            R1 = (ymax - self.yy)/thetay
            R2 = (self.yy - ymin)/thetay

            if (np.abs(R1)>=np.abs(R2)):
               R = R1
            else:
               R = R2

        R = np.abs(R) #avoids negative numbers
        bim = q * R

        # getting size

        xmin = self.xx - np.sqrt((R**2) * (np.cos(theta))**2 +
                         (bim**2) * (np.sin(theta))**2)

        xmax = self.xx + np.sqrt((R**2) * (np.cos(theta))**2 +
                         (bim**2) * (np.sin(theta))**2)

        ymin = self.yy - np.sqrt((R**2) * (np.sin(theta))**2 +
                         (bim**2) * (np.cos(theta))**2)

        ymax = self.yy + np.sqrt((R**2) * (np.sin(theta))**2 +
                         (bim**2) * (np.cos(theta))**2)

        mask = xmin < 1
        if mask.any():
            if isinstance(xmin,np.ndarray):
                xmin[mask] = 1
            else:
                xmin = 1

        mask = xmax > self.ncol
        if mask.any():
            if isinstance(xmax,np.ndarray):
                xmax[mask] = self.ncol
            else:
                xmax = self.ncol

        mask = ymin < 1
        if mask.any():
            if isinstance(ymin,np.ndarray):
                ymin[mask] = 1
            else:
                ymin = 1

        mask = ymax > self.nrow
        if mask.any():
            if isinstance(ymax,np.ndarray):
                ymax[mask] = self.nrow
            else:
                ymax =self.nrow

        xmin=np.int32(np.round(xmin))
        ymin=np.int32(np.round(ymin))
        xmax=np.int32(np.round(xmax))
        ymax=np.int32(np.round(ymax))


        return (xmin,xmax,ymin,ymax,np.int32(R))

    def GetRandBoxSky(self, Rinit, Rmax):
        '''compute mean, std, and median from random selected boxes'''

        (nrow,ncol)=self.img.shape

        # it obtains corners of Rinit
        (xmino, xmaxo, ymino, ymaxo) = self.GetSize(self.xx, self.yy, Rinit, self.thetadeg, self.q, self.ncol, self.nrow) 

        # it obtains corners of Rmax
        (xminf, xmaxf, yminf, ymaxf) = self.GetSize(self.xx, self.yy, Rmax, self.thetadeg, self.q, self.ncol, self.nrow) 
        
        Value=1 #  value of counts  of the  main target for  the mask image  
        self.maskimg = self.MakeKron(self.maskimg, Value, self.xx, self.yy, Rinit, self.thetadeg, self.q, xminf, xmaxf, yminf, ymaxf) 

        ########

        sky = np.array([])
        skystd = np.array([])
        skymed = np.array([])



        coordinates = self.MakeCoord(xmino-self.box,xmaxo+self.box,ymino-self.box,ymaxo+self.box,xmaxf,ymaxf)


        cont=self.num # this is the number of boxes to use 

        for idx,item in enumerate(range(cont)):

            flatimg,xinit,yinit=self.GetRandomPatch(self.img,self.maskimg,self.box,coordinates)
            xfin = xinit + self.box - 1
            yfin = yinit + self.box - 1 


            flatimg.sort()

            boxcont=0

            while( not(flatimg.any()) and (boxcont < 10)):

                if (boxcont == 0):
                    print("Picking another box ")

                flatimg,xinit,yinit=self.GetRandomPatch(self.img,self.maskimg,self.box,coordinates)
                xfin = xinit + self.box - 1
                yfin = yinit + self.box - 1


                flatimg.sort()

                boxcont+=1 # avoid eternal loop

            if (boxcont == 10):
                print("max. iteration reached. I couldn't find a box")
     
            tot=len(flatimg)

            top=round(.8*tot)
            bot=round(.2*tot)
            
            #imgpatch=flatimg[bot:top]
            imgpatch=flatimg#[bot:top]

            linebox = "Box:{}   xinit/fin: {}-{}; yinit/fin: {}-{}  ".format(item+1,xinit,xfin,yinit,yfin)

            mean=np.mean(imgpatch)
            std=np.std(imgpatch)

            median=np.median(imgpatch)

            linemean = "sky  mean = {:.2f}; std = {:.2f}; median = {:.2f}".format(mean,std,median)
            print(linemean)

            sky=np.append(sky,mean)
            skystd=np.append(skystd,std)
            skymed=np.append(skymed,median)

        return sky,skystd,skymed


    def GetRandBoxSky2(self, Rinit, Rmax):
        '''compute mean, std, and median from random selected boxes'''

        (nrow,ncol)=self.img.shape

        # it obtains corners of Rinit
        (xmino, xmaxo, ymino, ymaxo) = self.GetSize(self.xx, self.yy, Rinit, self.thetadeg, self.q, self.ncol, self.nrow) 

        # it obtains corners of Rmax
        (xminf, xmaxf, yminf, ymaxf) = self.GetSize(self.xx, self.yy, Rmax, self.thetadeg, self.q, self.ncol, self.nrow) 
        
        Value=1 #  value of counts  of the  main target for  the mask image  
        self.maskimg = self.MakeKron(self.maskimg, Value, self.xx, self.yy, Rinit, self.thetadeg, self.q, xminf, xmaxf, yminf, ymaxf) 

        ########

        sky = np.array([])
        skystd = np.array([])
        skymed = np.array([])

        boximg = [] #list not numpy array

        N= np.array([])


        coordinates = self.MakeCoord(xmino-self.box,xmaxo+self.box,ymino-self.box,ymaxo+self.box,xmaxf,ymaxf)


        cont=self.num # this is the number of boxes to use 



        for idx,item in enumerate(range(cont)):

            flatimg,xinit,yinit=self.GetRandomPatch(self.img,self.maskimg,self.box,coordinates)
            xfin = xinit + self.box - 1
            yfin = yinit + self.box - 1 


            flatimg.sort()

            boxcont=0

            while( not(flatimg.any()) and (boxcont < 10)):

                if (boxcont == 0):
                    print("Picking another box ")

                flatimg,xinit,yinit=self.GetRandomPatch(self.img,self.maskimg,self.box,coordinates)
                xfin = xinit + self.box - 1
                yfin = yinit + self.box - 1


                flatimg.sort()

                boxcont+=1 # avoid eternal loop

            if (boxcont == 10):
                print("max. iteration reached. I couldn't find a box") 
                return 0,0,0 # It couldn't found any box; terminating
    

            linebox = "Box:{}   xinit/fin: {}-{}; yinit/fin: {}-{}  ".format(item+1,xinit,xfin,yinit,yfin)

            tot=len(flatimg)

            top=round(.8*tot)
            bot=round(.2*tot)
            
            if self.outliers:   # eliminate top 80% and bottom 20%
                imgpatch=flatimg[bot:top]
            else:
                imgpatch=flatimg

            boximg.append(imgpatch) #save for later
            N=np.append(N,imgpatch.size) #save for later

            mean=np.mean(imgpatch)
            std=np.std(imgpatch)

            median=np.median(imgpatch)

            linemean = "sky  mean = {:.2f}; std = {:.2f}; median = {:.2f}".format(mean,std,median)
            print(linemean)

            sky=np.append(sky,mean)
            #skystd=np.append(skystd,std)
            skymed=np.append(skymed,median)



        totmean = np.mean(sky)
        # to compute standard deviation:
        for idx, item in enumerate(range(cont)):

            bimg=boximg[idx]

            # not quite the var of the box, but it is needed to compute the TOTAL std:
            skyvar =  ((bimg-totmean)**2).sum() 

            skystd=np.append(skystd,skyvar)


        return sky,skystd,skymed,N

    def MakeCoord(self,xmino,xmaxo,ymino,ymaxo,xmaxf,ymaxf):
        ''' creates (x,y) coordinates between the inner and outer box'''

        if xmino < 1:
            xmino = 1

        if ymino < 1:
            ymino = 1

        if xmaxo > xmaxf:
            xmaxo = xmaxf

        if ymaxo > ymaxf:
            ymaxo = ymaxf




        coordinates = [(x,y) for x in np.arange(0,xmaxf) for y in np.arange(0,ymaxf) if not((x >= xmino and x <= xmaxo) and ( y >= ymino and y <= ymaxo))]

        return coordinates


    def MakeKron(self,imagemat, idn, x, y, R, theta, q, xmin, xmax, ymin, ymax):
        "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

        xmin = np.int(xmin)
        xmax = np.int(xmax)
        ymin = np.int(ymin)
        ymax = np.int(ymax)

        bim = q * R

        theta = theta * np.pi / 180  # Rads!!!

        ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

        dx = xpos - x
        dy = ypos - y

        landa = np.arctan2(dy, dx)

        mask = landa < 0
        if mask.any():
            landa[mask] = landa[mask] + 2 * np.pi

        landa = landa - theta

        angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

        xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
            np.sin(angle) * np.sin(theta)
        yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
            np.sin(angle) * np.cos(theta)

        dell = np.sqrt((xell - x)**2 + (yell - y)**2)
        dist = np.sqrt(dx**2 + dy**2)

        mask = dist < dell
        imagemat[ypos[mask], xpos[mask]] = idn

        return imagemat

    def GetSize(self,x, y, R, theta, q, ncol, nrow):
        '''this subroutine get the maximun
        and minimim pixels for Kron and sky ellipse'''

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
            if isinstance(xmin,np.ndarray):
                xmin[mask] = 1
            else:
                xmin = 1

        mask = xmax > ncol

        if mask.any():
            if isinstance(xmax,np.ndarray):
                xmax[mask] = ncol
            else:
                xmax = ncol

        mask = ymin < 1
        if mask.any():
            if isinstance(ymin,np.ndarray):
                ymin[mask] = 1
            else:
                ymin = 1

        mask = ymax > nrow
        if mask.any():
            if isinstance(ymax,np.ndarray):
                ymax[mask] = nrow
            else:
                ymax = nrow

        xmin=np.int32(np.round(xmin))
        ymin=np.int32(np.round(ymin))
        xmax=np.int32(np.round(xmax))
        ymax=np.int32(np.round(ymax))


        return (xmin, xmax, ymin, ymax)

    def GetRandomPatch(self,imagemat,mimg,box,coordinates):
        '''get a random box patch of the imagemat'''



        xinit,yinit=self.GetRandomCoord(coordinates)

        xfin = xinit + box  - 1 
        yfin = yinit + box  - 1

        imagebox=imagemat[yinit - 1:yfin, xinit - 1:xfin]
        maskbox=mimg[yinit - 1:yfin, xinit - 1:xfin]

        invboxpatch=np.logical_not(maskbox)

        return imagebox[invboxpatch],xinit,yinit



    def GetRandomCoord(self,coordinates):
        '''choose xinit, and yinit from coordinates'''
        
        ridx = np.random.randint(0,len(coordinates)-1)


        xinit = coordinates[ridx][0] 
        yinit = coordinates[ridx][1] 

        return xinit,yinit 

        ######


    def GetEllipSky(self, ImageFile, MaskFile, xx, yy, thetadeg, q, Rinit, width,namering,ringmask,outliers=True):
        "Gradient sky method"

        self.xx = xx 
        self.yy = yy

        self.thetadeg = 90 + thetadeg
        self.q = q
        self.e = (1 - self.q)
        self.Rinit = Rinit
        self.width = width 

        self.NumRings = 5  # number of rings per loop # read this from function?

        self.ringmask = ringmask 

        self.outliers = outliers

        ###

        hdumask = fits.open(MaskFile)
        self.maskimg = hdumask[0].data
        hdumask.close()


        hdu = fits.open(ImageFile)
        self.img = hdu[0].data
        hdu.close()

        ####
       
        (self.nrow,self.ncol) = self.img.shape


        xmin,xmax,ymin,ymax,Rkron = self.GetXYRBorder()

        self.R = Rkron

        if self.Rinit > Rkron: # avoid radius greater than the border

            self.Rinit= Rkron/3 #is this number ok?
            print("Rinit is greater than image size")
            print("using Rinit = {:0.2f} ".format(self.Rinit))
            print("change this value with --radinit ")





        (xmin, xmax, ymin, ymax) = self.GetSize(self.xx, self.yy, Rkron, self.thetadeg, self.q, self.ncol, self.nrow) # obtain corners of R

        theta = self.thetadeg * np.pi / 180  # Rads!!!

        patch = self.maskimg[ymin - 1:ymax, xmin - 1:xmax] # logical patch mask image

        self.invpatch=np.logical_not(patch)


        Rings=np.arange(self.Rinit,self.R,self.width) # anillos de tamaño width

        sky = np.array([])
        skymed = np.array([])
        skystd = np.array([])
        radius = np.array([])


        count = 0
        idx=0


        ########################################
        #creating ring masks file

        masksky = np.zeros(np.shape(self.img))

        val = 1 
        for ridx, ritem in enumerate(Rings):


            bring = (Rings[ridx] + self.width) * self.q
            aring = Rings[ridx] + self.width 

            #incresing number of points per rad 
            points = 2*np.pi * np.sqrt(0.5*(aring**2 + bring**2)) #Aprox. 
            points = 2*points #doubling the number of points
            points = int(round(points))
           
            alpha = np.linspace(0,2*np.pi,points)
            for tidx, item in enumerate(range(self.width)):

                bim=(Rings[ridx]+item)*self.q

                tempxell = self.xx + (Rings[ridx]+item) * np.cos(alpha) * np.cos(theta) - bim * \
                  np.sin(alpha) * np.sin(theta)

                tempyell = self.yy + (Rings[ridx]+item) * np.cos(alpha) * np.sin(theta) + bim * \
                  np.sin(alpha) * np.cos(theta)

                tempxell = tempxell.round().astype("int")
                tempyell = tempyell.round().astype("int")

                tempxell,tempyell = self.CorSize(tempxell,tempyell)
                masksky[tempyell - 1, tempxell - 1] = val
                
            val += 1 


        hdu[0].data=masksky
        hdu.writeto(self.ringmask,overwrite=True) 

        ########################################
        ########################################



        #computing sky in every ring
        for ind, item in enumerate(Rings):


            maskring,idx=self.GetRingMask(masksky[ymin - 1:ymax, xmin - 1:xmax],idx)


            flatimg=self.img[ymin - 1:ymax, xmin - 1:xmax][maskring].flatten()  
            flatimg.sort()

            tot=len(flatimg)

            top=round(.8*tot)
            bot=round(.2*tot)


            if self.outliers:   # eliminate top 80% and bottom 20%
                imgpatch=flatimg[bot:top]
            else:
                imgpatch=flatimg

            mean=np.mean(imgpatch)
            std=np.std(imgpatch)

            median=np.median(imgpatch)

            sky=np.append(sky,mean)
            skymed=np.append(skymed,median)
            skystd=np.append(skystd,std)
            radius=np.append(radius,Rings[idx] + self.width/2)

            print("Ring = {}; rad = {:.2f}; sky mean = {:.2f}; sky std = {:.2f}; median: {:.2f} ".format(idx+1,Rings[idx] + self.width/2,mean,std, median))


            # calcular gradiente
            if (count >= self.NumRings):

                # [1:-1] avoiding the first and last element for gradient 
                gradmask = np.gradient(sky[1:-1]) >= 0 
               
                count = 0
                tempidx=np.where(np.gradient(sky[1:-1]) >= 0) 
                
                if (sky[1:-1][gradmask].any()): 
                    
                    savidx=tempidx[0][0]
                    maskring,none =self.GetRingMask(masksky[ymin - 1:ymax, xmin - 1:xmax],savidx)


                    print("sky computed in ring {} ".format(savidx+2))
                    
                    print("Ring radius = {:.2f} marked in {} ".format(radius[1:-1][savidx],namering))
                    print("the counts value within ring represent the long axis") 
                    self.img[ymin - 1:ymax, xmin - 1:xmax][maskring] = radius[1:-1][savidx] 
                    break

            count += 1
            idx +=1

            if idx == (len(Rings)-1): 
                print("The edge of image has been reached. Sky can not be computed")
                return 0,0,0,0


        hdu[0].data=self.img
        hdu.writeto(namering,overwrite=True) 

        finmean,finmedian,finstd,finRad = sky[1:-1][gradmask],skymed[1:-1][gradmask],skystd[1:-1][gradmask],radius[1:-1][gradmask]


        return finmean[0],finstd[0],finmedian[0],finRad[0]


    def GetRingMask(self,masksky,idx):
        ''' obtains the ring selected by index idx'''

        ring= idx + 1

        maskring = masksky == ring 

        maskring=maskring*self.invpatch

        ringcont=0

        while( not(maskring.any()) and (ringcont < 10)):

            if (ringcont == 0):
                    print("Selecting next ring ")

            idx += 1
            ring= idx + 1

            maskring = masksky == ring 

            maskring=maskring*self.invpatch


            ringcont+=1 # avoid eternal loop

        if (ringcont == 10):
            print("max. iteration reached. I couldn't find a ring") 
            return 0,0 # It couldn't found any ring ending 
 




        return maskring,idx


    def CorSize(self,xell,yell):
        '''Correct size for image borders'''
        masksx = xell < 1
        masksy = yell < 1

        xell[masksx] = 1
        yell[masksy] = 1

        masksx = xell > self.ncol
        masksy = yell > self.nrow 

        xell[masksx] = self.ncol 
        yell[masksy] = self.nrow

        return xell,yell

    ##### End of sky class ##################
    #########################################
    #########################################


