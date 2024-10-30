#! /usr/bin/env python

import numpy as np

from astropy.io import fits

import os


class SatBox:

    def __init__(self, params: object ):
        """
        This routine creates a Ds9 box region file where it contains the saturated 
        regions of the image.
        """




        if (params.run1 == 1 and params.run2 == 1):

            print ("creating {0} for ds9 ....\n".format(params.satfileout))
            self.Ds9SatBox(params.image,params.satfileout,params.output,params.satscale,params.satoffset,params.satlevel,params.minsatsize,params.satmethod,params.satq) 

            # crea archivo  de salida de reg
            print ("recomputing flags on objects which are inside saturated regions  ....\n")
            self.putflagsat(params.output,params.output2,params.satfileout)


        elif(params.run1 ==1):

            print ("creating {0} for ds9 ....\n".format(params.satfileout))
            self.Ds9SatBox(params.image,params.satfileout,"hot.cat",params.satscale,params.satoffset,params.satlevel,params.minsatsize,params.satmethod,params.satq) 

            # crea archivo  de salida de reg
            print ("recomputing flags on objects which are inside saturated regions  ....\n")
            self.putflagsat("hot.cat","hot2.cat",params.satfileout)

            os.rename("hot2.cat","hot.cat")



        elif(params.run2 == 1):

            print ("creating {0} for ds9 ....\n".format(params.satfileout))
            self.Ds9SatBox(params.image,params.satfileout,"cold.cat",params.satscale,params.satoffset,params.satlevel,params.minsatsize,params.satmethod,params.satq) 

            # crea archivo  de salida de reg
            print ("recomputing flags on objects which are inside saturated regions  ....\n")
            self.putflagsat("cold.cat","cold2.cat",params.satfileout)

            os.rename("cold2.cat","cold.cat")






    def Ds9SatBox (self,image,satfileout,sexcat,satscale,satoffset,satlevel,minsatsize,method,satq=0.5):
        "Creates a file for ds9 which selects bad saturated regions"


        (xx,yy,sx,sy,N)=self.GetSatBox(image,satfileout,sexcat,satlevel,minsatsize,satq)

        ### obtains maximun size
        (imaxx, imaxy) = self.GetAxis(image)

        #corrected for python array
        #imaxx= imaxx - 1
        #imaxy= imaxy - 1



        # increasing size of sat regions 
        sx = sx * satscale + satoffset
        sy = sy * satscale + satoffset


        (xx,yy,sx,sy,N)=self.AbsorbSmallBox(xx,yy,sx,sy,N) #remove small regions inside big ones 
        (xx,yy,sx,sy,N)=self.ResizeBox(xx,yy,sx,sy,N) #readjust size for neighbors regions
        #(xx,yy,sx,sy,N)=self.RemoveDupBox(xx,yy,sx,sy,N) #remove this one?
        (xx,yy,sx,sy,N)=self.AbsorbSmallBox(xx,yy,sx,sy,N) # final touch 

        
        f_out = open(satfileout, "w")

        for idx, item in enumerate(N):


            sidex=sx[idx]
            sidey=sy[idx]
            xcent=xx[idx]
            ycent=yy[idx]

            # increasing sides of saturated regions

            #sidex = sidex * satscale + satoffset
            #sidey = sidey * satscale + satoffset


            divflag = False
            
            #divides in two: horizontal and vertical box if axis ratio < satq
            if sidex >= sidey: #horizontal
                div = sidey/sidex
                if div <= satq:
                    divflag = True
                    sidey2=sidex*(satq) #reduction factor
                    sidex2=sidey
            else: #vertical
                div = sidex/sidey
                if div <= satq:
                    divflag = True
                    sidey2=sidex
                    sidex2=sidey*(satq) #reduction factor

            xmin,xmax,ymin,ymax = self.BoxSide2Corners(xcent, ycent, sidex, sidey)
            xmin,xmax,ymin,ymax = self.CorrectCorners(xmin,xmax,ymin,ymax,imaxx,imaxy)
            xcent, ycent, sidex, sidey  = self.BoxCorners2Side(xmin,xmax,ymin,ymax)


            line="box({0},{1},{2},{3},0) # color=red move=1 \n".format(xcent,ycent,sidex,sidey)
            f_out.write(line)


            if divflag:


               # here we have 4 methods for saturated regions divided in one vertical and one horizontal 
               # select the best method for your convinience. However the best method so far 
               # has been the third one, so this is one for default

                #method = 4
                
                posxmax,posymax = self.GetMaxCor(image,xmin,xmax,ymin,ymax)
                posxmax2,posymax2 = self.GetMaxCor2(image,xmin,xmax,ymin,ymax) # for method 4

                xmin,xmax,ymin,ymax = self.BoxSide2Corners(xcent, ycent, sidex2, sidey2)
                xmin,xmax,ymin,ymax = self.CorrectCorners(xmin,xmax,ymin,ymax,imaxx,imaxy)
 
                #method 1
                if method == 1:
                    # 5a. test moving the center of the new box or
                   xmin,xmax,ymin,ymax = self.BoxCenterInPoint(posxmax,posymax,xmin,xmax,ymin,ymax)
                    ###
                    
                #method 2
                elif method ==2:
                    # 5b. test increazing size of box.
                    xmin,xmax,ymin,ymax = self.IncludePointInBox(posxmax,posymax,xmin,xmax,ymin,ymax)

                #method 3
                elif method ==3:
                    #5c. combining methods
                    xmin,xmax,ymin,ymax = self.IncludeBoxCenterInPoint(posxmax,posymax,xmin,xmax,ymin,ymax)

                    ### 

                #method 4
                elif method ==4:
                    #5c. combining methods and avoids to select the greatest negative pixexl value as a center

                    xmin,xmax,ymin,ymax = self.IncludeBoxCenterInPoint(posxmax2,posymax2,xmin,xmax,ymin,ymax)

 


                xmin,xmax,ymin,ymax = self.CorrectCorners(xmin,xmax,ymin,ymax,imaxx,imaxy)
                xcent, ycent, sidex2, sidey2  = self.BoxCorners2Side(xmin,xmax,ymin,ymax)

                ##########

                line2="point({0},{1}) # color=red point=boxcircle font=\"times 10 bold\" text={{ {2} }} \n".format(xcent,ycent,int(N[idx]))

                f_out.write(line2)


                line="box({0},{1},{2},{3},0) # color=red move=1 \n".format(xcent,ycent,sidex2,sidey2)
                f_out.write(line)


        f_out.close()



    def IncludePointInBox(self,xx,yy,xmin,xmax,ymin,ymax):

            pflag=self.IsPointInBox(xx,yy,xmin,xmax,ymin,ymax)

            if not(pflag):

                if (xx <= xmin):
                    xmin = xx - 4 

                if (yy <= ymin):
                    ymin = yy - 4

                if (xx >= xmax):
                    xmax = xx + 4

                if (yy >= ymax):
                    ymax = yy + 4


            return xmin,xmax,ymin,ymax

    def BoxCenterInPoint(self,xx,yy,xmin,xmax,ymin,ymax):



            pflag=self.IsPointInBox(xx,yy,xmin,xmax,ymin,ymax)

            if not(pflag):

                xc, yc, sx,sy  = self.BoxCorners2Side(xmin,xmax,ymin,ymax)

                xc2 = xx
                yc2 = yy


                (xmin,xmax,ymin,ymax)=self.BoxSide2Corners(xc2,yc2,sx,sy)


            return xmin,xmax,ymin,ymax


    def IncludeBoxCenterInPoint(self,xx,yy,xmin,xmax,ymin,ymax):
            '''combines the two methods'''

            pflag=self.IsPointInBox(xx,yy,xmin,xmax,ymin,ymax)

            if not(pflag):

                if (xx <= xmin):
                    xmin = xx - 4 

                if (yy <= ymin):
                    ymin = yy - 4

                if (xx >= xmax):
                    xmax = xx + 4

                if (yy >= ymax):
                    ymax = yy + 4

                xc, yc, sx,sy  = self.BoxCorners2Side(xmin,xmax,ymin,ymax)

                xc2 = xx
                yc2 = yy

                (xmin,xmax,ymin,ymax)=self.BoxSide2Corners(xc2,yc2,sx,sy)



            return xmin,xmax,ymin,ymax





    def GetMaxCor(self,image,xmin,xmax,ymin,ymax):
            '''Get coordinate (x,y) where the max value in counts'''
        
            hdu = fits.open(image)
            imgdat = hdu[0].data
            hdu.close()

            chunkimg = np.abs(imgdat[ymin-1:ymax-1,xmin-1:xmax-1])

            mask = chunkimg < 0
            chunkimg[mask]=chunkimg[mask]*-1

            p = np.where(chunkimg == np.amax(chunkimg))

            return  p[1][0] + xmin, p[0][0] + ymin
 
    def GetMaxCor2(self,image,xmin,xmax,ymin,ymax):
            '''Get coordinate (x,y) where the max value in counts. Avoids negative pixels'''
        
            hdu = fits.open(image)
            imgdat = hdu[0].data
            hdu.close()

            chunkimg = imgdat[ymin-1:ymax-1,xmin-1:xmax-1]

            #mask = chunkimg < 0
            #chunkimg[mask]=chunkimg[mask]*-1

            p = np.where(chunkimg == np.amax(chunkimg))

            return  p[1][0] + xmin, p[0][0] + ymin
 



    def ResizeBox(self,xx,yy,sx,sy,N):

        xc, yc, sidx, sidy = np.array([]), np.array([]), np.array([]), np.array([])
        nn=np.array([])


        for idx, item in enumerate(N):

            xx1,yy1,sx1,sy1,nn1 = xx[idx],yy[idx],sx[idx],sy[idx],N[idx]
            erflag = True

            (xmin,xmax,ymin,ymax)=self.BoxSide2Corners(xx1,yy1,sx1,sy1)


            for idx2, item2 in enumerate(N):

                xx2,yy2,sx2,sy2 = xx[idx2],yy[idx2],sx[idx2],sy[idx2] # remove this line after tests

                if idx!=idx2:

                    (xmin2,xmax2,ymin2,ymax2)=self.BoxSide2Corners(xx2,yy2,sx2,sy2)
                    
                    pflag=self.IsPointInBox(xx2,yy2,xmin,xmax,ymin,ymax)

                    if pflag:

                        if (xmin2 <= xmin):
                            xmin = xmin2 - 4 

                        if (ymin2 <= ymin):
                            ymin = ymin2 - 4

                        if (xmax2 >= xmax):
                            xmax = xmax2 + 4

                        if (ymax2 >= ymax):
                            ymax = ymax2 + 4


                        xx1, yy1, sx1,sy1  = self.BoxCorners2Side(xmin,xmax,ymin,ymax)
                        #update new size in old array 
                        xx[idx],yy[idx],sx[idx],sy[idx] =  xx1, yy1, sx1,sy1

            xc=np.append(xc,xx1)
            yc=np.append(yc,yy1)
            sidx=np.append(sidx,sx1)
            sidy=np.append(sidy,sy1)
            nn=np.append(nn,nn1)

        return(xc,yc,sidx,sidy,nn)

    def AbsorbSmallBox(self,xx,yy,sx,sy,N):
        "Removes small boxes inside big boxes"
        xc, yc, sidx, sidy = np.array([]), np.array([]), np.array([]), np.array([])
        nn=np.array([])


        for idx, item in enumerate(N):

            xx1,yy1,sx1,sy1,nn1 = xx[idx],yy[idx],sx[idx],sy[idx],N[idx]
            foflag = False 

            for idx2, item2 in enumerate(N):

                xx2,yy2,sx2,sy2 = xx[idx2],yy[idx2],sx[idx2],sy[idx2]

                if idx!=idx2:
                    
                    bflag=self.IsBoxInBox(xx1,yy1,sx1,sy1,xx2,yy2,sx2,sy2)


                    if bflag == True:
                        foflag = True
                        break


            if foflag == False:

                xc=np.append(xc,xx1)
                yc=np.append(yc,yy1)
                sidx=np.append(sidx,sx1)
                sidy=np.append(sidy,sy1)
                nn=np.append(nn,nn1)


        return(xc,yc,sidx,sidy,nn)

    def RemoveDupBox(self,xx,yy,sx,sy,N):
        "Removes duplicated boxes and small boxes inside big boxes"

        xc, yc, sidx, sidy = np.array([]), np.array([]), np.array([]), np.array([])
        nn=np.array([])


        for idx, item in enumerate(N):

            xx1,yy1,sx1,sy1,nn1 = xx[idx],yy[idx],sx[idx],sy[idx],N[idx]

            foflag = False 

            for idx2, item2 in enumerate(nn):

                xx2,yy2,sx2,sy2 = xc[idx2],yc[idx2],sidx[idx2],sidy[idx2]
                    
                bflag=self.AreTwoBoxSame(xx1,yy1,sx1,sy1,xx2,yy2,sx2,sy2)

                if (bflag == True):
                    foflag = True
                    break


            if foflag == False:

                xc=np.append(xc,xx1)
                yc=np.append(yc,yy1)
                sidx=np.append(sidx,sx1)
                sidy=np.append(sidy,sy1)
                nn=np.append(nn,nn1)


        return(xc,yc,sidx,sidy,nn)




    def BoxSide2Corners(self,X,Y,sidex,sidey):

             # -1 correction for python array
            #xmin= X - sidex / 2 - 1
            #xmax= X + sidex / 2  -1
            #ymin= Y - sidey / 2  -1
            #ymax= Y + sidey / 2  -1

            # -1 correction for python array
            xmin= X - sidex / 2 
            xmax= X + sidex / 2
            ymin= Y - sidey / 2
            ymax= Y + sidey / 2

            xmin=int(np.round(xmin))
            xmax=int(np.round(xmax))
            ymin=int(np.round(ymin))
            ymax=int(np.round(ymax))


            return (xmin,xmax,ymin,ymax)


    def CorrectCorners(self,xmin,xmax,ymin,ymax,imaxx,imaxy):


            # correct if xmin, xmax, ymin, ymax have values outsides of image limits

            if xmin < 1:
                xmin=1
            if xmax > imaxx:
                xmax=imaxx

            if ymin < 1:
                ymin=1
            if ymax > imaxy:
                ymax=imaxy

            return (xmin,xmax,ymin,ymax)



    def BoxCorners2Side(self,xmin,xmax,ymin,ymax):

            sidex= xmax - xmin
            sidey= ymax - ymin

            # increase +1 for conversion to DS9
            #xcent = xmin + sidex/2 + 1
            #ycent = ymin + sidey/2 + 1

            xcent = xmin + sidex/2 
            ycent = ymin + sidey/2


            xcent=int(np.round(xcent))
            ycent=int(np.round(ycent))
    

            return (xcent,ycent,sidex,sidey)


    def putflagsat(self,sexfile,sexfile2,regfile):
        "Put flags on objects which are inside saturated regions"


        f_out= open(sexfile2, "w")

        scale = 1
        offset=0


        flagsat=4      ## flag value when object is saturated (or close to)


        N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(sexfile,delimiter="",unpack=True)

        for idx, item in enumerate(N):

            Rkron = scale * Ai[idx] * Kr[idx] + offset

            if Rkron == 0:

                Rkron = 1

            q = (1 - E)
            bim = q * Rkron


            check=self.CheckFlag(Flg[idx],flagsat)  ## check if object doesn't has saturated regions
    #        regflag=CheckSatReg(X[idx],Y[idx],Rkron,Theta[idx],E[idx],regfile)    ## check if object is inside of a saturaded box region indicated by user in ds9
            regflag=self.CheckSatReg2(X[idx],Y[idx],regfile)    ## check if object is inside of a saturaded box region indicated by user in ds9


            if  (check == False ) and ( regflag == True) :

                Flg[idx] = Flg[idx] + 4


                line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

                f_out.write(line)


            else:

                line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


                f_out.write(line)



        f_out.close()

    def GetAxis(self,Image):
        # k Check
        "Get number of rows and columns from the image"

        hdu = fits.open(Image)
        ncol = hdu[0].header["NAXIS1"]
        nrow = hdu[0].header["NAXIS2"]
        hdu.close()
        return ncol, nrow



    def CheckFlag(self,val,check):
       "Check for flag contained in $val, returns True if found "

       flag = False
       mod = 1
       fmax=128


       while (mod != 0):


           res = int(val/fmax)

           if (fmax == check and res == 1 ):

               flag=True


           mod = val % fmax

           val = mod
           fmax = fmax/2


       return flag


    def CheckSatReg2(self,x,y,filein):
       "Check if object is inside of saturated region. returns True if at least one pixel is inside"
    ## check if object is inside of
    ## saturaded region as indicated by ds9 box region
    ## returns True if object center is in saturaded region

       flag = False

       with open(filein) as f_in:

           lines = (line.rstrip() for line in f_in) # All lines including the blank ones
           lines = (line.split('#', 1)[0] for line in lines) # remove comments
           lines = (line.rstrip() for line in lines)   # remove lines containing only comments
           lines = (line for line in lines if line) # Non-blank lines

           for line in lines:

               if (line != "image"):

                   (box,info)=line.split('(')

                   if(box == "box"):

                       (xpos,ypos,xlong,ylong,trash)=info.split(',')

                       xpos=float(xpos)
                       ypos=float(ypos)
                       xlong=float(xlong)
                       ylong=float(ylong)


                       xlo = xpos - xlong/2
                       xhi = xpos + xlong/2

                       ylo = ypos - ylong/2
                       yhi = ypos + ylong/2

                       if ( (x > xlo and x < xhi) and (y > ylo and y < yhi) ):
                           flag=True
                           break


       return flag


    def GetSatBox (self,image,satfileout,sexcat,satlevel,minsatsize,satq):
        "Obtains the saturated boxes from the sextractor catalog "


        flagsat=4      ## flag value when object is saturated (or close to)


        ### read image
        hdu = fits.open(image)
        imgdat = hdu[0].data
        hdu.close()

        ### obtains maximun size
        (imaxx, imaxy) = self.GetAxis(image)

        #corrected for python array
        #imaxx= imaxx -1
        #imaxy= imaxy -1


        f_out = open(satfileout, "w")

        N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(sexcat,delimiter="",unpack=True)


        line="image \n"
        f_out.write(line)


        xx, yy, sx, sy = np.array([]), np.array([]), np.array([]), np.array([])
        nn=np.array([])

        for idx, item in enumerate(N):


            check=self.CheckFlag(Flg[idx],flagsat)  # check if object has saturated region

            if (check):

                bi=Ai[idx]*(1-E[idx])

                Theta[idx] = Theta[idx] * np.pi /180  #rads!!!

                Rkronx =  2 * Ai[idx] * Kr[idx] * np.cos(Theta[idx])
                Rkrony =  2 * bi * Kr[idx] * np.sin(Theta[idx])

                if Rkronx <= minsatsize:
                    Rkronx = minsatsize

                if Rkrony <= minsatsize:
                    Rkrony = minsatsize

                Rkronx=int(np.round(Rkronx))
                Rkrony=int(np.round(Rkrony))


                xmin,xmax,ymin,ymax = self.BoxSide2Corners(X[idx],Y[idx], Rkronx, Rkrony)
                xmin,xmax,ymin,ymax = self.CorrectCorners(xmin,xmax,ymin,ymax,imaxx,imaxy) #correct corners
                #chunkimg = imgdat[ymin:ymax,xmin:xmax]
                chunkimg = imgdat[ymin-1:ymax-1,xmin-1:xmax-1]


                satm= (chunkimg > satlevel) + (chunkimg < (-1)*satlevel)
                satind= np.where(satm)

                y,x=satind

                errmsg="Can't found sat pixels in [{}:{},{}:{}]. Try to increase MinSatSize value \n".format(xmin+1,xmax+1,ymin+1,ymax+1)
                assert x.size != 0, errmsg
                
                #return to original coordinates
                #satymin = y.min() + ymin
                #satxmin = x.min() + xmin
                #satymax = y.max() + ymin
                #satxmax = x.max() + xmin

                #return to original coordinates
                satymin = y.min() + ymin + 1
                satxmin = x.min() + xmin + 1
                satymax = y.max() + ymin + 1
                satxmax = x.max() + xmin + 1

               
                xcent, ycent, sidex, sidey  = self.BoxCorners2Side(satxmin,satxmax,satymin,satymax)


                xx=np.append(xx,xcent)
                yy=np.append(yy,ycent)
                sx=np.append(sx,sidex)
                sy=np.append(sy,sidey)
                nn=np.append(nn,N[idx])



        f_out.close()


        return (xx,yy,sx,sy,nn)
    

    def IsPointInBox(self,xx,yy,xmin,xmax,ymin,ymax):

        flag=False

        if ((xmin < xx < xmax ) and ( ymin < yy < ymax)):
            flag = True


        return flag


    def IsBoxInBox(self,xx1,yy1,sx1,sy1,xx2,yy2,sx2,sy2):
        "Check if box is inside another box"
        flag=False


        (xmin,xmax,ymin,ymax)=self.BoxSide2Corners(xx1,yy1,sx1,sy1)
        (xmin2,xmax2,ymin2,ymax2)=self.BoxSide2Corners(xx2,yy2,sx2,sy2)


    
        if ((xmin2 < xmin) and (xmax < xmax2) ):
            if ((ymin2 < ymin) and (ymax < ymax2) ):
                flag = True
 
        return flag


    def AreTwoBoxSame(self,xx1,yy1,sx1,sy1,xx2,yy2,sx2,sy2):
        "Check if two boxes are the same"

        flag=False


        if ((xx1  == xx2 ) and (yy1  == yy2) and (sx1  == sx2 ) and (sy1  == sy2)): 
                flag = True

        return flag





