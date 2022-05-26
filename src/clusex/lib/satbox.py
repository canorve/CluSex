



def Ds9SatBox (image,satfileout,output,satscale,satoffset,satlevel,minsatsize,satq):
    "Creates a file for ds9 which selects bad saturated regions"

    scaleflag=1
    offsetflag=1
    regfileflag=1
    magflag=1
    clasflag=1

    flagsat=4      ## flag value when object is saturated (or close to)
    check=0
    regflag = 0    ## flag for saturaded regions


### read image

    hdu = fits.open(image)
    imgdat = hdu[0].data
    hdu.close()

###
    (imaxx, imaxy) = GetAxis(image)

    #corrected for python array
    imaxx= imaxx -1
    imaxy= imaxy -1


    f_out = open(satfileout, "w")

    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(output,delimiter="",unpack=True)


    line="image \n"
    f_out.write(line)


    for idx, item in enumerate(N):


        bi=Ai[idx]*(1-E[idx])

        Theta[idx] = Theta[idx] * np.pi /180  #rads!!!

        Rkronx =  2 * Ai[idx] * Kr[idx] * np.cos(Theta[idx])
        Rkrony =  2 * bi * Kr[idx] * np.sin(Theta[idx])

        if Rkronx <= minsatsize:
            Rkronx = minsatsize

        if Rkrony <= minsatsize:
            Rkrony = minsatsize

        Rkronx=np.int(np.round(Rkronx))
        Rkrony=np.int(np.round(Rkrony))

        check=CheckFlag(Flg[idx],flagsat)  # check if object has saturated region

        if (check):

            # -1 correction for python array
            xmin= X[idx] - Rkronx / 2 - 1
            xmax= X[idx] + Rkronx / 2  -1
            ymin= Y[idx] - Rkrony / 2  -1
            ymax= Y[idx] + Rkrony / 2  -1

            xmin=np.int(np.round(xmin))
            xmax=np.int(np.round(xmax))
            ymin=np.int(np.round(ymin))
            ymax=np.int(np.round(ymax))

            # check if xmin, xmax, ymin, ymax have values outsides of image limits
            if xmin < 0:
                xmin=0
            if xmax > imaxx:
                xmax=imaxx

            if ymin < 0:
                ymin=0
            if ymax > imaxy:
                ymax=imaxy

            chunkimg = imgdat[ymin:ymax,xmin:xmax]


            satm= (chunkimg > satlevel) + (chunkimg < (-1)*satlevel)
            satind= np.where(satm)

            y,x=satind

            errmsg="Can't found sat pixels in [{}:{},{}:{}]. Try to increase MinSatSize value \n".format(xmin+1,xmax+1,ymin+1,ymax+1)
            assert x.size != 0, errmsg

            satymin = y.min() +  ymin
            satxmin = x.min() +  xmin
            satymax = y.max() +  ymin
            satxmax = x.max() +  xmin

            sidex= x.max() - x.min()
            sidey= y.max() - y.min()

            xcent= satxmin + sidex/2 + 1
            ycent= satymin + sidey/2 + 1

            xcent=np.int(np.round(xcent))
            ycent=np.int(np.round(ycent))


            # increasing sides of saturated regions
            sidex=sidex*satscale + satoffset
            sidey=sidey*satscale + satoffset



            if sidex >= sidey:
                div = sidey/sidex
                if div < satq:
                    sidey=sidex*satq
                    sidey=np.int(np.round(sidey))
            else:
                div = sidex/sidey
                if div < satq:
                    sidex=sidey*satq
                    sidex=np.int(np.round(sidex))





#            line="box({0},{1},{2},{3},0) # color=red move=0 \n".format(X[idx],Y[idx],sidex,sidey)
            line="box({0},{1},{2},{3},0) # color=red move=0 \n".format(xcent,ycent,sidex,sidey)

            f_out.write(line)

#            line2="point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={{ {2} }} \n".format(X[idx],Y[idx],np.int(N[idx]))
            line2="point({0},{1}) # color=red point=boxcircle font=\"times 10 bold\" text={{ {2} }} \n".format(xcent,ycent,np.int(N[idx]))

            f_out.write(line2)


    f_out.close()


# this need an update
def putflagsat(sexfile,sexfile2,regfile):
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


        check=CheckFlag(Flg[idx],flagsat)  ## check if object doesn't has saturated regions
#        regflag=CheckSatReg(X[idx],Y[idx],Rkron,Theta[idx],E[idx],regfile)    ## check if object is inside of a saturaded box region indicated by user in ds9
        regflag=CheckSatReg2(X[idx],Y[idx],regfile)    ## check if object is inside of a saturaded box region indicated by user in ds9


        if  (check == False ) and ( regflag == True) :

            Flg[idx] = Flg[idx] + 4


            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

            f_out.write(line)


        else:

            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])


            f_out.write(line)



    f_out.close()



