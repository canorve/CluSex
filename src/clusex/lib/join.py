#! /usr/bin/env python

import numpy as np
import argparse
import os

from clusex.lib.check  import CheckKron

from clusex.lib.check import CheckFlag 
from clusex.lib.check import CheckSatReg2

def joincat():
    """ joins two sextractor catalogs"""

    parser = argparse.ArgumentParser(description="Joincat: quickly combines two Sextractor catalogs")

    # required arguments
    parser.add_argument("FirstCatalog",help="First sextractor catalog")
    parser.add_argument("SecondCatalog",help="Second sextractor catalog")


    parser.add_argument("-s","--joinscale", type=float, help="factor that multiplies the radius of the First catalog objects. Objects outside of this catalog will be added. Default = 1",default=1)

    parser.add_argument("-o","--output", type=str, help="output catalog ",default='out.cat')

    parser.add_argument("-sf","--SatFile", type=str, help="Saturation DS9 reg file")

    args = parser.parse_args()

    firstsex = args.FirstCatalog
    secondsex = args.SecondCatalog
    joinscale = args.joinscale
    output = args.output
    satfile = args.SatFile


    ##
    line="joining {} with {} using a scale of {}".format(firstsex,secondsex,joinscale)
    print(line)


    joinsexcat(firstsex,secondsex,output,joinscale)


    if satfile:

        print("recomputing saturation flags for catalog")
        putFlagSat(output,"temp.cat",satfile)

        os.rename("temp.cat",output)

    else:

        print("no saturation DS9 reg file")

    line="joincat finished output file  {} created".format(output)
    print(line)




def joinsexcat (maincat,secondcat,output,JoinScale):
    "merges two Sextractor catalogs"

    f_out = open(output, "w")

    KronScale2 = 1

    #maincat
    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(maincat,delimiter="",unpack=True)

    AR         = 1 - E
    RKron      = JoinScale * Ai * Kr



    maskron = RKron <= 0
    RKron[maskron]=1

    maskar = AR <= 0.005

    AR[maskar]=0.005

    for idx, item in enumerate(N):

        line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])



        f_out.write(line)


    total =len(N)

    NewN=total + 1


    #second cat
    N2,Alpha2,Delta2,X2,Y2,Mg2,Kr2,Fluxr2,Isoa2,Ai2,E2,Theta2,Bkgd2,Idx2,Flg2=np.genfromtxt(secondcat,delimiter="",unpack=True)

    AR2         = 1 - E2
    RKron2      = KronScale2 * Ai2 * Kr2


    maskar2 = AR2 <= 0.005
    AR2[maskar2]=0.005


    count=0
    count2=0

    flag1 = False
    flag2 = False

    for idx2, item2 in enumerate(N2):

        flagf =False
        for idx, item in enumerate(N):


            flag1=CheckKron(X2[idx2],Y2[idx2],X[idx],Y[idx],RKron[idx],Theta[idx],AR[idx])
            #flag2=CheckKron(X[idx],Y[idx],X2[idx2],Y2[idx2],RKron2[idx2],Theta2[idx2],AR2[idx2])
            flagf=flag1 or flag2

            if flagf:   # boolean value
                break

        if not flagf:

            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(NewN, Alpha2[idx2], Delta2[idx], X2[idx2], Y2[idx2], Mg2[idx2], Kr2[idx2], Fluxr2[idx2], Isoa2[idx2], Ai2[idx2], E2[idx2], Theta2[idx2], Bkgd2[idx2], Idx2[idx2], Flg2[idx2])

            f_out.write(line)

            NewN+=1
            count2+=1
        else:
            count+=1       
    f_out.close()

    linout="{} objects from second run rejected ".format(count)
    print(linout)

    linout="{} objects were added from second run ".format(count2)
    print(linout)



def putFlagSat(sexfile,sexfile2,regfile):
    """Put flags on objects which are inside saturated regions"""


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


