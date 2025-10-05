#! /usr/bin/env python

import numpy as np
import argparse
import os
import sys

from clusex.lib.check  import CheckKron

from clusex.lib.check import CheckFlag 
from clusex.lib.check import CheckSatReg2



def joinsexcat (maincat,secondcat,output,JoinScale,incFlag=False, red=0.1,minrad = 5):
    "merges two Sextractor catalogs"

    f_out = open(output, "w")

    KronScale2 = 1

    try:
    #maincat
        N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(maincat,delimiter="",unpack=True)
    except IOError:
        print("bad sex file or wrong number of columns ")
        sys.exit(1)
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

    try: 
        #second cat
        N2,Alpha2,Delta2,X2,Y2,Mg2,Kr2,Fluxr2,Isoa2,Ai2,E2,Theta2,Bkgd2,Idx2,Flg2=np.genfromtxt(secondcat,delimiter="",unpack=True)
    except:
        print("bad sextractor file or wrong number of columns ")
        sys.exit(1)


    AR2         = 1 - E2
    RKron2      = KronScale2 * Ai2 * Kr2


    maskar2 = AR2 <= 0.005
    AR2[maskar2]=0.005


    count=0
    count2=0

    flag1 = False
    flag2 = False
 
    distmax = 5 # dist max to compare with incFlag=True

    dist = 0

    if incFlag:
        print("Including all the galaxies from the second catalog that are not in the first catalog" )
    else:
        print("Including all the galaxies from the second catalog that are outside from ellipse objects from the first catalog" )

    for idx2, item2 in enumerate(N2):

        if incFlag:
            flagf = True
            for idx, item in enumerate(N):
                
                dx = X[idx]-X2[idx2]
                dy = Y[idx]-Y2[idx2]
                dist=np.sqrt( dx**2  + dy**2 )

                if dist < distmax:   
                    flagf = False 
                    break

            if flagf:

                Kr2[idx2] = red * Kr2[idx2]

                if Kr2[idx2] * Ai2[idx2] < minrad:    

                    Kr2[idx2] = minrad/Ai2[idx2]
 

                line="{0:.0f} {1} {2} {3} {4} {5} {6:.2f} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(NewN, Alpha2[idx2], Delta2[idx], X2[idx2], Y2[idx2], Mg2[idx2], Kr2[idx2], Fluxr2[idx2], Isoa2[idx2], Ai2[idx2], E2[idx2], Theta2[idx2], Bkgd2[idx2], Idx2[idx2], Flg2[idx2])

                f_out.write(line)

                NewN+=1
                count2+=1
            else:
                count+=1       


        else:

            flagf =False
            for idx, item in enumerate(N):


                flag1=CheckKron(X2[idx2],Y2[idx2],X[idx],Y[idx],RKron[idx],Theta[idx],AR[idx])
                #flag2=CheckKron(X[idx],Y[idx],X2[idx2],Y2[idx2],RKron2[idx2],Theta2[idx2],AR2[idx2])
                flagf=flag1 or flag2

                if flagf:   # boolean value
                    break

            if not flagf:

                line="{0:.0f} {1} {2} {3} {4} {5} {6:.2f} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(NewN, Alpha2[idx2], Delta2[idx2], X2[idx2], Y2[idx2], Mg2[idx2], Kr2[idx2], Fluxr2[idx2], Isoa2[idx2], Ai2[idx2], E2[idx2], Theta2[idx2], Bkgd2[idx2], Idx2[idx2], Flg2[idx2])

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


