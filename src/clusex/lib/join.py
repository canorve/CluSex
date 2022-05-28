#! /usr/bin/env python

import numpy as np



def joinsexcat (maincat,secondcat,output,KronScale,KronScale2):
    "merges two Sextractor catalogs"

    f_out = open(output, "w")

#maincat
    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(maincat,delimiter="",unpack=True)

    AR         = 1 - E
    RKron      = KronScale * Ai * Kr



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



    for idx2, item2 in enumerate(N2):

        flagf =False
        for idx, item in enumerate(N):


            flag1=CheckKron(X2[idx2],Y2[idx2],X[idx],Y[idx],RKron[idx],Theta[idx],AR[idx])
            flag2=CheckKron(X[idx],Y[idx],X2[idx2],Y2[idx2],RKron2[idx2],Theta2[idx2],AR2[idx2])
            flagf=flag1 or flag2

            if flagf:   # boolean value
                break

        if not flagf:

            line="{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(NewN, Alpha2[idx2], Delta2[idx], X2[idx2], Y2[idx2], Mg2[idx2], Kr2[idx2], Fluxr2[idx2], Isoa2[idx2], Ai2[idx2], E2[idx2], Theta2[idx2], Bkgd2[idx2], Idx2[idx2], Flg2[idx2])

            f_out.write(line)

            NewN+=1
        else:
            linout="object {} from cold run rejected ".format(np.int(N2[idx2]))
            print(linout)
    f_out.close()


def CheckKron (xpos,ypos,x,y,R,theta,q):
    "check if position is inside of the Kron Ellipse saturaded region returns True if object center is in Ellipse region"


    bim = q * R

    theta = theta * np.pi /180  ## Rads!!!


    flag =False


    dx = xpos - x
    dy = ypos - y


    landa = np.arctan2( dy,dx )

    if landa < 0:
        landa=landa + 2 * np.pi


    landa = landa - theta


    angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

    xell =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
    yell =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)

    dell = np.sqrt ( (xell - x)**2 + (yell - y)**2 )
    dist = np.sqrt ( dx**2 + dy**2 )


    if  dist < dell:
	      flag=True


    return flag




