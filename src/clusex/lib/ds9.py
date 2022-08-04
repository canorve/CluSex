#! /usr/bin/env python

import numpy as np


from clusex.lib.check import CheckFlag 

import argparse

def sex2ds9():
    """creates a ds9 reg file from sextractor catalog"""

    parser = argparse.ArgumentParser(description="sex2ds9: Creates a ds9 reg file from sextractor catalog")

    # required arguments
    parser.add_argument("SexCatalog",help="sextractor catalog")

    #optional arguments
    parser.add_argument("-s","--scale", type=float, help="factor that multiplies the radius of the catalog objects. Default = 1",default=1)
    
    parser.add_argument("-off","--offset", type=float, help="factor that it is added to the scale times radius of the catalog objects. Default = 0",default=0)


    parser.add_argument("-o","--outreg", type=str, help="name of the output ds9 reg file ",default='ds9.reg')


    args = parser.parse_args()

    sexcatalog = args.SexCatalog
    scale = args.scale
    offset = args.offset
    regoutfile = args.outreg
 
    print("Creating ds9 reg file: ",regoutfile)

    ds9kron(sexcatalog,regoutfile,scale,offset)

    print("done.")
    


def ds9kron(sexfile,regfile,scale,offset):
    "Creates ds9 region file to check output catalog "



    f_out= open(regfile, "w")

    #    scale = 1
    #offset=0


    flagsat=4      ## flag value when object is saturated (or close to)


    #print OUT "image \n";


    line="image \n"
    f_out.write(line)

    count=0

    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Star,Flg=np.genfromtxt(sexfile,delimiter="",unpack=True)

    for idx, item in enumerate(N):

        Rkron = scale * Ai[idx] * Kr[idx] + offset

       

        if Rkron == 0:

            Rkron = 1

        q = (1 - E)
        bim = q * Rkron


        check=CheckFlag(Flg[idx],flagsat)  ## check if object doesn't has saturated regions


        if  (check == False) :


            line="ellipse({0},{1},{2},{3},{4}) # color=blue move=0 \n".format(X[idx],Y[idx],Rkron,bim[idx],Theta[idx])

            f_out.write(line)


            line2="point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={2} {3} {4} \n".format(X[idx],Y[idx],"{",np.int(N[idx]),"}")

            f_out.write(line2)

        else:
            count +=1

            line="ellipse({0},{1},{2},{3},{4}) # color=red move=0 \n".format(X[idx],Y[idx],Rkron,bim[idx],Theta[idx])

            #f_out.write(line)


            line2="point({0},{1}) # point=boxcircle color=red font=\"times 10 bold\" text={2} {3} {4} \n".format(X[idx],Y[idx],"{",np.int(N[idx]),"}")

            f_out.write(line2)




    print ("{} objects have  at least one saturated pixels  \n".format(count))


    #        f_out.write(line)



    f_out.close()




