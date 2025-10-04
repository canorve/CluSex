#! /usr/bin/env python

import numpy as np




class RadMod:

    def __init__(self, sexcat1: str, sexcat2: str, newcat: str, tol=1, red=1, 
                minrad=10,scalecor=1):
        """
        This routine modifies wrong large estimated radius for Sextractor catalogs. 
        It does so by comparing the radius of two catalogs of the same image. It 
        kept the smallest radius of the two catalogs.
        The objects radius in sexcat1 will be changed for the ones if sexcat2 if the
        tolerance is achieved. newcat will be the modified catalog
        """

        #first cat
        N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(sexcat1,delimiter="",unpack=True)


        #second cat
        N2,Alpha2,Delta2,X2,Y2,Mg2,Kr2,Fluxr2,Isoa2,Ai2,E2,Theta2,Bkgd2,Idx2,Flg2=np.genfromtxt(sexcat2,delimiter="",unpack=True)


        dmax= 3 # max distance to match the same object in different catalogs 

        count = 0
        count2 = 0

        for idx, item in enumerate(N):

            foundflag = False

            dx = X2 - X[idx]
            dy = Y2 - Y[idx]
            dist = np.sqrt(dx**2 + dy**2)

            distmin = dist.min()
            idx2 = dist.argmin()

            if distmin <= dmax:

                foundflag = True
                #cflag1 = self.CheckFlag(Flg[idx],4) #check for saturated flag 
                #cflag2 = self.CheckFlag(Flg2[idx2],4) #check for saturated flag 

                cflag1 = False 
                cflag2 = False 


                #the comparison is based in Fluxr parameter

                #rad = Fluxr[idx]
                #rad2 = Fluxr2[idx2]

                # or  you can use kr * ai

                rad = Kr[idx] * Ai[idx]
                rad2 =  Kr2[idx2] *  Ai2[idx2]

                kr1 = Kr[idx]
                kr2 = Kr2[idx2]

                ai1 = Ai[idx]
                ai2 = Ai2[idx2]


                if rad2 < .1: 
                    rad2 = .1 #avoids division by zero


                den = rad2 

                num = rad - rad2

                comp = num/den

                if ((cflag1  and cflag2) == False): 

                    if comp > tol:

                        count +=1 

                        # reduction factor included
                        #Kr[idx],Fluxr[idx],Isoa[idx],Ai[idx] = red * Kr2[idx2],Fluxr2[idx2],Isoa2[idx2],Ai2[idx2]
                        Kr[idx],Fluxr[idx],Isoa[idx],Ai[idx] =  Kr2[idx2],Fluxr2[idx2],Isoa2[idx2],Ai2[idx2]

                        #Ai[idx] =  Ai2[idx2]
                        #Kr[idx] =  Kr2[idx2]


                        #this is to avoid galaxies with very low radius after interchanging
                        #radius. However we don't want to increase it to a larger 
                        #radius than the previous one:

                        #kr2 is the smallest radii of the two 
                        

                        comp = rad/(scalecor*rad2)  - 1
                        if comp > tol:
                            Kr[idx] = scalecor*Kr[idx]
                            #Ai[idx] = scalecor*Ai[idx]
                        

            if foundflag == False:

                # If object was not found in the other catalog, then
                # the Kr is reduced by a factor "red". Final radius
                # can not be less than minrad
                # This is done to avoid faint large galaxies.
 

                Kr[idx]= red * Kr[idx]

                if Kr[idx]*Ai[idx] < minrad:    

                    Kr[idx] = minrad/Ai[idx]
        
                count2 +=1 

        

        line = "catalog {}: {} objects with modified radius ".format(sexcat1,count)
        print(line)


        line = "catalog {}: {} objects with reduced radius ".format(sexcat1,count2)
        print(line)


        #writing catalogs

        fout = open(newcat, "w")

        for idx, item in enumerate(N):

            line="{0:.0f} {1} {2} {3} {4} {5} {6:.2f} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

            fout.write(line)

        fout.close()



    def CheckFlag(self,val: int, check: int) -> bool:
        "Check for flag contained in val, returns True if found "

        flag = False
        mod = 1
        maxx=128


        while (mod != 0):

            res = int(val / maxx)

            if (maxx == check and res == 1):

                flag = True

            mod = val % maxx

            val = mod
            maxx = maxx / 2

        return flag


