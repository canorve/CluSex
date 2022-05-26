

def ds9kron(sexfile,regfile,scale):
    "Creates ds9 region file to check output catalog "




    f_out= open(regfile, "w")

#    scale = 1
    offset=0


    flagsat=4      ## flag value when object is saturated (or close to)


#print OUT "image \n";


    line="image \n"
    f_out.write(line)


    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Star,Flg=np.genfromtxt(sexfile,delimiter="",unpack=True)

    for idx, item in enumerate(N):

        Rkron = scale * Ai[idx] * Kr[idx] + offset

#        print (Rkron)

        if Rkron == 0:

            Rkron = 1

        q = (1 - E)
        bim = q * Rkron


        check=CheckFlag(Flg[idx],flagsat)  ## check if object doesn't has saturated regions


        if  (check == False ) :


            line="ellipse({0},{1},{2},{3},{4}) # color=blue move=0 \n".format(X[idx],Y[idx],Rkron,bim[idx],Theta[idx])

            f_out.write(line)


            line2="point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={2} {3} {4} \n".format(X[idx],Y[idx],"{",np.int(N[idx]),"}")

            f_out.write(line2)

        else:


            print ("Skipping object {} one or more pixels are saturated \n".format(np.int(N[idx])))


    #        f_out.write(line)



    f_out.close()




