##############################################################
##### lectura y modificacion de los archivos de Sextractor ########
#########################################################


# hot.sex

    flaghead = False

    f_out = open(outhot, "w")

    with open(sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

#every line of the file:
        for line2 in lines:

        #print line2
            (params)=line2.split()

            if params[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(dn1)

            if params[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(dm1)

            if params[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(at1)

            if params[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(dt1)

            if params[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(da1)

            if params[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME hot.cat"

            if params[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(bs1)

            if params[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(bf1)

######  check output catalog
            if params[0] ==  "CATALOG_TYPE":
                if params[1] !=  "ASCII":
	                   flaghead = True
######

            line2 = line2+"\n"
            f_out.write(line2)


    f_out.close()


# cold.sex

    f_out = open(outcold, "w")

    with open(sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

# every line of the file:
        for line2 in lines:

            (params)=line2.split()

            if params[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(dn2)

            if params[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(dm2)

            if params[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(at2)

            if params[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(dt2)

            if params[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(da2)

            if params[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME cold.cat"

            if params[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(bs2)

            if params[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(bf2)


            f_out.write(line2+"\n")



    f_out.close()

    if flaghead:
        print ("CATALOG_TYPE must be ASCII in default.sex. Ending program.. \n") # change here to assert function
        sys.exit()

##############################################################
##### ejecucion de sextractor ################### 
#########################################################



######  Running Sextractor  #######

    if (run1 == 1):
        print("Running hot.sex   \n")

        runcmd="sextractor -c {} {} ".format(outhot,image)
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT


    if (run2 == 1):
        print("Running  cold.sex  \n")

        runcmd="sextractor -c {} {} ".format(outcold,image)
        err2 = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT


