#! /usr/bin/env python

import sys
import subprocess as sp



def wsex(params):

    # hot.sex

    flaghead = False

    f_out = open(params.outhot, "w")

    with open(params.sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

        for line2 in lines:

            (linparams)=line2.split()

            if linparams[0] ==  "MAG_ZEROPOINT":
	            line2= "MAG_ZEROPOINT "+str(params.zpt)
            if linparams[0] ==  "GAIN":
	            line2= "GAIN "+str(params.gain)
            if linparams[0] ==  "PIXEL_SCALE":
	            line2= "PIXEL_SCALE "+str(params.plate)
            if linparams[0] ==  "SEEING_FWHM":
	            line2= "SEEING_FWHM "+str(params.seeing)



            if linparams[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(params.dn1)

            if linparams[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(params.dm1)

            if linparams[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(params.at1)

            if linparams[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(params.dt1)

            if linparams[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(params.da1)

            if linparams[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME hot.cat"

            if linparams[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(params.bs1)

            if linparams[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(params.bf1)

            if linparams[0] ==  "CATALOG_TYPE":
                if linparams[1] !=  "ASCII":
	                   flaghead = True



            line2 = line2+"\n"
            f_out.write(line2)


    f_out.close()


    # cold.sex

    f_out = open(params.outcold, "w")

    with open(params.sexfile) as f_in2:

        lines = (line.rstrip() for line in f_in2) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

# every line of the file:
        for line2 in lines:

            (linparams)=line2.split()

            if linparams[0] ==  "MAG_ZEROPOINT":
	            line2= "MAG_ZEROPOINT "+str(params.zpt)
            if linparams[0] ==  "GAIN":
	            line2= "GAIN "+str(params.gain)
            if linparams[0] ==  "PIXEL_SCALE":
	            line2= "PIXEL_SCALE "+str(params.plate)
            if linparams[0] ==  "SEEING_FWHM":
	            line2= "SEEING_FWHM "+str(params.seeing)


            if linparams[0] == "DEBLEND_NTHRESH":
                line2="DEBLEND_NTHRESH "+str(params.dn2)

            if linparams[0] == "DEBLEND_MINCONT":
                line2= "DEBLEND_MINCONT "+ str(params.dm2)

            if linparams[0] ==  "ANALYSIS_THRESH":
                line2= "ANALYSIS_THRESH "+str(params.at2)

            if linparams[0] ==  "DETECT_THRESH":
   	            line2= "DETECT_THRESH "+str(params.dt2)

            if linparams[0] ==  "DETECT_MINAREA":
	            line2= "DETECT_MINAREA "+str(params.da2)

            if linparams[0] ==  "CATALOG_NAME":
	            line2= "CATALOG_NAME cold.cat"

            if linparams[0] ==  "BACK_SIZE":
	            line2= "BACK_SIZE "+str(params.bs2)

            if linparams[0] ==  "BACK_FILTERSIZE":
	            line2= "BACK_FILTERSIZE "+str(params.bf2)


            f_out.write(line2+"\n")


    f_out.close()

    if flaghead:
        print ("CATALOG_TYPE must be ASCII in default.sex. Ending program.. \n") # change here to assert function
        sys.exit()


def runsex(params):

######  Running Sextractor  #######

    if (params.run1 == 1):
        print("Running hot.sex ")
        print("source-extractor -c {} {} \n".format(params.outhot,params.image))
        runcmd="source-extractor -c {} {} ".format(params.outhot,params.image)
        # Run Sextractor 
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  
        if err.returncode != 0:
            print("error: ",err.returncode,err.stderr)
            sys.exit("There was a problem running sextractor")

    if (params.run2 == 1):
        print("Running cold.sex ")
        print("source-extractor -c {} {} \n".format(params.outcold,params.image))
        runcmd="source-extractor -c {} {} ".format(params.outcold,params.image)
        # Run Sextractor 
        err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True) 
        if err.returncode != 0:
            print("error: ",err.returncode,err.stderr)
            sys.exit("There was a problem running sextractor")











