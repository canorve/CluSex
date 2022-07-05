#! /usr/bin/env python

import sys
import os.path




def readcon(params,confile):

    if not os.path.exists(confile):
        print ('%s: filename does not exist!' %sys.argv[1])
        sys.exit()


    with open(confile) as f_in:

        lines = (line.rstrip() for line in f_in) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) #remove comments
        lines = (line for line in lines if line) # Non-blank lines


        for line2 in lines:

            (param,val)=line2.split()

    # param 1

            if param == "FirstRun":

                (run1) = val.split()[0]
                params.run1=int(run1)


            if param == "DEBLEND_NTHRESH1":
                (dn1) = val.split()[0]
                params.dn1=int(dn1)

            if param == "DEBLEND_MINCONT1":
                (dm1) = val.split()[0]
                params.dm1=float(dm1)


            if param == "ANALYSIS_THRESH1":

                (at1) = val.split()[0]
                params.at1=float(at1)


            if param == "DETECT_THRESH1":

                (dt1) = val.split()[0]
                params.dt1= float(dt1)

            if param == "DETECT_MINAREA1":

                (da1) = val.split()[0]
                params.da1=int(da1)

            if param == "BACK_SIZE1":

                (bs1) = val.split()[0]
                params.bs1= float(bs1)

            if param == "BACK_FILTERSIZE1":

                (bf1) = val.split()[0]
                params.bf1=int(bf1)

    #  param2

            if param == "SecondRun":

                (run2) = val.split()[0]
                params.run2=int(run2)


            if param == "DEBLEND_NTHRESH2":

                (dn2) = val.split()[0]
                params.dn2=float(dn2)


            if param == "DEBLEND_MINCONT2":

                (dm2) = val.split()[0]
                params.dm2=float(dm2)


            if param == "ANALYSIS_THRESH2":

                (at2) = val.split()[0]
                params.at2=float(at2)


            if param == "DETECT_THRESH2":

                (dt2) = val.split()[0]
                params.dt2=float(dt2)

            if param == "DETECT_MINAREA2":

                (da2) = val.split()[0]
                params.da2=int(da2)

            if param == "BACK_SIZE2":

                (bs2) = val.split()[0]
                params.bs2= float(bs2)

    # other options

            if param == "SatDs9":
                (satfileout) = val.split()[0]
                params.satfileout=str(satfileout)

            if param == "SatScale":
                (satscale) = val.split()[0]
                params.satscale=float(satscale)

            if param == "SatOffset":

                (satoffset) = val.split()[0]
                params.satoffset=int(satoffset)

            if param == "SATUR_LEVEL":

                (satlevel) = val.split()[0]
                params.satlevel=int(satlevel)


            if param == "MakeMask":

                (create) = val.split()[0]
                params.create=int(create)

            if param == "Scale":

                (scale) = val.split()[0]
                params.scale=float(scale)

            if param == "Offset":

                (offset) = val.split()[0]
                params.offset=int(offset)



            if param == "Scale2":

                (scale2) = val.split()[0]
                params.scale2=float(scale2)


            if param == "OutCatalog":

                (output2) = val.split()[0]
                params.output2=str(output2)

            if param == "MinSatSize":   # min size for sat regions

                (minsatsize) = val.split()[0]
                params.minsatsize=int(minsatsize)

            if param == "SatQ":

                (satq) = val.split()[0]
                params.satq=float(satq)


            if param == "image":

                (image) = val.split()[0]
                params.image = str(image)


            if param == "PropTol":   # proportion tolerance

                (tol) = val.split()[0]
                params.tol=float(tol)


            if param == "RedFact":  #reduction factor

                (red) = val.split()[0]
                params.red =float(red)

            if param == "MinRad":  #reduction factor

                (minrad) = val.split()[0]
                params.minrad =float(minrad)
                params.flagminrad = True



            if param == "SatMethod":   #   method for detecting very bright saturated stars 

                (satmethod) = val.split()[0]
                params.satmethod=int(satmethod)


        # SExtractor configuration parameters

            if param == "MAG_ZEROPOINT":  # zpt for sextractor  

                (zpt) = val.split()[0]
                params.zpt=float(zpt)


            if param == "GAIN":   

                (gain) = val.split()[0]
                params.gain=float(gain)

            if param == "PIXEL_SCALE":   

                (plate) = val.split()[0]
                params.plate=float(plate)


            if param == "SEEING_FWHM":   

                (seeing) = val.split()[0]
                params.seeing=float(seeing)



            if param == "JoinScale":
                (joinscale) = val.split()[0]
                params.joinscale=float(joinscale)






