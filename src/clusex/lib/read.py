#########################################
##### lectura del archivo de parametros #######
#########################################




    if not os.path.exists(ConfigFile):
        print ('%s: filename does not exist!' %sys.argv[1])
        sys.exit()

    count=0


    with open(ConfigFile) as f_in:

        lines = (line.rstrip() for line in f_in) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) #remove comments
        lines = (line for line in lines if line) # Non-blank lines



        for line2 in lines:

            (param,val)=line2.split()

# param 1

            if param == "FirstRun":

                (run1) = val.split()[0]
                run1=int(run1)


            if param == "DEBLEND_NTHRESH1":
                (dn1) = val.split()[0]
                dn1=int(dn1)

            if param == "DEBLEND_MINCONT1":
                (dm1) = val.split()[0]
                dm1=float(dm1)


            if param == "ANALYSIS_THRESH1":

                (at1) = val.split()[0]
                at1=float(at1)


            if param == "DETECT_THRESH1":

                (dt1) = val.split()[0]
                dt1= float(dt1)

            if param == "DETECT_MINAREA1":

                (da1) = val.split()[0]
                da1=int(da1)

            if param == "BACK_SIZE1":

                (bs1) = val.split()[0]
                bs1= float(bs1)

            if param == "BACK_FILTERSIZE1":

                (bf1) = val.split()[0]
                bf1=int(bf1)


#  param2


            if param == "SecondRun":

                (run2) = val.split()[0]
                run2=int(run2)


            if param == "DEBLEND_NTHRESH2":

                (dn2) = val.split()[0]
                dn2=float(dn2)


            if param == "DEBLEND_MINCONT2":

                (dm2) = val.split()[0]
                dm2=float(dm2)


            if param == "ANALYSIS_THRESH2":

                (at2) = val.split()[0]
                at2=float(at2)


            if param == "DETECT_THRESH2":

                (dt2) = val.split()[0]
                dt2=float(dt2)

            if param == "DETECT_MINAREA2":

                (da2) = val.split()[0]
                da2=int(da2)

            if param == "BACK_SIZE2":

                (bs2) = val.split()[0]
                bs2= float(bs2)

# other options

            if param == "SatDs9":
                (satfileout) = val.split()[0]
                satfileout=str(satfileout)

            if param == "SatScale":
                (satscale) = val.split()[0]
                satscale=float(satscale)

            if param == "SatOffset":
                (satoffset) = val.split()[0]
                satoffset=int(satoffset)

            if param == "SatLevel":
                (satlevel) = val.split()[0]
                satlevel=int(satlevel)


            if param == "MakeMask":

                (create) = val.split()[0]
                create=int(create)

            if param == "Scale":

                (scale) = val.split()[0]
                scale=float(scale)

            if param == "Scale2":

                (scale2) = val.split()[0]
                scale2=float(scale2)



            if param == "OutCatalog":
                (output2) = val.split()[0]
                output2=str(output2)

            if param == "MinSatSize":   # min size for sat regions

                (minsatsize) = val.split()[0]
                minsatsize=int(minsatsize)

            if param == "SatQ":

                (satq) = val.split()[0]
                satq=float(satq)


# f_in.close() #not need if with is used


