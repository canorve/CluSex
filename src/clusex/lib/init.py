
import clusex

class Params:

    # class variables
    # default values


    image = "none.fits"


    outhot  = "hot.sex"
    outcold = "cold.sex"
    sexfile= "default.sex"
    output="hc.cat"
    output2="hotcold.cat"
    scale=1
    offset = 0


    run1=run2=1

    #create = 0 #has moved to makemask


    satscale=1
    minsatsize=20
    satq=0.1
    satlevel=50000
    

    regoutfile = "hotcold.reg"

    #maskfile="mask.fits"

    SexArSort="sortar.cat"

    dn1 = dn2 = 1
    dm1 = dm2 = 1
    at1 = at2 = 1
    dt1 = dt2 = 1
    da1 = da2 = 1
    bs1 = bs2 = 1
    bf1 = bf2 = 1


    satfileout="sat.reg"


    satoffset  = 1

    satlevel = 100000

    create = 0


    output2 = "hotcold.cat"

    minsatsize = 10

    satq = 0.1


    tol = 0.5  #Fraction tolerance for differences between the two catalogs (proportion)  
    red = 0.3  # reduction factor or reduction coefficient 
    minrad = 10 # minimum radius to apply for reduction factor
    flagminrad = False # if minrad is not readed from conf file, then fwhm is set as minrad

    satmethod = 3 


    zpt = 25

    gain = 10


    plate = 1


    seeing = 2

    joinscale = 0.7


    scalecor = 1

    flagDs9 = 1 




def printWelcome():


    print("CluSex: Sextractor on Cluster Galaxies")

    print("Version:",clusex.__version__)

    url = "https://github.com/canorve/CluSex"
    
    print("webpage: "+url+"\n")


    #print("CluSex joins two sextractor catalogs. ")
    #print("It creates masks, finds saturated regions and computes sky background")




