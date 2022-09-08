
.. contents::
   :depth: 3
..

Parameter Configuration File
===============================

Here you will find the description of each parameter
of the configuration file needed to run CluSex.

To run CluSex with the configuration 
file is as easy as typing: 

::

   clusex ConfigFile 


An example of a configuration file is shown below::


  

  FirstRun  1  # Enable first run (1 = run)

  SecondRun 1 # enable second run   (1 = run)

  Image A1656.fits


  MAG_ZEROPOINT   28.32
  GAIN            5.4
  PIXEL_SCALE     0.68
  SATUR_LEVEL     30000
  SEEING_FWHM     1.5



  DEBLEND_NTHRESH1 64          # Number of deblending sub-thresholds
  DEBLEND_MINCONT1 0.001         # Minimum contrast parameter for deblending

  ANALYSIS_THRESH1 5        # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
  DETECT_THRESH1   5          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
  DETECT_MINAREA1  20          # minimum number of pixels above threshold

  BACK_SIZE1      100
  BACK_FILTERSIZE1 11


  # params for second run
  # run with high deblend number and low SNR

  DEBLEND_NTHRESH2 32           # Number of deblending sub-thresholds
  DEBLEND_MINCONT2 .01         # Minimum contrast parapymeter for deblending

  ANALYSIS_THRESH2 1.5         # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
  DETECT_THRESH2   1.5         # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
  DETECT_MINAREA2  20      # minimum number of pixels above threshold

  BACK_SIZE2       11
  BACK_FILTERSIZE2 10 

  Scale  1.5   # factor scale which ellipses are enlarged

  Offset 5


  SatDs9 sat.reg

  SatScale 1 

  SatOffset  20

  OutCatalog  hotcold.cat

  RegDs9   hotcold.reg

  MinSatSize 20      # min size for sat regions

  SatQ 0.7

  SatMethod  3 

  ReduCoef 0.2

  FracTol 0.5


  JoinScale 2


  ScaleCor 0.7 

  ScaleCor 1.5 


Description of the Configuration file parameters
--------------------------------------------------

Let's see the explanation of each of the parameters above one by one. 


FirstRun  
        Enables to run the first run (Enable it 1; disable with  = 0)

SecondRun 
        Enables to run the second run (Enable it 1; disable with  = 0)


The previous two routines enable that first catalog 
or second catalog can run. Enable only one run it is
the same as just using sextractor alone. This will 
aid to visualize if the setup of one of the runs is
working as desire.


Image 
    The FITS cluster image.



Capital parameters are the same parameters 
as Sextractor parameters:

MAG_ZEROPOINT   28.32
GAIN            5.4
PIXEL_SCALE     0.68
SATUR_LEVEL     30000
SEEING_FWHM     1.5

They are needed for Sextractor 
to run.

For the capital parameters ending 
with 1 or 2 (like the ones below) 
refer to the same parameter that
sextractor has, but the difference 
is that 1 refer to the first run
and 2 for the second. 

Parameters of the first run:

DEBLEND_NTHRESH1 64       
DEBLEND_MINCONT1 0.001   

ANALYSIS_THRESH1 5      
DETECT_THRESH1   5     
DETECT_MINAREA1  20   

BACK_SIZE1      100
BACK_FILTERSIZE1 11


Parameters of the second run:


DEBLEND_NTHRESH2 32           # Number of deblending sub-thresholds
DEBLEND_MINCONT2 .01         # Minimum contrast parapymeter for deblending

ANALYSIS_THRESH2 1.5         # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
DETECT_THRESH2   1.5         # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
DETECT_MINAREA2  20      # minimum number of pixels above threshold

BACK_SIZE2       11
BACK_FILTERSIZE2 10 



Scale and Offset  

    CluSex defines the size of every galaxy 
    drawing a concentric ellipse. The major
    axis of this ellipse is defined by 
    Scale * Kron radius  * Ai + Offset. 
    Kron radius and Ai are parameters given
    by the output of Sextractor. 


SatDs9 
    The name of the saturation Ds9 region file. CluSex
    creates a box Ds9 saturation region file 
    where contains the saturated or bad regions of 
    the image.

SatScale and SatOffset 

    Same as the Scale and Offset parameters but 
    for the saturated regions

OutCatalog 

    The name of the output CluSex catalog

RegDs9  

    The name of the output Ds9 region file catalog.
    CluSex creates an ds9 region file from the final catalog. 
    Consequently, user can visualize the detected objects
    and their respective sizes. 

MinSatSize 

    In case is needed, user can establish 
    minimum size for the saturated region. Saturated
    regions are represented by boxes, hence the value
    of this parameters represent the side of the box. 

SatQ 

    The value of this parameter set a limit
    for the axis ratio of the saturated box. Boxes
    with axis ratio lower than this value will be 
    break it one horizontal and vertical box to
    cover the most part of the saturated regions.
    

SatMethod
    CluSex have 4 methods to identify the size of
    saturated regions. Best method is 3 (which is
    combination of methods 1 and 2), if this 
    doesn't work for you, try 4. 

ReduCoef

    This value is multiplied to the size of the objects
    if those objects are just found only in one of the
    run catalogs. A value of .2 means that the size 
    is reduced 20%
  

FracTol 

    CluSex compares the size of the same object found
    in the two run catalogs. If the difference is greater
    than this value, CluSex will modify the object size 
    keeping the smaller size of the two catalogs. A value
    of FracTol of 0.5 means that only a difference of 50% in
    radius is allowed.

ScaleCor 

    This parameter is related with the previous one, if the 
    object size is modified, the value of ScaleCor is multiplied
    by the final object size. Use it only if you think it is needed.

JoinScale 

    This parameter is the same as Scale, but this is only 
    used when CluSex will join the two catalogs and it is not
    used anymore.


Note: Not all the parameters must be in the configuration file.

