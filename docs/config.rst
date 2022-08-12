
.. contents::
   :depth: 3
..

Configuration File
=====

Here you will find explanation of each of the parameters
of the config file needed to run CluSex.

First of all, to run CluSex with 
the configuration file is as easy as typing: 

::

   clusex ConfigFile 


An example of a configuration file is shown below::




  # run using low deblend parameters and high SNR

  #FirstRun  0  # Enable first run (1 = run)

  #SecondRun 0 # enable second run   (1 = run)

  image A1656.fits


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

  MakeMask  0

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
--------------

Let''s see the explanaiton of the parameters one by one:

# params for first run of Sextractor # run with low deblend number and
high SNR

FirstRun 1 # Enable first run (1 = run)

ANALYSIS_THRESH1 20 # or , in mag.arcsec-2

DETECT_THRESH1 20 # or , in mag.arcsec-2

DETECT_MINAREA1 10 # minimum number of pixels above threshold

DEBLEND_NTHRESH1 64 # Number of deblending sub-thresholds

DEBLEND_MINCONT1 0.001 # Minimum contrast parameter for deblending

BACK_SIZE1 100

BACK_FILTERSIZE1 11

# params for second run of Sextractor # run with high deblend number and
low SNR

SecondRun 1 # enable second run (1 = run)

ANALYSIS_THRESH2 1.5 # or , in mag.arcsec-2

DETECT_THRESH2 1.5 # or , in mag.arcsec-2

DETECT_MINAREA2 10 # minimum number of pixels above threshold

DEBLEND_NTHRESH2 16 # Number of deblending sub-thresholds

DEBLEND_MINCONT2 0.01 # Minimum contrast parapymeter for deblending

BACK_SIZE2 10

BACK_FILTERSIZE2 2

# General parameters:

Scale 1 # factor scale which ellipses are enlarged

SatDs9 sat.reg

SatScale 3

SatOffset 1

MakeMask 0

OutCatalog hotcold.cat

RegDs9 hotcold.reg


