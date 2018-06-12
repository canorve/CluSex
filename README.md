# Pysex

## What is pysex?

Pysex is a script written in python
to combine two Sextractor catalogs.  
Pysex runs Sextractor twice with two different configurations.

In this way, this combination allows
to extract and estimate the size of almost all possible sources in an image.

## Installation

Copy or clone this code. This code is
written for python 3.

The python libraries need to be installed are:
- numpy
- sys
- os
- subprocess
- astropy
- scipy


The main files for this code are:
- pysex.py
- config.txt
- default.conv
- default.nnw
- default.sex
- sex.param



## Quickstart

To run the code just type in the command line:
```
./pysex.py ConfigFile ImageFile
```
Where ConfigFile is the configuration parameters filename for pysex, and ImageFile is *obviously* the image.


#### Example of Configuration filename
Check the *config.txt* file that comes with the code. It is displayed here:

\# params for first run of Sextractor
\# run with low deblend number and high SNR

FirstRun  1   # Enable first run (1 = run)

ANALYSIS_THRESH1 20          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

DETECT_THRESH1   20          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

DETECT_MINAREA1  10          # minimum number of pixels above threshold

DEBLEND_NTHRESH1 64          # Number of deblending sub-thresholds

DEBLEND_MINCONT1 0.001         # Minimum contrast parameter for deblending

BACK_SIZE1      100

BACK_FILTERSIZE1 11


\# params for second run of Sextractor
\# run with high deblend number and low SNR

SecondRun 1  # enable second run   (1 = run)

ANALYSIS_THRESH2 1.5          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

DETECT_THRESH2   1.5          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

DETECT_MINAREA2  10      # minimum number of pixels above threshold

DEBLEND_NTHRESH2 16           # Number of deblending sub-thresholds

DEBLEND_MINCONT2 0.01         # Minimum contrast parapymeter for deblending

BACK_SIZE2       10

BACK_FILTERSIZE2 2

\# General parameters:

Scale  1  # factor scale which ellipses are enlarged

SatDs9 sat.reg

SatScale 3

SatOffset  1

MakeMask  0


OutCatalog  hotcold.cat

RegDs9   hotcold.reg



## How pysex works?

Pysex adds all the objects in the catalog
of the first Sextractor Run. Later, it adds
the objects of the second Sextractor Run with
the following condition: The object center
of the second run must not be inside within the ellipse of the objects of the first run.


To make pysex works properly, the first run
must be configurated with a low deblend number and high SNR, and,  on the other hand,
the second run with a high deblend number and low SNR (check manual for details to how to do this).



## Additional code
In addition, this code comes with pysex3.py
and pysexbcg.py  which basically do the same
as pysex.py but they run Sextractor 3 times for more  reliability.


## More information
Check the manual (pysex.pdf) that comes with this distribution.

#### Questions?
Do you have questions or suggestions?
Please send an email to canorve [at] gmail [dot] com

## License
This code is under the license of **GNU**
