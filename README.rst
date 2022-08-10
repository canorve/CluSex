.. contents::
   :depth: 3
..

CluSex
=====

What is CluSex?
--------------

CluSex is a set of routines that aids Sextractor 
to perform on images of cluster galaxies (or high 
density of objects).  

More specifically, it joins two `sextractor`_ catalogs,
 creates masks, finds saturated regions and computes 
sky background. 

.. _sextractor:https://www.astromatic.net/software/sextractor/

Why CluSex?
------------
It is hard to find an effective Sextractor configuration for
all the objects in the image (check Haussler 2007). Using 
a configuration to detect large galaxies is unable
to detect small dim galaxies, and vice versa. In addition
for regions of high density of objets, Sextractor overestimates 
the size of low surface brightness galaxies. To see those 
effects check the image below: 

.. insert images of examples


After using CluSex:


.. insert images of CluSex

In addition, performance on computation of sky 
background, creation of masks, and estimation of 
the area of saturated stars are done poorly in Sextractor. 

How it works
~~~~~~~~~~

In order to solve these problems, CluSex runs Sextractor with two
different configuration parameters: the first run detects large bright  
saturated galaxies  and the second run detects small dim galaxies. 
The combination of the two catalogs gives a better representation
of all objects of the image. It also estimates the size of saturated 
stars in the image. 

To estimate the true size of low surface brightness objects, CluSex 
compares the sizes of the detected object in each of the two catalogs.
If the object was detected only for one catalog, it is reduced by 
a constant factor introduced by the user.

As said by the Sextractor creator, output Check images are bad
to be used as a masks. In CluSex, masks are created using the 
data given by sextractor catalog. Ellipses are created for every 
detected object and these can be enlarged (or shortened) by the user.
To see the masks included the saturated stars, check the 
image below. 


.. insert images of examples

Sky background can be done poorly if objects's sizes are wrongly 
estimated or not detected at all. Also it is known (Haussler 2007)
that Sextractor overestimates the sky background. A wrong sky background 
value will produce a bad estimation of the Sersic index of the galaxy.

CluSex uses two different methods to compute sky background: 1) gradient sky
and 2) random boxes around the objects.

Gradient sky method computes the background sky in a ring around the object (check galapagos)
To locate this ring, Clusex creates concentring rings around the object and 
detect the background in every ring. When the gradient of ring sky values turns positive,
clusex stops and measure the sky in that ring. 

In the case of the random box method, clusex creates boxes of the same size located 
randomly around the object. After a given number of boxes, clusex computes the 
sky background. 

Requirements
------------

- astropy
- numpy

Installation
------------

just download the code and run

::

   pip install . 



Quickstart
----------

To run the code just type in the command line:

::

   clusex ConfigFile 

Where ConfigFile is the configuration parameters filename for pysex


.. make another page to create the configfile



Example of Configuration filename
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the *config.txt* file that comes with the code. It is displayed
here:

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

How CluSex works?
----------------

Pysex adds all the objects in the catalog of the first Sextractor Run.
Later, it adds the objects of the second Sextractor Run with the
following condition: The object center of the second run must not be
inside within the ellipse of the objects of the first run.

To make pysex works properly, the first run must be configurated with a
low deblend number and high SNR, and, on the other hand, the second run
with a high deblend number and low SNR (check manual for details to how
to do this).




Additional features 
-------------------

Full explanations of the commands below are found in

.. make additional page


Joincat 
~~~~~~~


MakeMask
~~~~~~~

Sky background
~~~~~~~~~~~~~~


sex2ds9
~~~~~~~



.. joincat,  makemask, sex2ds9, compsky   

NOTES
----
Since CluSex is designed to give 
the input catalog for my other project, actually, 
CluSex only works only for the 14 output sextractor columns below:

.. insert columns

Additional columns will be added in the future releases.

API
----

.. make another webpage for the API 


Questions?
~~~~~~~~~~

Code is far from perfect, so if you have suggestions or questions
Please send an email to canorve [at] gmail [dot] com

License
-------

This code is under the license of **GNU**
