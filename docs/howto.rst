
.. contents::
   :depth: 3
..

--------------

**How to run**
=========================

Joincat 
~~~~~~~

Joincat is a small CluSex version. It just joins two 
existent sextractor catalogs. The aim is that a sextractor 
catalog can be merged with the output of CluSex. The aim is to 
detect those objects that were unable to be detected 
by CluSex. 

The principle is the same as CluSex: objects of the second catalog
will be added to the first one only if their center is outside the 
ellipse of the objects of the first catalog. Use it only if it is necessary. 


::

    usage: joincat [-h] [-s JOINSCALE] [-o OUTPUT] [-sf SATFILE] FirstCatalog SecondCatalog
    joincat: error: the following arguments are required: FirstCatalog, SecondCatalog



MakeMask
~~~~~~~

This routine creates an image containing ellipse masks for every object. 
It needs the CluSex output catalog and saturated ds9 regions (created by
CluSex as well)

::

    usage: makemask [-h] [-s SCALE] [-off OFFSET] [-o OUTMASK] [-sf SATFILE] [-n] SexCatalog Image
    makemask: error: the following arguments are required: SexCatalog, Image
    


Sky background
~~~~~~~~~~~~~~

This routine use two methods (gradient sky and random box) to compute
sky background for every detected object by CluSex. Output catalog
is the same as the input catalog but with the background column changed
to the new values

::

    usage: compsky [-h] [-s SCALERADIUS] [-w WIDTH] [-b BOX] [-nb NUMBOX] [-sm SCALERADMAX] [-m METHOD] [-o OUTCAT] SexCatalog Image MaskFile
    compsky: error: the following arguments are required: SexCatalog, Image, MaskFile



sex2ds9
~~~~~~~

Creates a ds9 region file from the sextractor output catalog

::
  
  usage: sex2ds9 [-h] [-s SCALE] [-off OFFSET] [-o OUTREG] SexCatalog
  sex2ds9: error: the following arguments are required: SexCatalog


