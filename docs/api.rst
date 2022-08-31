
.. contents::
   :depth: 3
..

--------------

**API**
===================


joinsexcat
------------
::

    from clusex import joinsexcat 
    joinsexcat(firstsex,secondsex,output,joinscale)

firstsex
      This is the first sextractor catalog
secondsex 
      This is the second sextractor catalog
output
      name of the output catalog
joinscale
      the scale used for the radius of the first catalog to join 
      the two catalogs

CatArSort
------------
::

    from clusex import CatArSort 
 
    Total = CatArSort(sexcatalog,scale,offset,sexarsort,ncol,nrow)

sexcatalog
    the sextractor catalog
scale
    the scale used to enlarge (or decrease) the ellipse's axis ratio
offset
    the value added to the ellipse's axis ratio 
sexarsort
   the output sextractor catalog ordered by area in decreased order
ncol
  The number of columns (X-axis) in the image
nrow
  The number of rows  (Y-axis) in the image
returns:
Total 
    the total number of objects in the catalog


MakeMask
---------
::

    from clusex import MakeMask
    MakeMask(maskimage, catfile, scale, offset, regfile)

maskimage
    An empty FITS image. This is the output, it will overwrite it with the new
    mask
catfile
    The sextractor catalog file. The one created by CatArSort prefereably 
    but it can be a simple sextractor catalog.
scale
    the scale used to enlarge (or decrease) the ellipse's axis ratio
offset
    the value added to the ellipse's axis ratio 
regfile
    the region ds9 file containing the saturated objects

MakeSatBox
-----------
::

    from clusex import MakeSatBox 
    MakeSatBox(maskimage, region, val, ncol, nrow)

maskimage
    the mask fits image
region
   the ds9 saturated region file
val
  the value  that will be filled for those saturated regions in the mask image 
ncol
  The number of columns (X-axis) in the image
nrow
  The number of rows  (Y-axis) in the image


ds9kron
--------

::

    from clusex import ds9kron
    ds9kron(sexcatalog,regoutfile,scale,offset)
sexcatalog
    the sextractor catalog
regoutfile
    the region ds9 file containing the saturated objects
scale
    the scale used to enlarge (or decrease) the ellipse's axis ratio
offset
    the value added to the ellipse's axis ratio 


SkyCal
---------

Gradient sky computation method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    from clusex import SkyCal
    mean,std, median,rad = SkyCal().GetEllipSky(datimg,mask,x,y,
                                                        thetadeg,q,Rinit,width,
                                                        ringfile,ringmask)


datimg
  the matrix data image (not the FITS file)
mask
  the matrix data mask (not the FITS file)
x, y
  the x,y position of the objects
thetadeg
  the angle of the object. thetadeg = 0 means that major axis is aligned with the Y-axis
q
  the axis ratio of the object 
Rinit
  the initial radius where it will start to measure the gradient
width
  the width of the ring
ringfile
  the same input image but a ring is shown where the sky was measured  
ringmask
  output: the mask used to identify and mask the rings

returns:
  mean, std, median, of the sky background objects
  rad
    the radius of the major axis of the ring where sky was measured


Random box sky computation method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from clusex import SkyCal
    mean,std, median = SkyCal().RandBox(datimg,maskimg,xx,yy,
                                                thetadeg,q,Rinit,box,numbox,Rmax)

datimg
  the matrix data image (not the FITS file)
mask
  the matrix data mask (not the FITS file)
x, y
  the x,y position of the objects
thetadeg
  the angle of the object. thetadeg = 0 means that major axis is aligned with the Y-axis
q
  the axis ratio of the object 
Rinit
  the minimum radius where the random boxes will start to be positioned 
box
  the size of sides of the boxes
numbox
  the number of boxes used
Rmax
  the maximum radius where the random boxes will start to be positioned


