Extra features
**************

Create a customized mask
========================

An image taken with the CCD SITE and with own reduction, it is necessary to use a mask that prevents the appearance of artifacts in the limits of the circular aperture of the detector.

It is also a good way to limit detection to a specific area of the image.

To do this, use the python script located in CFC_configuration / masks mask_creator.py with the following parameters: image.fits is the image to be calibrated so that the mask has the same dimensions, X0 is the value on the horizontal axis of the pixel in the center of the mask, y0 is the value on the vertical axis of the pixel in the center of the mask, r is the radius of the mask (only circular openings, although it is easily editable to obtain other shapes).

With the following syntax:

.. code-block:: python 

   python mask_creator image.fits x0 y0 r

Verbosity
=========

To change the verbosity use a prefix command before the CFC command:

.. code-block:: bash 

   verbosity=NORMAL sh CFC.sh filabres_tree.csv

The avalaible values are FULL, NORMAL and QUIET.


How to not overwrite
====================

Due to the SExtractor and PSFEx runs take up most of the time, it may be of interest not to run the entire program again, To not the overwrite the SExtractor final catalogs use a prefix command before the CFC command:

.. code-block:: bash 

   overwrite=no sh CFC.sh filabres_tree.csv

The avalaible values are no and yes (used by default).