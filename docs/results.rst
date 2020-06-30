Results
*******

Catalogs
========

CFC produces three types of catalogs:

SExtractor catalogs in plane text format with (.cat) termination.

CFC catalogs with informaton about diffrent parameters of each accepted source in the image.

CFC catalogs with informaton about diffrent parameters of each source (accepted sources and saturated or artifacts) in the image. (in process)

Figures
=======

CFC returns some figures about the calibration process for each image.

This figures are:

The FLUX_MAX vs FLUX_PSF curve where the saturated sources have red color, rejected sources with a excessively large FLUX_MAX/FLUX_PSF relation have blue color, accepted sources as real sources well measured have green color.

The SDSS magnitude vs MAG_SEX curve with the calibration relation, sources with bad morphology parameters have black color, sources excluded by a sigma clipping or without class=6 & qmode=1 if the is in SDSS field have blue color, accepted sources for the calibration have red color.

The SDSS magnitude substracted by calibrated magnitude vs SDSS magnitude scatter plot and histogram.

Logout
======

CFC produces a data table with the calibration parameters for each image and if its image was calibretad or rejected because its Pearson r parameter was lower than 0.98.