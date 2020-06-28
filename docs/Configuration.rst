CFC configuration
*****************

Description of the configuration files for SExtractor and PSFEx.

CFC performs a step by step process to identify the image sources and estimate their instrumental parameters, using `SExtractor <https://sextractor.readthedocs.io/en/latest/>`_ and `PSFEx <https://psfex.readthedocs.io/en/latest/>`_ with the configuration parameters located in CFC_configuration/sextractor_configuration_files.

The method is based on a first iteration of SExtractor, which is used by PSFEx to estimate the PSF of the image sources, which is used in a second iteration of SExtractor to obtain the instrumental parameters of these sources.

Operation files
===============

sextractor1.sex
---------------

Is the configuration file of the first SExtractor iteration, producing a sources catalog in FITS_LDAC format.

sextractor2.sex
---------------

Is the configuration file of the second SExtractor iteration, producing a sources catalog in ASCII_HEAD format.

psfex_config.psfex
------------------

Is the PSFEx configuration file, producing a PSF model of the image sources.

Parameters files
================

prepsfex.param
--------------

The parameters file of the firts SExtractor iteration. The parameters are indicated to use PSFEx.


postpsfex.param
---------------

The parameters file of the second SExtractor iteration. The parameters are necessary for the CFC selection criteria and calibration.
