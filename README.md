# CFC (Calibrador Fotométrico de CAFOS)

Bash and Python subroutines that perform automatic photometric calibration of images from CAFOS instrument of the 2.2m telescope at the Calar Alto Observatory (CAHA). 

# Requirements

-*[SExtractor][1]* (Source-Extractor) 2.25.0 version or later.

-*[PSFEx][2]* (PSF Extractor) 3.17.1 version or later.

-*[Astroquery][3]*.

-*[Astropy][4]*.

[1]: https://github.com/astromatic/sextractor
[2]: https://www.astromatic.net/software/psfex
[3]: https://astroquery.readthedocs.io/en/latest/
[4]: https://www.astropy.org/

# Operation mode

Description of how use the program in two different modes, one adapted to the output of filabres and the other for self-calibrated images (currently works only with images of SDSS filters).

## Filabres output

it is recommended to use *[filabres][5]* for a correct reduction and calibration of the images.

[5]: https://github.com/nicocardiel/filabres

The program works by depositing the CFC.sh file, the CFC_configuration folder and a csv file with the columns: caha_id, filter_name, program and reduction_file (for example "filabres_tree.csv") and using the following command in the UNIX terminal:

sh CFC.sh filabres_tree.csv

## Images with own reduction and astrometric calibration

The program works by depositing the CFC.sh file, the CFC_configuration folder and a folder with the images to calibrate (for example "files") and using the following command in the UNIX terminal:

sh CFC.sh files

The images have to be reduced and calibrated astrometrically.

# Image processing

-Estimation of the image dimensions and application of a mask (if the image has been trimed in a peculiar way, it is recommended to use the mask.py program to create a custom mask.) and identification of the filter used (at the moment it only works for SDSS filters).

-Use of SExtractor and PSFEx for the identification of the objects in the image and the estimation of their parameters.

-Application of different selection criteria to identify real objects and discard artifacts.

-Selection of the objects with the highest photometric quality and calibration by a comparison with the SDSS DR12 or APASS DR9 catalogs.

-Elaboration of a catalog associated with each image with the objects calibrated in magnitudes and different quality parameters.

# Merge catalogs

To merge the obtained catalogs, you can use the command:

sh CFC.sh merge

This produces a catalog with lines of all catalogs obtained previously in a new folder called catalogs_folder/merge_catalogs.

# Results

-Catalogues of objects calibrated in magnitude associated with each image of CAFOS.

-Images of source selection curves and magnitude calibration.

-Elaboration of a summary table (data_table.csv) of the photometric parameters of each image and if it has been calibrated correctly.

see the official documentation at: *[https://readthedocs.org/projects/cafos-photometry-calibrator][6]*

see the official repository at: *[https://github.com/piexpiex/CFC][7]*

[6]: https://readthedocs.org/projects/cafos-photometry-calibrator/

[7]: https://github.com/piexpiex/CFC
