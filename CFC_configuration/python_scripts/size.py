import numpy as np
from astropy.io import fits
from read_files import *
import sys

fichero=sys.argv[1]

hdulist = fits.open(fichero)

data = hdulist[0].data

#######################
### imagen circular ###
#######################

if np.shape(data)==(1650,1650):
	_=write_files('CFC_configuration/sextractor_configuration_files/sextractor1.sex','WEIGHT_TYPE','MAP_WEIGHT'+'\n'+'WEIGHT_IMAGE CFC_configuration/masks/mask_1650.fits'+'\n')
	_=write_files('CFC_configuration/sextractor_configuration_files/sextractor2.sex','WEIGHT_TYPE','MAP_WEIGHT'+'\n'+'WEIGHT_IMAGE CFC_configuration/masks/mask_1650.fits'+'\n')
elif np.shape(data)==(1601,1601):
	_=write_files('CFC_configuration/sextractor_configuration_files/sextractor1.sex','WEIGHT_TYPE','MAP_WEIGHT'+'\n'+'WEIGHT_IMAGE CFC_configuration/masks/mask_1601.fits'+'\n')
	_=write_files('CFC_configuration/sextractor_configuration_files/sextractor2.sex','WEIGHT_TYPE','MAP_WEIGHT'+'\n'+'WEIGHT_IMAGE CFC_configuration/masks/mask_1601.fits'+'\n')

#######################
### imagen cuadrada ###
#######################

else:
	_=write_files('CFC_configuration/sextractor_configuration_files/sextractor1.sex','WEIGHT_TYPE','NONE'+'\n')
	_=write_files('CFC_configuration/sextractor_configuration_files/sextractor2.sex','WEIGHT_TYPE','NONE'+'\n')
