# -*-coding: utf-8 -*-

import numpy as np
from astropy.io import fits
from read_files import *
import sys
import os

try:
	os.remove('CFC_configuration/sextractor_result_files/test.cat')
except:
	pass

fichero_fits=sys.argv[1]

fichero=delete_folder_name(fichero_fits)
try:
	hdulist = fits.open(fichero_fits)
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Image not valid for SExtractor'+'\n')
	images_table.close
except:
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Image not found'+'\n')
	images_table.close
	exit()

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
