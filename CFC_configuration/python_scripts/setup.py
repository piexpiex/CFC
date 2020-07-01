import os
import sys
from read_files import *
#################
### Verbosity ###
#################



try:
	verb=str(sys.argv[1])
except:
	verb="NORMAL"

verbosity="NORMAL"
if verb=="NORMAL" or verb=="normal":
	verbosity="NORMAL"
if verb=="FULL" or verb=="full":
	verbosity="FULL"
if verb=="QUIET" or verb=="quiet":
	verbosity="QUIET"


_=write_files('CFC_configuration/sextractor_configuration_files/sextractor1.sex','VERBOSE_TYPE',verbosity+'\n')
_=write_files('CFC_configuration/sextractor_configuration_files/sextractor2.sex','VERBOSE_TYPE',verbosity+'\n')
_=write_files('CFC_configuration/sextractor_configuration_files/psfex_config.psfex','VERBOSE_TYPE',verbosity+'\n')

########################
### Folders creation ###
########################

print('making directories')
try:
	os.mkdir('logouts_folder')
except:
	pass
try:
	os.mkdir('figures_folder')
except:
	pass
try:
	os.mkdir('catalogs_folder')
except:
	pass
try:
	os.mkdir('catalogs_folder/SExtractor_catalogs')
except:
	pass
try:
	os.mkdir('catalogs_folder/CFC_catalogs')
except:
	pass
try:
	os.mkdir('catalogs_folder/CFC_sources')
except:
	pass
try:
	os.mkdir('CFC_configuration/Intermediate_files')
except:
	pass
try:
	os.mkdir('CFC_configuration/sextractor_result_files')
except:
	pass	
print('starting the analysis')

######################################
### Creating the logout data table ###
######################################


images_table=open('logouts_folder/data_table.csv','w')
images_table.write('Description: This table shows the results obtained for each image, including the final \n status of the image based on the quality of its magnitude versus magnitude curve calibration\n')
images_table.write('Image name,-->,Name of the image\n')
images_table.write('number of sources,-->,Number of sources for the calibration\n')
images_table.write('A,-->,Zeropoint of the calibration curve in magnitudes\n')
images_table.write('e_A,-->,Uncertainty of the zeropoint of the calibration curve in magnitudes\n')
images_table.write('B,-->,gradient of the calibration curve\n')
images_table.write('e_B,-->,Uncertainty of the gradient of the calibration curve\n')
images_table.write('r,-->,Pearson correlation coefficient\n')
images_table.write('status,-->,Status of the image, if the status is calibrated a corresponding catalog was maked\n')

images_table.write('Image name' +','+'number of sources' +','+ 'A (mag)' +','+ 'e_A (mag)' +','+   'B' +','+'e_B'+','+   'r' +','+  'status'+ '\n'+'\n'+'\n')

images_table.close
