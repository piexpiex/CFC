import os

########################
### Folders creation ###
########################

try:
	os.mkdir('logouts_folder')
	os.mkdir('figures_folder')
	os.mkdir('catalogs_folder')
	os.mkdir('CFC_configuration/Intermediate_files')
	os.mkdir('CFC_configuration/sextractor_result_files')
	print('making directories')
	print('starting the analysis')
except:
	print('starting the analysis')

######################################
### Creating the logout data table ###
######################################

images_table=open('logouts_folder/data_table.txt','w')


images_table.write('Image name' +','+'number of sources' +','+ 'A (mag)' +','+  'B' +','+ 'r' +','+  'state'+ '\n'+'\n'+'\n')

images_table.close
