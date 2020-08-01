from astropy.io import fits
from astropy.table import Table
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

try:
	os.mkdir('catalogs_folder/merge_catalogs')
except:
	pass


catalogs=os.listdir('catalogs_folder/CFC_catalogs')
#print(catalogs)
print('Reading catalogs')
for j in range(len(catalogs)):
	print(catalogs[j])
	if catalogs[j][len(catalogs[j])-4:len(catalogs[j])]=='fits':
		hdulist = fits.open('catalogs_folder/CFC_catalogs/'+catalogs[j])
		catalog = hdulist['catalog'].data
		#print(catalog['cl_SDSS'])
k=0
print('Merging catalogs')
medidor=0
for j in range(len(catalogs)):
	if j==round(len(catalogs)*0.1):print(10,"%")
	if j==round(len(catalogs)*0.2):print(20,"%")
	if j==round(len(catalogs)*0.3):print(30,"%")
	if j==round(len(catalogs)*0.4):print(40,"%")
	if j==round(len(catalogs)*0.5):print(50,"%")
	if j==round(len(catalogs)*0.6):print(60,"%")
	if j==round(len(catalogs)*0.7):print(70,"%")
	if j==round(len(catalogs)*0.8):print(80,"%")
	if j==round(len(catalogs)*0.9):print(90,"%")
	if j==round(len(catalogs)):print(100,"%")
	if catalogs[j][len(catalogs[j])-4:len(catalogs[j])]=='fits':
		if medidor==0:
			hdu11= fits.open('catalogs_folder/CFC_catalogs/'+catalogs[j])
			nrows = hdu11[1].data.shape[0]
			medidor=1
		else:
			hdul2=fits.open('catalogs_folder/CFC_catalogs/'+catalogs[j])
			nrows2 = hdul2[1].data.shape[0]
			nrows = nrows + nrows2
			if medidor==2:
				hdu = fits.BinTableHDU.from_columns(hdu2.columns, nrows=nrows)
			if medidor==1:
				hdu = fits.BinTableHDU.from_columns(hdu11[1].columns, nrows=nrows)
				medidor=2
			
			for colname in hdu11[1].columns.names:
				hdu.data[colname][nrows-nrows2:] = hdul2[1].data[colname]
			hdu2=hdu
			

#print(hdu.data)

sys.stdout = open(os.devnull, 'w')
hdu.writeto('catalogs_folder/merge_catalogs/total_catalog.fits',overwrite=True)
votable2=Table.read('catalogs_folder/merge_catalogs/total_catalog.fits')
votable2.write('catalogs_folder/merge_catalogs/total_catalog.xml',table_id='table_id',format='votable',overwrite=True)
sys.stdout = sys.__stdout__

