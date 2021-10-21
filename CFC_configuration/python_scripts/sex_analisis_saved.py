# -*-coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle
from matplotlib import rcParams
from math import *
from sigma_c import *
from selector import *
from filter_identificator import *
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.patches as patches
from astropy.wcs import WCS
from astropy.coordinates import Angle
import astropy.visualization as ap_vis
import numpy as np
from matplotlib.offsetbox import AnchoredText
import os
import sys

sys.stdout = open(os.devnull, 'w')
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astroquery.simbad import Simbad
sys.stdout = sys.__stdout__

try:
	import warnings
	warnings.filterwarnings('ignore')
except:
	print('warnings active')


#############################
## Lectura de las imagenes ##
#############################

fichero=sys.argv[1]

try:
	hdulist = fits.open(fichero)
except:
	fichero=delete_folder_name(fichero)
	print(fichero+'   No avalaible')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Image not found'+'\n')
	images_table.close
	exit()

try:
	id_table=open(sys.argv[2])

	lista=[]
	images=[]
	for linea in id_table:
		medidor=0
		cuenta=0
		numero=''
		for k in range(len(linea)):
			if linea[k]==',':
				medidor=0
			elif linea[k]!=',':
				medidor=1
				numero=numero+linea[k]
			if medidor==0 and numero!='' or k==len(linea)-1:
				numero=str(numero)
				lista.append(delete_space(numero))
				numero=''
		if len(lista)>1:
				images.append(lista)
		lista=[]
	images=np.array(images)
	caha_id=images[:,0]
	filter_id=images[:,1]
	programa_id=images[:,2] 
	fich_reducido_id=images[:,3]
	
	fichero=delete_folder_name(fichero)
	name_filter=filter_id[np.where(fich_reducido_id==fichero+'\n')]#hdulist[0].header['INSFLNAM']
	name_filter=name_filter[0]
	CAHA_ID=caha_id[np.where(fich_reducido_id==fichero+'\n')]

except:
	try:
		pathtoimages=''
		medidor_path=0
		for k in range(len(fichero)):
			if fichero[k]=='/':
				midealgo=k
		pathtoimages=fichero[0:midealgo]
		id_table=open(pathtoimages+'/id.csv')
		


		lista=[]
		images=[]
		for linea in id_table:
			medidor=0
			cuenta=0
			numero=''
			for k in range(len(linea)):
				if linea[k]==',':
					medidor=0
				elif linea[k]!=',':
					medidor=1
					numero=numero+linea[k]
				if medidor==0 and numero!='' or k==len(linea)-1:
					numero=str(numero)
					lista.append(delete_space(numero))
					numero=''
			if len(lista)>1:
					images.append(lista)
			lista=[]
		images=np.array(images)
		caha_id=images[:,0]
		filter_id=images[:,1]
		#url_id=images[:,3] 
		fich_reducido_id=images[:,2]
	
		fichero=delete_folder_name(fichero)
		if len(images[0])==3:
			name_filter=filter_id[np.where(fich_reducido_id==fichero+'\n')]#hdulist[0].header['INSFLNAM']
			name_filter=name_filter[0]
			CAHA_ID=caha_id[np.where(fich_reducido_id==fichero+'\n')]	
		else:
			name_filter=filter_id[np.where(fich_reducido_id==fichero)]#hdulist[0].header['INSFLNAM']
			name_filter=name_filter[0]
			CAHA_ID=caha_id[np.where(fich_reducido_id==fichero)]
			
	except:
		fichero=delete_folder_name(fichero)
		name_filter=hdulist[0].header['INSFLNAM']
		CAHA_ID='X'

MJD=hdulist[0].header['MJD-OBS']
try:
	color=search_name(name_filter)
except:
	name_filter=hdulist[0].header['INSFLNAM']
	color=search_name(name_filter)
data = hdulist[0].data

############################################
## quitar última línea del data_table.csv ##
############################################

data_table_to_write=[]
data_table = open('logouts_folder/data_table.csv')

for linea in data_table:
	if linea[0:len(fichero)-5]==fichero[0:len(fichero)-5]:
		linea=''
	else:
		pass
	data_table_to_write.append([linea])

new_data_table = open('logouts_folder/data_table.csv','w')
for j in range(len(data_table_to_write)):
	new_data_table.write(data_table_to_write[j][0])
############################################

if color=='X':
	print('no SDSS filter')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'No SDSS filter'+'\n')
	images_table.close
	exit()
#hdulist = fits.open("CFC_configuration/sextractor_result_files/bkg.fits")
#data_bkg= hdulist[0].data
#data_sub=data-data_bkg
# show the image


########################
### Reading test.cat ###
########################

try:
	catalogo = open('catalogs_folder/SExtractor_catalogs/'+fichero[0:len(fichero)-5]+'_sex.cat')
except:
	print('no saved Sextractor catalog for'+ fichero)
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'No saved Sextractor catalog'+'\n')
	images_table.close
	exit()


NUMBER=0       #   1 NUMBER                 Running object number                                     
FLAGS=1        #   2 FLAGS                  Extraction flags                                          
FLAGS_WEIGHT=2 #   3 FLAGS_WEIGHT           Weighted extraction flags                                 
SNR_WIN=3      #   4 SNR_WIN                Gaussian-weighted SNR                                     
X_IMAGE=4      #   5 X_IMAGE                Object position along x                                    [pixel]
Y_IMAGE=5      #   6 Y_IMAGE                Object position along y                                    [pixel]
A_IMAGE=6      #   7 A_IMAGE                Profile RMS along major axis                               [pixel]
B_IMAGE=7      #   8 B_IMAGE                Profile RMS along minor axis                               [pixel]
THETA_IMAGE=8  #   9 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
FWHM_IMAGE=9   #  10 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
FLUX_MAX=10    #  11 FLUX_MAX               Peak flux above background                                 [count]
FLUX_RADIUS=11 #  12 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
ALPHA_J2000=12 #  13 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
DELTA_J2000=13 #  14 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
FWHM_WORLD=14  #  15 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
FLUX_PSF=15    #  16 FLUX_PSF               Flux from PSF-fitting                                      [count]
FLUXERR_PSF=16 #  17 FLUXERR_PSF            RMS flux error for PSF-fitting                             [count]
MAG_PSF=17     #  18 MAG_PSF                Magnitude from PSF-fitting                                 [mag]
MAGERR_PSF=18  #  19 MAGERR_PSF             RMS magnitude error from PSF-fitting                       [mag]
ELONGATION=19  #  20 ELONGATION             A_IMAGE/B_IMAGE                                           
ELLIPTICITY=20 #  21 ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
CLASS_STAR=21  #  22 CLASS_STAR             S/G classifier output                                     
SPREAD_MODEL=22#  23 SPREAD_MODEL           Spread parameter from model-fitting   
X_WORLD=23##  24 X_WORLD                Barycenter position along world x axis                     [deg]
Y_WORLD=24#  25 Y_WORLD                Barycenter position along world y axis                     [deg]
ERRX2_WORLD=25#  26 ERRX2_WORLD            Variance of position along X-WORLD (alpha)                 [deg**2]
ERRY2_WORLD=26#  27 ERRY2_WORLD            Variance of position along Y-WORLD (delta)                 [deg**2]
XWIN_IMAGE=27
YWIN_IMAGE=28
XMIN_IMAGE=29
YMIN_IMAGE=30
XMAX_IMAGE=31
YMAX_IMAGE=32
FLUX_AUTO=33
FLUXERR_AUTO=34
MAG_AUTO=35
MAGERR_AUTO=36



lista=[]
objects=[]

for linea in catalogo:

	if linea[0]=='#':
		continue
	else:
		medidor=0
		cuenta=0
		numero=''
		for k in range(len(linea)):
			if linea[k]==' ':
				medidor=0
			elif linea[k]!=' ':
				medidor=1
				numero=numero+linea[k]
			if medidor==0 and numero!='' or k==len(linea)-1:

				numero=float(numero)
				lista.append(numero)
				numero=''
		objects.append(lista)
		lista=[]
objects=np.array(objects)
total_objects=objects.copy()

len_objects_key=0
try:
	len_objects=len(objects[:,0])
	len_objects_key=1
except:
	pass

if len_objects_key==1:
	if len_objects<6:
		print('Not enough objects in the image')
		images_table=open('logouts_folder/data_table.csv','a')
		images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects in the image'+'\n')
		images_table.close
		exit()
if len_objects_key==0:
	print('Not enough objects in the image')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects in the image'+'\n')
	images_table.close
	exit()
Nobjetos=0

listaok=np.array([1.0]*len(objects))

for i in range(len(objects)):
	if objects[i,FLUX_MAX]<=0 or objects[i,FLUX_RADIUS]<=0 or objects[i,FWHM_IMAGE]<=0 or objects[i,FLUX_PSF]<=0:
		listaok[i]=0
	if objects[i,MAG_PSF]==99.0:
		listaok[i]=0
	if objects[i,MAGERR_PSF]==0 or objects[i,MAGERR_PSF]>1 or objects[i,MAGERR_PSF]==99:
		listaok[i]=0
	if objects[i,FLAGS_WEIGHT]==2:
		listaok[i]=0
	if objects[i,SNR_WIN]<5 or objects[i,SNR_WIN]>0.99*10**30:
		listaok[i]=0
	if np.shape(data)==(1650,1650):
		if ((abs(objects[i,X_IMAGE]-825)+objects[i,A_IMAGE])**2+ (abs(objects[i,Y_IMAGE]-825)+objects[i,A_IMAGE])**2)**0.5>814: 
			listaok[i]=0
	if np.shape(data)==(1601,1601):
		if ((abs(objects[i,X_IMAGE]-808)+objects[i,A_IMAGE])**2+ (abs(objects[i,Y_IMAGE]-817)+objects[i,A_IMAGE])**2)**0.5>814:
			listaok[i]=0
	if listaok[i]==1.0:
		Nobjetos=Nobjetos+1
		listaok[i]=1.0

flux_coef=objects[:,FLUX_MAX]/objects[:,FLUX_PSF]
C_1=np.where(listaok==1)

flux_mean=np.mean(flux_coef[C_1])
flux_sigma=np.std(flux_coef[C_1])

C_2=np.where((flux_coef<flux_mean+2*flux_sigma) & (listaok==1))

flux_mean2=np.mean(flux_coef[C_2])
flux_sigma2=np.std(flux_coef[C_2])

for i in range(len(objects)):
	if listaok[i]==1:
		if flux_coef[i]>flux_mean2+3*flux_sigma2:
			listaok[i]=-1

max_flux=objects[np.where(listaok==1),FLUX_MAX]
psf_flux=objects[np.where(listaok==1),FLUX_PSF]
FLUX_step=29000
FLUX_step_d=2000
FLUX_SATURATION_1=60000
while abs(FLUX_step_d)>100 and FLUX_step<65000: # Pearson r criteria
	FLUX_step=FLUX_step+FLUX_step_d
	r=ajuste_lineal(psf_flux[np.where(max_flux<FLUX_step)],max_flux[np.where(max_flux<FLUX_step)])[4]
	if FLUX_step>61000:
		FLUX_step=60000
		break
	if r>=0.98:
		continue
	if r<0.98:
		FLUX_SATURATION_1=FLUX_step
		FLUX_step=FLUX_step-FLUX_step_d
		FLUX_step_d=FLUX_step_d/2

FLUX_SATURATION_1=FLUX_step-FLUX_step_d	
FLUX_step=29000
FLUX_step_d=2000
FLUX_SATURATION_2=60000
A0=ajuste_lineal(psf_flux[np.where(max_flux<FLUX_step)],max_flux[np.where(max_flux<FLUX_step)])[0]
while abs(FLUX_step_d)>100 and FLUX_step<65000: # Slope criteria
	FLUX_step=FLUX_step+FLUX_step_d
	A=ajuste_lineal(psf_flux[np.where(max_flux<FLUX_step)],max_flux[np.where(max_flux<FLUX_step)])[0]
	if FLUX_step>61000:
		FLUX_step=60000
		break
	if A>=A0/1.01:
		continue
	if A<A0/1.01:
		FLUX_SATURATION_2=FLUX_step
		FLUX_step=FLUX_step-FLUX_step_d
		FLUX_step_d=FLUX_step_d/2

FLUX_SATURATION_2=FLUX_step-FLUX_step_d	
FLUX_SATURATION=min([FLUX_SATURATION_1,FLUX_SATURATION_2])-1*flux_sigma2

for i in range(len(objects)):
	if objects[i][FLUX_MAX]>FLUX_SATURATION and listaok[i]==1.0:
		listaok[i]=2
	elif (int(float(objects[i][FLAGS]/4)) % 2)==1 and listaok[i]==1.0:
		listaok[i]=2
		
#Double check for no linear weak sources (they were taken as real sources)
##########################################################################
#C_1=np.where(listaok==1)
#
#flux_mean=np.mean(flux_coef[C_1])
#flux_sigma=np.std(flux_coef[C_1])
#
#C_2=np.where((flux_coef<flux_mean+2*flux_sigma) & (listaok==1))
#
#flux_mean2=np.mean(flux_coef[C_2])
#flux_sigma2=np.std(flux_coef[C_2])
#
#for i in range(len(objects)):
#	if listaok[i]==1:
#		if flux_coef[i]>flux_mean2+3*flux_sigma2:
#			listaok[i]=-1
##########################################################################

plt.figure(figsize=(22.0,7.0))

plt.plot(objects[np.where(listaok==2),FLUX_PSF], objects[np.where(listaok==2),FLUX_MAX],'r.')
plt.plot(objects[np.where(listaok==1),FLUX_PSF], objects[np.where(listaok==1),FLUX_MAX],'g.')
plt.plot(objects[np.where(listaok==-1),FLUX_PSF], objects[np.where(listaok==-1),FLUX_MAX],'b.')
plt.plot(objects[np.where(listaok==2),FLUX_PSF][0], objects[np.where(listaok==2),FLUX_MAX][0],'r.',label='Saturated')
plt.plot(objects[np.where(listaok==1),FLUX_PSF][0], objects[np.where(listaok==1),FLUX_MAX][0],'g.',label='Valid sources')
plt.plot(objects[np.where(listaok==-1),FLUX_PSF][0], objects[np.where(listaok==-1),FLUX_MAX][0],'b.',label='Artifacts')

#plt.plot([1*10**4,10**5,max(f1)],ffit([1*10**4,10**5,max(f1)]),'k')
plt.suptitle('FLUX MAX vs FLUX PSF (SNR >=5)')
plt.ylabel('FLUX MAX [count]')
plt.xlabel('FLUX PSF [count]')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig('figures_folder/'+fichero[0:len(fichero)-5]+'_flux_selection.pdf')


############################
######### Skymatch #########
############################
numeros=objects[:,NUMBER]
objects=objects[np.where(listaok==1)]
source_flag=[3.0]*len(listaok)
source_flag=np.array(source_flag)
source_flag[np.where(listaok==0)]=0
source_flag[np.where(listaok==-1)]=1
source_flag[np.where(listaok==2)]=2

SM_flag=objects[:,SPREAD_MODEL]
#for i in range(len(objects)):
#	if objects[i,SPREAD_MODEL]<-0.05:
#		SM_flag[i]=0.0
cl_SDSS=np.array([0.0]*len(objects[:,SPREAD_MODEL]))
final_objects=objects

objects=objects[np.where(SM_flag>-0.05)] 

if len(objects[:,0])<6:
	print('Not enough objects in the image')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects in the image'+'\n')
	images_table.close
	exit()
	
c1 = fits.Column(name='NUMBER', array=objects[:,NUMBER], format='E')
c2 = fits.Column(name='SNR_WIN',array=objects[:,SNR_WIN], format='E')
c3 = fits.Column(name='FLUX_MAX',array=objects[:,FLUX_MAX], format='E')
c4 = fits.Column(name='ALPHA',array=objects[:,ALPHA_J2000], format='E')
c5 = fits.Column(name='DELTA',array=objects[:,DELTA_J2000], format='E')
c6 = fits.Column(name='FWHM_WORLD',array=objects[:,FWHM_WORLD], format='E')
c7 = fits.Column(name='ELONGATION',array=objects[:,ELONGATION], format='E')
c8 = fits.Column(name='ELLIPTICITY',array=objects[:,ELLIPTICITY], format='E')
c9 = fits.Column(name='FLUX_PSF',array=objects[:,FLUX_PSF], format='E')
c10 = fits.Column(name='FLUXERR_PSF',array=objects[:,FLUXERR_PSF], format='E')
c11 = fits.Column(name='MAG_PSF',array=objects[:,MAG_PSF], format='E')
c12 = fits.Column(name='MAGERR_PSF',array=objects[:,MAGERR_PSF], format='E')
c13 = fits.Column(name='SPREAD_MODEL',array=objects[:,SPREAD_MODEL], format='E')

alpha_find_sources=objects[:,ALPHA_J2000]
delta_find_sources=objects[:,DELTA_J2000]

#t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8],name='valores')
#t.writeto('a.fits',overwrite=True)

t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13],name='valores')
t.writeto('CFC_configuration/Intermediate_files/parametros.fits',overwrite=True)

hdul = fits.open('CFC_configuration/Intermediate_files/parametros.fits')
delta=hdul['VALORES']
t = Table.read(delta)
t.write('CFC_configuration/Intermediate_files/parametros.vot', table_id='updated_table', format='votable',overwrite=True)

sdss_key=0
catalog = XMatch.query(cat1=open('CFC_configuration/Intermediate_files/parametros.vot'),
                         cat2='vizier:V/147/sdss12',
                         max_distance=2 * u.arcsec,
                         colRA1='ALPHA', colDec1='DELTA')

if len(catalog)==0:
	print('No SSDS field, working with APASS DR9 catalog')
	catalog = XMatch.query(cat1=open('CFC_configuration/Intermediate_files/parametros.vot'),
                         cat2='vizier:II/336/apass9',
                         max_distance=2 * u.arcsec,
                         colRA1='ALPHA', colDec1='DELTA')
	sdss_key=1


if len(catalog)<6:
	if len(catalog)==0:
		motivo='no SDSS/APASS coverage'
		print('no SDSS/APASS coverage')
	else:
		motivo='Not avalaible cross-match catalog in this skyfield'
		print('Not avalaible cross-match catalog in this skyfield')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ motivo +'\n')
	images_table.close
	
	exit()	

#Columns selection
NUMBER_XMATCH=catalog['NUMBER']
for j in range(len(NUMBER_XMATCH)):
	source_flag[int(NUMBER_XMATCH[j]-1)]=4
	
mag_sex=catalog['MAG_PSF']
magerr_sex=catalog['MAGERR_PSF']
ellongation=catalog['ELONGATION']
ellipticity=catalog['ELLIPTICITY']
FWHM=catalog['FWHM_WORLD']
SPREAD_VALUE=catalog['SPREAD_MODEL']
limit_detection=25
limit_sat=0
if sdss_key==0:
	alpha_2_find_sources=catalog['RAdeg']
	delta_2_find_sources=catalog['DEdeg']
	q_mode=catalog['q_mode']
	class_sdss=catalog['class']
	mode=catalog['mode']
	if color=='u'or color=='U':
		pmag=catalog['umag']
		e_pmag=catalog['e_umag']
		name_mag='SDSS umag'
		limit_detection=22.15
		limit_sat=13
	if color=='r'or color=='R':
		pmag=catalog['rmag']
		e_pmag=catalog['e_rmag']
		name_mag='SDSS rmag'
		limit_detection=22.70
		limit_sat=14
	if color=='g'or color=='G':
		pmag=catalog['gmag']
		e_pmag=catalog['e_gmag']
		name_mag='SDSS gmag'
		limit_detection=23.13
		limit_sat=14
	if color=='i'or color=='I':
		pmag=catalog['imag']
		e_pmag=catalog['e_imag']
		name_mag='SDSS imag'
		limit_detection=22.2
		limit_sat=14
	if color=='z'or color=='Z':
		pmag=catalog['zmag']
		e_pmag=catalog['e_zmag']
		name_mag='SDSS zmag'
		limit_detection=20.71
		limit_sat=12
if sdss_key==1:
	alpha_2_find_sources=catalog['RAJ2000']
	delta_2_find_sources=catalog['DEJ2000']
	if color=='r'or color=='R':
		pmag=catalog['rpmag']
		e_pmag=catalog['e_rpmag']
		name_mag='APASS rmag'
		limit_detection=17 #Vmag
		limit_sat=7 #Vmag
		Vmag=catalog['Vmag']
		e_Vmag=catalog['e_Vmag']
	elif color=='v'or color=='V':
		exit()
		pmag=catalog['Vmag']
		e_pmag=catalog['e_Vmag']
		name_mag='APASS Vmag'
		limit_detection=25
		limit_sat=0
	elif color=='b'or color=='B':
		exit()
		pmag=catalog['Bmag']
		e_pmag=catalog['e_Bmag']
		name_mag='APASS Bmag'
		limit_detection=25
		limit_sat=0
	elif color=='g'or color=='G':
		pmag=catalog['gpmag']
		e_pmag=catalog['e_gpmag']
		name_mag='APASS gmag'
		limit_detection=17 #Vmag
		limit_sat=7 #Vmag
		Vmag=catalog['Vmag']
		e_Vmag=catalog['e_Vmag']
	elif color=='i'or color=='I':
		pmag=catalog['ipmag']
		e_pmag=catalog['e_ipmag']
		name_mag='APASS imag'
		limit_detection=17 #Vmag
		limit_sat=7 #Vmag
		Vmag=catalog['Vmag']
		e_Vmag=catalog['e_Vmag']
	else:
		print('no SDSS/APASS coverage')
		images_table=open('logouts_folder/data_table.csv','a')
		images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'no SDSS/APASS coverage' +'\n')
		images_table.close
		exit()

if sdss_key==0:
	lista=[mag_sex,magerr_sex,ellongation,ellipticity,FWHM,pmag,e_pmag,SPREAD_VALUE,NUMBER_XMATCH,class_sdss,q_mode,mode]
if sdss_key==1:
	lista=[mag_sex,magerr_sex,ellongation,ellipticity,FWHM,pmag,e_pmag,SPREAD_VALUE,NUMBER_XMATCH,Vmag,e_Vmag]

#search with 5 digits
#alpha_find_sources=np.around(alpha_find_sources,5)
#delta_find_sources=np.around(delta_find_sources,5)

#bestmatch
#Numbers,Numbers_ok=find_sources(alpha_find_sources,delta_find_sources,alpha_2_find_sources,delta_2_find_sources)

#for k in range(len(lista)):
#	lista[k]=lista[k][Numbers]
#	lista[k]=lista[k][np.where(Numbers_ok==1)]

mag_sex=lista[0]
magerr_sex=lista[1]
ellongation=lista[2]
ellipticity=lista[3]
FWHM=lista[4]
pmag=lista[5]
e_pmag=lista[6]
SPREAD_VALUE=lista[7]
NUMBER_XMATCH=lista[8]
if sdss_key==0:
	class_sdss=lista[9]
	q_mode=lista[10]
	mode=lista[11]
if sdss_key==1:
	Vmag=lista[9]
	e_Vmag=lista[10]


if len(pmag)<6:
	print('Not enough objects for the calibration')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects for the calibration'+'\n')
	images_table.close
	exit()
print('number of skymatch objects',len(pmag))

if sdss_key==0:
	for j in range(len(class_sdss)):
		cl_SDSS[np.where(final_objects[:,0]==NUMBER_XMATCH[j])]=class_sdss[j]

plt.figure(figsize=(22.0,7.0))
plt.suptitle('Magnitude calibration')
plt.xlabel('MAG PSF ('+name_filter+')')
plt.ylabel(name_mag)
N_A=len(pmag)
plt.plot(mag_sex,pmag,'k.',label='skymatch objects('+str(N_A)+')')

###############################
### Photometric calibration ###
###############################

#error mag selection
if sdss_key==0:
	for k in range(len(lista)):
		lista[k]=lista[k][np.where((magerr_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]
else:
	for k in range(len(lista)):
		lista[k]=lista[k][np.where((magerr_sex<0.2) & (np.isnan(pmag)==False) & (Vmag<limit_detection) & (Vmag>limit_sat) & (e_pmag<0.2))]

#if sdss_key==0:
#	class_sdss=class_sdss[np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]
#	q_mode=q_mode[np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]

mag_sex=lista[0]
magerr_sex=lista[1]
ellongation=lista[2]
ellipticity=lista[3]
FWHM=lista[4]
pmag=lista[5]
e_pmag=lista[6]
SPREAD_VALUE=lista[7]
NUMBER_XMATCH=lista[8]
if sdss_key==0:
	class_sdss=lista[9]
	q_mode=lista[10]
	mode=lista[11]
#Morphology selection
if len(pmag)<6:
	print('Not enough objects for the calibration')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ ' ' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects for the calibration'+'\n')
	images_table.close
	exit()
median_FWHM=np.median(FWHM)
Rq_FWHM=np.percentile(FWHM,75)-np.percentile(FWHM,25)
median_ellongation=np.median(ellongation)
Rq_ellongation=np.percentile(ellongation,75)-np.percentile(ellongation,25)
median_ellipticity=np.median(ellipticity)
Rq_ellipticity=np.percentile(ellipticity,75)-np.percentile(ellipticity,25)

for k in range(len(lista)):
	lista[k]=lista[k][np.where((FWHM<median_FWHM+2*Rq_FWHM) & (FWHM>median_FWHM-2*Rq_FWHM)
	&(ellongation<median_ellongation+2*Rq_ellongation) & (ellongation>median_ellongation-2*Rq_ellongation)
	&(ellipticity<median_ellipticity+2*Rq_ellipticity) &(ellipticity>median_ellipticity-2*Rq_ellipticity))]
#if sdss_key==0:
#	class_sdss=class_sdss[np.where((FWHM<median_FWHM+2*Rq_FWHM) & (FWHM>median_FWHM-2*Rq_FWHM)
#	&(ellongation<median_ellongation+2*Rq_ellongation) & (ellongation>median_ellongation-2*Rq_ellongation)
#	&(ellipticity<median_ellipticity+2*Rq_ellipticity) &(ellipticity>median_ellipticity-2*Rq_ellipticity))]
#	q_mode=q_mode[np.where((FWHM<median_FWHM+2*Rq_FWHM) & (FWHM>median_FWHM-2*Rq_FWHM)
#	&(ellongation<median_ellongation+2*Rq_ellongation) & (ellongation>median_ellongation-2*Rq_ellongation)
#	&(ellipticity<median_ellipticity+2*Rq_ellipticity) &(ellipticity>median_ellipticity-2*Rq_ellipticity))]

mag_sex=lista[0]
magerr_sex=lista[1]
pmag=lista[5]
e_pmag=lista[6]
SPREAD_VALUE=lista[7]
NUMBER_XMATCH=lista[8]
if sdss_key==0:
	class_sdss=lista[9]
	q_mode=lista[10]
	mode=lista[11]
for j in range(len(NUMBER_XMATCH)):
	source_flag[int(NUMBER_XMATCH[j]-1)]=5

print('number of objects after morphology criteria',len(mag_sex))

N_B=len(pmag)
plt.plot(mag_sex,pmag,'b.',label='objects with morphology criteria ('+str(N_B)+')')

if sdss_key==0:
	for k in range(len(lista)):
		try:
			lista[k]=lista[k][np.where((class_sdss==6) & (q_mode=='+') & (mode==1.0))]
		except:
			lista[k]=lista[k][np.where((class_sdss==6) & (q_mode==1.0) & (mode==1.0))]

	mag_sex=lista[0]
	magerr_sex=lista[1]
	pmag=lista[5]
	e_pmag=lista[6]
	SPREAD_VALUE=lista[7]
	NUMBER_XMATCH=lista[8]
	class_sdss=lista[9]
	for j in range(len(class_sdss)):
		cl_SDSS[np.where(final_objects[:,0]==NUMBER_XMATCH[j])]=6
	print('number of objects with SDSS class=6 & q mode=1.0 ',len(mag_sex))

if len(mag_sex)<6:
	print('Insufficient number of objects for photometric calibration')
	if sdss_key==0:
		images_table=open('logouts_folder/data_table.csv','a')
		images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ 'SDSS' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects for the calibration'+'\n')
		images_table.close
	else:
		images_table=open('logouts_folder/data_table.csv','a')
		images_table.write(fichero[0:len(fichero)-5] +','+  ' ' +','+ 'APASS' +','+  ' ' +','+  ' ' +','+  ' ' +','+ ' ' +','+ ' ' +','+ 'rejected'+','+ 'Not enough objects for the calibration'+'\n')
		images_table.close	
	exit()
	


#Photometry calibration

#Method of mean two extreme curves
#X_up,Y_up,Z_up=sigma_c(X=mag_sex,Y=pmag+e_pmag,n_sigma=2)
#X_down,Y_down,Z_down=sigma_c(X=mag_sex,Y=pmag-e_pmag,n_sigma=2)
#X=(X_up+X_down)*0.5
#Y=(Y_up+Y_down)*0.5
#Z=(Z_up+Z_down)*0.5
#e_A=abs(Z_up[2]-Z_down[2])/2 
#e_B=abs(Z_up[3]-Z_down[3])/2 

#Classic method with Monte Carlo uncertainties
X,Y,NUMBER_XMATCH,Z,lista_sigma_c=sigma_c(X=mag_sex,Y=pmag,idfs=NUMBER_XMATCH,n_sigma=2)
n_times=500
X_MC=np.repeat(X,n_times)+np.repeat(magerr_sex[lista_sigma_c],n_times)*np.random.normal(size=np.size(X)*n_times)
Y_MC=np.repeat(Y,n_times)+np.repeat(e_pmag[lista_sigma_c],n_times)*np.random.normal(size=np.size(Y)*n_times)
Z_MC=ajuste_lineal(X_MC,Y_MC,W=0)
e_A=abs(Z_MC[2])
e_B=abs(Z_MC[3])

for j in range(len(NUMBER_XMATCH)):
	source_flag[int(NUMBER_XMATCH[j]-1)]=6

N_C=len(X)
plt.plot(X,Y,'r.',label='Calibration objects ('+str(N_C)+')')
plt.plot(mag_sex,Z[1]+Z[0]*mag_sex,'r',label=name_mag+'={0:.3g}'.format(Z[0])+' × MAG PSF + {0:.3g}'.format(Z[1]) + ' (r={0:.3g}'.format(Z[4])+')')
plt.legend()
#plt.ylim(limit_sat,limit_detection+1)

plt.savefig('figures_folder/'+fichero[0:len(fichero)-5]+'_magnitude_calibration.pdf')


###########################################
### Comparison between both magnitudes ###
###########################################

comp_mag=Y - Z[1]-Z[0]*X
plt.figure(figsize=(22.0,7.0))
plt.subplot(1,2,1)
for j in range(len(comp_mag)):
	plt.plot(pmag[j], comp_mag[j],'r.',markersize=6)
plt.ylabel(name_mag+' - MAG_PSF ('+name_filter+')')
plt.xlabel(name_mag)
plt.subplot(1,2,2)
n,bins,patches=plt.hist(comp_mag,bins=np.arange(-1,1.05,0.05),density=False,facecolor='r')
plt.ylim(0,max(n)*1.2)
plt.xlabel(name_mag+' - MAG_PSF ('+name_filter+')')
plt.savefig('figures_folder/'+fichero[0:len(fichero)-5]+'_magnitude_comparison.pdf')



##############
### logout ###
##############

state='rejected'
semi_state='r Pearson coefficient < 0.98'
if Z[4]>=0.98:
	state='calibrated'
	semi_state=' '
images_table=open('logouts_folder/data_table.csv','a')

if sdss_key==0:
	images_table.write(fichero[0:len(fichero)-5] +','+  str(N_C) +','+'SDSS'+','+  str(round(Z[1],5))+','+  str(round(e_B,5)) +','+  str(round(Z[0],5))+','+  str(round(e_A,5)) +','+ str(round(Z[4],3)) +','+  state+ ',' + semi_state+'\n')
if sdss_key==1:
	images_table.write(fichero[0:len(fichero)-5] +','+  str(N_B) +','+'APASS'+','+  str(round(Z[1],5))+','+  str(round(e_B,5)) +','+  str(round(Z[0],5))+','+  str(round(e_A,5)) +','+ str(round(Z[4],3)) +','+  state+ ',' + semi_state+'\n')

images_table.close
###########
if Z[4]>=0.98:
	calibration_mag=Z[1]+Z[0]*final_objects[:,17]
	calibration_mag_error=(abs(e_B)**2+abs(e_A*final_objects[:,17])**2+abs(Z[0]*final_objects[:,18])**2)**0.5
	min_pmag=min(Y)
	max_pmag=max(Y)

	extrapolation_mag=np.array(['A']*len(final_objects))
	for i in range(len(final_objects)):
		if calibration_mag[i]>max_pmag:
			extrapolation_mag[i]='C'
		elif calibration_mag[i]<min_pmag:
			extrapolation_mag[i]='B'
	if CAHA_ID=='X':
		c1 = fits.Column(name='Image_identifier', array=np.array(len(final_objects[:,0])*['CAHA_CAFOS_BBI_DR1']), format='50A')
	else:
		c1 = fits.Column(name='Image_identifier', array=np.array(len(final_objects[:,0])*['CAHA_CAFOS_BBI_DR1_'+str(CAHA_ID[0])]), format='50A')
	
	
	NUMBER_ID=np.arange(1,1+len(final_objects[:,0]),1).astype(np.str)
	DETECTION_ID=np.array(len(NUMBER_ID)*['CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_0000'])
	if CAHA_ID=='X':
		for j in range(len(NUMBER_ID)):
			DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+'0'*(3-int(np.log10(1+j)))+NUMBER_ID[j]
	else:
		for j in range(len(NUMBER_ID)):
			DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_'+'0'*(3-int(np.log10(1+j)))+NUMBER_ID[j]
	c2 = fits.Column(name='Detection_ID', array=DETECTION_ID, format='50A')
	MJD_array=np.zeros(len(final_objects[:,0]))
	MJD_array[:]=MJD
	c3 = fits.Column(name='MJD', array=MJD_array, format='D')
	c4 = fits.Column(name='SNR_WIN',array=final_objects[:,SNR_WIN], format='E')
	#c5 = fits.Column(name='RAJ2000', unit='deg',array=np.around(final_objects[:,ALPHA_J2000],5), format='E')
	#c6 = fits.Column(name='DEJ2000', unit='deg',array=np.around(final_objects[:,DELTA_J2000],5), format='E')
	#c7 = fits.Column(name='e_RAJ2000', unit='arcsec',array=np.around(3600*final_objects[:,ERRX2_WORLD]**0.5,5), format='E')
	#c8 = fits.Column(name='e_DEJ2000', unit='arcsec',array=np.around(3600*final_objects[:,ERRY2_WORLD]**0.5,5), format='E')
	c5 = fits.Column(name='RAJ2000', unit='deg',array=final_objects[:,ALPHA_J2000], format='E')
	c6 = fits.Column(name='DEJ2000', unit='deg',array=final_objects[:,DELTA_J2000], format='E')
	c7 = fits.Column(name='e_RAJ2000', unit='arcsec',array=3600*final_objects[:,ERRX2_WORLD]**0.5, format='E')
	c8 = fits.Column(name='e_DEJ2000', unit='arcsec',array=3600*final_objects[:,ERRY2_WORLD]**0.5, format='E')
	RA=Angle(final_objects[:,ALPHA_J2000]* u.deg)
	DEC=Angle(final_objects[:,DELTA_J2000]* u.deg)
	e_RA=Angle(final_objects[:,ERRX2_WORLD]**0.5* u.deg)
	e_DEC=Angle(final_objects[:,ERRY2_WORLD]**0.5* u.deg)
	c9 = fits.Column(name='RA_hms', unit='hh:mm:ss', array=RA.to_string(unit=u.hourangle, sep=(':',':')), format='20A')
	c10 = fits.Column(name='DE_dms', unit='dd:mm:ss', array=DEC.to_string(unit=u.deg, sep=(':',':')), format='20A')
	c11 = fits.Column(name='e_RA_hms', unit='hh:mm:ss', array=e_RA.to_string(unit=u.hourangle, sep=(':',':')), format='20A')
	c12 = fits.Column(name='e_DE_dms', unit='dd:mm:ss', array=e_DEC.to_string(unit=u.deg, sep=(':',':')), format='20A')
	c13 = fits.Column(name='MAG',array=np.around(calibration_mag,3), format='E')
	c14 = fits.Column(name='e_MAG',array=np.around(calibration_mag_error,3), format='E')
	c15 = fits.Column(name='MAG_sex',array=np.around(final_objects[:,MAG_PSF],3), format='E')
	c16 = fits.Column(name='e_MAG_sex',array=np.around(final_objects[:,MAGERR_PSF],3), format='E')
	c17 = fits.Column(name='cl_SDSS',array=cl_SDSS, format='E')
	c18 = fits.Column(name='SPREAD_MODEL',array=SM_flag, format='E')
	c19 = fits.Column(name='flag_calib',array=extrapolation_mag, format='3A')
	c20 = fits.Column(name='Filter',array=np.array([name_filter]*len(final_objects[:,0])), format='10A')
	c21 = fits.Column(name='Elongation',array=np.around(final_objects[:,ELONGATION],2), format='E')
	c22 = fits.Column(name='Ellipticity',array=np.around(final_objects[:,ELLIPTICITY],2), format='E')
	c23 = fits.Column(name='FWHM', unit='arcsec',array=np.around(3600*final_objects[:,FWHM_WORLD],2), format='E')
	c24 = fits.Column(name='FLAGS',array=final_objects[:,FLAGS], format='E')
	c25 = fits.Column(name='FLAGS_WEIGHT',array=final_objects[:,FLAGS_WEIGHT], format='E')
	c26 = fits.Column(name='XWIN_IMAGE',array=final_objects[:,XWIN_IMAGE], format='E')
	c27 = fits.Column(name='YWIN_IMAGE',array=final_objects[:,YWIN_IMAGE], format='E')
	c28 = fits.Column(name='XMIN_IMAGE',array=final_objects[:,XMIN_IMAGE], format='E')
	c29 = fits.Column(name='YMIN_IMAGE',array=final_objects[:,YMIN_IMAGE], format='E')
	c30 = fits.Column(name='XMAX_IMAGE',array=final_objects[:,XMAX_IMAGE], format='E')
	c31 = fits.Column(name='YMAX_IMAGE',array=final_objects[:,YMAX_IMAGE], format='E')
	c32 = fits.Column(name='X_IMAGE',array=final_objects[:,X_IMAGE], format='E')
	c33 = fits.Column(name='Y_IMAGE',array=final_objects[:,Y_IMAGE], format='E')
	c34 = fits.Column(name='FLUX_AUTO',array=final_objects[:,FLUX_AUTO], format='E')
	c35 = fits.Column(name='FLUXERR_AUTO',array=final_objects[:,FLUXERR_AUTO], format='E')
	c36 = fits.Column(name='MAG_AUTO',array=final_objects[:,MAG_AUTO], format='E')
	c37 = fits.Column(name='MAGERR_AUTO',array=final_objects[:,MAGERR_AUTO], format='E')
	c38 = fits.Column(name='FLUX_MAX',array=final_objects[:,FLUX_MAX], format='E')
	c39 = fits.Column(name='FLUX_PSF',array=final_objects[:,FLUX_PSF], format='E')

	t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39],name='catalog')
	t.writeto('catalogs_folder/CFC_catalogs/'+fichero[0:len(fichero)-5]+'_catalog.fits',overwrite=True)
	votable1=Table.read('catalogs_folder/CFC_catalogs/'+fichero[0:len(fichero)-5]+'_catalog.fits')
	votable1.write('catalogs_folder/CFC_catalogs/'+fichero[0:len(fichero)-5]+'_catalog.xml',table_id='table_id',format='votable',overwrite=True)



	if CAHA_ID=='X':
		c1 = fits.Column(name='Image_identifier', array=np.array(len(total_objects[:,0])*['CAHA_CAFOS_BBI_DR1']), format='50A')
	else:
		c1 = fits.Column(name='Image_identifier', array=np.array(len(total_objects[:,0])*['CAHA_CAFOS_BBI_DR1_'+str(CAHA_ID[0])]), format='50A')
	DETECTION_ID=np.array(len(total_objects)*['CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_0000'])
	contador_A=1+len(final_objects)
	contador_B=1
	if CAHA_ID=='X':
		for j in range(len(total_objects)):
			if listaok[j]!=1:
				DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+'0'*(3-int(np.log10(contador_A)))+str(contador_A)
				contador_A=contador_A+1
			if listaok[j]==1:
				DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+'0'*(3-int(np.log10(contador_B)))+str(contador_B)
				contador_B=contador_B+1	
	else:
		for j in range(len(total_objects)):
			if listaok[j]!=1:
				DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_'+'0'*(3-int(np.log10(contador_A)))+str(contador_A)
				contador_A=contador_A+1
			if listaok[j]==1:
				DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_'+'0'*(3-int(np.log10(contador_B)))+str(contador_B)
				contador_B=contador_B+1
	c2 = fits.Column(name='Detection_ID', array=DETECTION_ID, format='50A')
	c3 = fits.Column(name='MJD', array=MJD_array, format='D')
	c4 = fits.Column(name='SNR_WIN',array=total_objects[:,SNR_WIN], format='E')
	#c5 = fits.Column(name='RAJ2000', unit='deg',array=np.around(total_objects[:,ALPHA_J2000],5), format='E')
	#c6 = fits.Column(name='DEJ2000', unit='deg',array=np.around(total_objects[:,DELTA_J2000],5), format='E')
	#c7 = fits.Column(name='e_RAJ2000', unit='arcsec',array=np.around(3600*total_objects[:,ERRX2_WORLD]**0.5,5), format='E')
	#c8 = fits.Column(name='e_DEJ2000', unit='arcsec',array=np.around(3600*total_objects[:,ERRY2_WORLD]**0.5,5), format='E')
	c5 = fits.Column(name='RAJ2000', unit='deg',array=total_objects[:,ALPHA_J2000], format='E')
	c6 = fits.Column(name='DEJ2000', unit='deg',array=total_objects[:,DELTA_J2000], format='E')
	c7 = fits.Column(name='e_RAJ2000', unit='arcsec',array=3600*total_objects[:,ERRX2_WORLD]**0.5, format='E')
	c8 = fits.Column(name='e_DEJ2000', unit='arcsec',array=3600*total_objects[:,ERRY2_WORLD]**0.5, format='E')
	RA=Angle(total_objects[:,ALPHA_J2000]* u.deg)
	DEC=Angle(total_objects[:,DELTA_J2000]* u.deg)
	e_RA=Angle(total_objects[:,ERRX2_WORLD]**0.5* u.deg)
	e_DEC=Angle(total_objects[:,ERRY2_WORLD]**0.5* u.deg)
	c9 = fits.Column(name='RA_hms', unit='hh:mm:ss', array=RA.to_string(unit=u.hourangle, sep=(':',':')), format='20A')
	c10 = fits.Column(name='DE_dms', unit='dd:mm:ss', array=DEC.to_string(unit=u.deg, sep=(':',':')), format='20A')
	c11 = fits.Column(name='e_RA_hms', unit='hh:mm:ss', array=e_RA.to_string(unit=u.hourangle, sep=(':',':')), format='20A')
	c12 = fits.Column(name='e_DE_dms', unit='dd:mm:ss', array=e_DEC.to_string(unit=u.deg, sep=(':',':')), format='20A')
	c13 = fits.Column(name='MAG',array=np.around(Z[1]+Z[0]*total_objects[:,17],3), format='E')
	c14 = fits.Column(name='e_MAG',array=np.around((abs(e_B)**2+abs(e_A*total_objects[:,17])**2+abs(Z[0]*total_objects[:,18])**2)**0.5,3), format='E')
	c15 = fits.Column(name='MAG_sex',array=np.around(total_objects[:,MAG_PSF],3), format='E')
	c16 = fits.Column(name='e_MAG_sex',array=np.around(total_objects[:,MAGERR_PSF],3), format='E')
	c17 = fits.Column(name='SPREAD_MODEL',array=total_objects[:,SPREAD_MODEL], format='E')
	c18 = fits.Column(name='Filter',array=np.array([name_filter]*len(total_objects[:,0])), format='10A')
	c19 = fits.Column(name='Elongation',array=np.around(total_objects[:,ELONGATION],2), format='E')
	c20 = fits.Column(name='Ellipticity',array=np.around(total_objects[:,ELLIPTICITY],2), format='E')
	c21 = fits.Column(name='FWHM', unit='arcsec',array=np.around(3600*total_objects[:,FWHM_WORLD],2), format='E')
	c22 = fits.Column(name='source_type',array=source_flag, format='E')
	c23 = fits.Column(name='FLAGS',array=total_objects[:,FLAGS], format='E')
	c24 = fits.Column(name='FLAGS_WEIGHT',array=total_objects[:,FLAGS_WEIGHT], format='E')
	c25 = fits.Column(name='XWIN_IMAGE',array=total_objects[:,XWIN_IMAGE], format='E')
	c26 = fits.Column(name='YWIN_IMAGE',array=total_objects[:,YWIN_IMAGE], format='E')
	c27 = fits.Column(name='XMIN_IMAGE',array=total_objects[:,XMIN_IMAGE], format='E')
	c28 = fits.Column(name='YMIN_IMAGE',array=total_objects[:,YMIN_IMAGE], format='E')
	c29 = fits.Column(name='XMAX_IMAGE',array=total_objects[:,XMAX_IMAGE], format='E')
	c30 = fits.Column(name='YMAX_IMAGE',array=total_objects[:,YMAX_IMAGE], format='E')
	c31 = fits.Column(name='X_IMAGE',array=total_objects[:,X_IMAGE], format='E')
	c32 = fits.Column(name='Y_IMAGE',array=total_objects[:,Y_IMAGE], format='E')
	c33 = fits.Column(name='FLUX_AUTO',array=total_objects[:,FLUX_AUTO], format='E')
	c34 = fits.Column(name='FLUXERR_AUTO',array=total_objects[:,FLUXERR_AUTO], format='E')
	c35 = fits.Column(name='MAG_AUTO',array=total_objects[:,MAG_AUTO], format='E')
	c36 = fits.Column(name='MAGERR_AUTO',array=total_objects[:,MAGERR_AUTO], format='E')
	c37 = fits.Column(name='FLUX_MAX',array=total_objects[:,FLUX_MAX], format='E')
	c38 = fits.Column(name='FLUX_PSF',array=total_objects[:,FLUX_PSF], format='E')

	t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38],name='catalog')
	t.writeto('catalogs_folder/CFC_sources/'+fichero[0:len(fichero)-5]+'_sources.fits',overwrite=True)
	votable2=Table.read('catalogs_folder/CFC_sources/'+fichero[0:len(fichero)-5]+'_sources.fits')
	votable2.write('catalogs_folder/CFC_sources/'+fichero[0:len(fichero)-5]+'_sources.xml',table_id='table_id',format='votable',overwrite=True)
if Z[4]<0.98:
	print('r pearson coefficient lower than 0.98')
