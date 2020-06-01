import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle
from matplotlib import rcParams
from math import *
from sigma_c import *
from filter_identificator import *
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
#from astropy.samp import SAMPIntegratedClient
from astropy.wcs import WCS
#from astropy.stats import sigma_clipped_stats
import astropy.visualization as ap_vis
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astroquery.simbad import Simbad
import numpy as np
from matplotlib.offsetbox import AnchoredText
import os
import sys
import warnings
warnings.filterwarnings('ignore')

########################
### Folders creation ###
########################

try:
	os.mkdir('logouts_folder')
	os.mkdir('figures_folder')
	os.mkdir('catalogs_folder')
	print('making directories')
	print('starting the analysis')
except:
	print('starting the analysis')

#############################
## Lectura de las imagenes ##
#############################

fichero=sys.argv[1]
hdulist = fits.open(fichero)
name_filter=hdulist[0].header['INSFLNAM']
color=search_name(name_filter)
data = hdulist[0].data

#hdulist = fits.open("background_substracted.fits") #caf-20170225-21_51_59-sci-krek.fits
hdulist = fits.open("CFC_configuration/sextractor_result_files/bkg.fits")
data_bkg= hdulist[0].data


data_sub=data-data_bkg
# show the image

#leyendo test.cat
catalogo = open('CFC_configuration/sextractor_result_files/test.cat')
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

lista=[]
objects=[]

hj=0
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
		hj=hj+1
objects=np.array(objects)

#Object detection
# plot background-subtracted image
fig, ax = plt.subplots()
#m, s = np.mean(data_sub), np.std(data_sub)
m, s = np.mean(data_sub), np.std(data_sub)
im = ax.imshow(data_sub, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
Nobjetos=0
# plot an ellipse for each object
listaok=np.array([1.0]*len(objects))

#print('objetos al principio',len(objects))

for i in range(len(objects)):
	if objects[i,FLUX_MAX]<=0 or objects[i,FLUX_RADIUS]<=0 or objects[i,FWHM_IMAGE]<=0:
		listaok[i]=0
	if objects[i,MAG_PSF]==99.0:
		listaok[i]=0
	if objects[i,FLAGS_WEIGHT]==2:
		listaok[i]=0
	if abs(objects[i,SNR_WIN])<5 or abs(objects[i,SNR_WIN])>0.99*10**30:
		listaok[i]=0
	if np.shape(data_sub)==(1650,1650):
		if ((abs(objects[i,X_IMAGE]-825)+objects[i,A_IMAGE])**2+ (abs(objects[i,Y_IMAGE]-825)+objects[i,A_IMAGE])**2)**0.5>814: 
			listaok[i]=0
	if np.shape(data_sub)==(1601,1601):
		if ((abs(objects[i,X_IMAGE]-808)+objects[i,A_IMAGE])**2+ (abs(objects[i,Y_IMAGE]-817)+objects[i,A_IMAGE])**2)**0.5>814:
			listaok[i]=0
	if listaok[i]==1.0:
		Nobjetos=Nobjetos+1
		listaok[i]=1.0
		

print('number of identify objects',len(listaok[np.where(listaok==1.0)]))

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

for i in range(len(objects)):
	if objects[i][FLUX_MAX]>0.9*max(objects[:,FLUX_MAX]):
		listaok[i]=2

plt.figure(figsize=(22.0,7.0))
for i in range(len(objects)):
	if listaok[i]==2:
		plt.plot(objects[i,FLUX_PSF], objects[i,FLUX_MAX],'r.')
	elif listaok[i]==1:
		plt.plot(objects[i,FLUX_PSF], objects[i,FLUX_MAX],'g.')
	elif listaok[i]==-1:
		plt.plot(objects[i,FLUX_PSF], objects[i,FLUX_MAX],'b.')
#plt.plot([1*10**4,10**5,max(f1)],ffit([1*10**4,10**5,max(f1)]),'k')
plt.suptitle('FLUX MAX vs FLUX PSF (SNR >=5)')
plt.ylabel('FLUX MAX [count]')
plt.xlabel('FLUX PSF [count]')
plt.yscale('log')
plt.xscale('log')
plt.savefig('figures_folder/'+fichero[6:len(fichero)-5]+'flux_selection.pdf')


###########################################
######### Calibración fotometrica #########
###########################################

objects=objects[np.where(listaok==1)]

PSF_FIT=np.array([1.0]*len(objects))
for i in range(len(objects)):
	if objects[i,SPREAD_MODEL]<-0.05:
		PSF_FIT[i]=0.0

final_objects=objects
objects=objects[np.where(PSF_FIT==1.0)]

c1 = fits.Column(name='NUMBER', array=objects[:,0], format='E')
c2 = fits.Column(name='SNR_WIN',array=objects[:,3], format='E')
c3 = fits.Column(name='FLUX_MAX',array=objects[:,10], format='E')
c4 = fits.Column(name='ALPHA',array=objects[:,12], format='E')
c5 = fits.Column(name='DELTA',array=objects[:,13], format='E')
c6 = fits.Column(name='FWHM_WORLD',array=objects[:,14], format='E')
c7 = fits.Column(name='ELLONGATION',array=objects[:,19], format='E')
c8 = fits.Column(name='ELLIPTICITY',array=objects[:,20], format='E')
c9 = fits.Column(name='FLUX_PSF',array=objects[:,15], format='E')
c10 = fits.Column(name='FLUXERR_PSF',array=objects[:,16], format='E')
c11 = fits.Column(name='MAG_PSF',array=objects[:,17], format='E')
c12 = fits.Column(name='MAGERR_PSF',array=objects[:,18], format='E')
c13 = fits.Column(name='SPREAD_MODEL',array=objects[:,22], format='E')

#t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8],name='valores')
#t.writeto('a.fits',overwrite=True)

t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13],name='valores')
t.writeto('CFC_configuration/Intermediate_files/parametros.fits',overwrite=True)

#objects[i,MAG_ISO] #magnitud sextractor
#INSFLID #numero del filtro
#INSFLNAM #nombre del filtro (en este caso SDSS r)
# magnitudes de APASS V,B,g,r,i seguido de mag o err
#rmag y rerr en este caso

hdul = fits.open('CFC_configuration/Intermediate_files/parametros.fits')
delta=hdul['VALORES']
t = Table.read(delta)
t.write('CFC_configuration/Intermediate_files/parametros.vot', table_id='updated_table', format='votable',overwrite=True)

#cross match (en open nuestra tabla)
#V/147/sdss12 #II/336/apass9 #hay que ponerlo para que intente SDSS y si no APASS
sdss_key=0
try:
	catalog = XMatch.query(cat1=open('CFC_configuration/Intermediate_files/parametros.vot'),
                         cat2='vizier:V/147/sdss12',
                         max_distance=2 * u.arcsec,
                         colRA1='ALPHA', colDec1='DELTA')
	sdss_key=0
except:
	catalog = XMatch.query(cat1=open('CFC_configuration/Intermediate_files/parametros.vot'),
                         cat2='vizier:II/336/apass9',
                         max_distance=2 * u.arcsec,
                         colRA1='ALPHA', colDec1='DELTA')
	sdss_key=1
#seleccion de columnas







mag_sex=catalog['MAG_PSF']
magerr_sex=catalog['MAGERR_PSF']
ellongation=catalog['ELLONGATION']
ellipticity=catalog['ELLIPTICITY']
FWHM=catalog['FWHM_WORLD']
SPREAD_VALUE=catalog['SPREAD_MODEL']
limit_detection=25
limit_sat=0
if sdss_key==0:
	q_mode=catalog['q_mode']
	class_sdss=catalog['class']
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
	if color=='r'or color=='R':
		pmag=catalog['rpmag']
		e_pmag=catalog['e_rpmag']
		name_mag='APASS rmag'
		limit_detection=25
		limit_sat=0
	elif color=='v'or color=='V':
		pmag=catalog['Vmag']
		e_pmag=catalog['e_Vmag']
		name_mag='APASS Vmag'
		limit_detection=25
		limit_sat=0
	elif color=='b'or color=='B':
		pmag=catalog['Bmag']
		e_pmag=catalog['e_Bmag']
		name_mag='APASS Bmag'
		limit_detection=25
		limit_sat=0
	elif color=='g'or color=='G':
		pmag=catalog['gpmag']
		e_pmag=catalog['e_gpmag']
		name_mag='APASS gmag'
		limit_detection=25
		limit_sat=0
	elif color=='i'or color=='I':
		pmag=catalog['ipmag']
		e_pmag=catalog['e_ipmag']
		name_mag='APASS imag'
		limit_detection=25
		limit_sat=0

lista=[mag_sex,magerr_sex,ellongation,ellipticity,FWHM,pmag,e_pmag,SPREAD_VALUE]
print('number of skymatch objects',len(pmag))
N_A=len(pmag)


for k in range(len(lista)):
	lista[k]=lista[k][np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]
if sdss_key==0:
	class_sdss=class_sdss[np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]
	q_mode=q_mode[np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]



##########################################
## criterios de seleccion en fotometria ##
##########################################

mag_sex=lista[0]
magerr_sex=lista[1]
ellongation=lista[2]
ellipticity=lista[3]
FWHM=lista[4]
pmag=lista[5]
e_pmag=lista[6]
SPREAD_VALUE=lista[7]


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
if sdss_key==0:
	class_sdss=class_sdss[np.where((FWHM<median_FWHM+2*Rq_FWHM) & (FWHM>median_FWHM-2*Rq_FWHM)
	&(ellongation<median_ellongation+2*Rq_ellongation) & (ellongation>median_ellongation-2*Rq_ellongation)
	&(ellipticity<median_ellipticity+2*Rq_ellipticity) &(ellipticity>median_ellipticity-2*Rq_ellipticity))]
	q_mode=q_mode[np.where((FWHM<median_FWHM+2*Rq_FWHM) & (FWHM>median_FWHM-2*Rq_FWHM)
	&(ellongation<median_ellongation+2*Rq_ellongation) & (ellongation>median_ellongation-2*Rq_ellongation)
	&(ellipticity<median_ellipticity+2*Rq_ellipticity) &(ellipticity>median_ellipticity-2*Rq_ellipticity))]

mag_sex=lista[0]
magerr_sex=lista[1]
pmag=lista[5]
e_pmag=lista[6]
SPREAD_VALUE=lista[7]



print('number of objects after selection criteria',len(mag_sex))
N_B=len(pmag)
if sdss_key==0:
	for k in range(len(lista)):
		lista[k]=lista[k][np.where((class_sdss==6) & (q_mode==1.0))]
	mag_sex=lista[0]
	magerr_sex=lista[1]
	pmag=lista[5]
	e_pmag=lista[6]
	SPREAD_VALUE=lista[7]
	print('number of objects with SDSS class=6 & q mode=1.0 ',len(mag_sex))
	N_C=len(pmag)

if len(mag_sex)<6:
	print('Insufficient number of objects for photometric calibration')
	exit()
if abs(max(mag_sex)-min(mag_sex))<2:
	print('Small magnitude range,photometric calibration may not be correct')

#ajuste=ajuste_lineal(mag_sex,pmag)

plt.figure(figsize=(22.0,7.0))
if sdss_key==0:
	plt.suptitle('Ajuste de magnitudes (n='+str(N_A)+'/'+str(N_B)+'/'+str(N_C)+')')
if sdss_key==1:
	plt.suptitle('Ajuste de magnitudes (n='+str(N_A)+'/'+str(N_B)+')')
plt.xlabel('MAG PSF ('+name_filter+')')
plt.ylabel(name_mag)

plt.plot(mag_sex,pmag,'k.')

X,Y,Z=sigma_c(X=mag_sex,Y=pmag,n_sigma=3)
plt.plot(X,Y,'b.')
X,Y,Z=sigma_c(X=mag_sex,Y=pmag,n_sigma=2)
plt.plot(X,Y,'r.')
plt.plot(mag_sex,Z[1]+Z[0]*mag_sex,'r',label=name_mag+'={0:.3g}'.format(Z[0])+' × MAG SEX + {0:.3g}'.format(Z[1]) + ' (r={0:.3g}'.format(Z[4])+')')
plt.legend()
plt.ylim(limit_sat,limit_detection+1)

plt.savefig('figures_folder/'+fichero[6:len(fichero)-5]+'magnitude_calibration.pdf')

c1 = fits.Column(name='NUMBER', array=final_objects[:,0], format='E')
c2 = fits.Column(name='SNR_WIN',array=final_objects[:,3], format='E')
c3 = fits.Column(name='ALPHA',array=final_objects[:,12], format='E')
c4 = fits.Column(name='DELTA',array=final_objects[:,13], format='E')
c5 = fits.Column(name='MAG',array=Z[1]+Z[0]*final_objects[:,17], format='E')
c6 = fits.Column(name='MAG error',array=Z[3]+Z[2]*final_objects[:,17]+Z[0]*final_objects[:,18], format='E')
c7 = fits.Column(name='PSF fit',array=PSF_FIT, format='E')


t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7],name='catalog')
t.writeto('catalogs_folder/'+fichero[6:len(fichero)-5]+'catalog.fits',overwrite=True)

