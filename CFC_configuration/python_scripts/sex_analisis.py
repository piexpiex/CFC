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
	hdulist = fits.open(fichero)
	fichero=delete_folder_name(fichero)
	name_filter=filter_id[np.where(fich_reducido_id==fichero+'\n')]#hdulist[0].header['INSFLNAM']
	name_filter=name_filter[0]
	CAHA_ID=caha_id[np.where(fich_reducido_id==fichero+'\n')]

except:
	hdulist = fits.open(fichero)
	name_filter=hdulist[0].header['INSFLNAM']
	CAHA_ID='X'
MJD=hdulist[0].header['MJD-OBS']
color=search_name(name_filter)

data = hdulist[0].data
if color=='X':
	print('no SDSS filter')
	images_table=open('logouts_folder/data_table.csv','a')
	images_table.write(fichero[0:len(fichero)-5]+','+'   ----No SDSS filter---')
	exit()
#hdulist = fits.open("CFC_configuration/sextractor_result_files/bkg.fits")
#data_bkg= hdulist[0].data
#data_sub=data-data_bkg
# show the image

########################
### Reading test.cat ###
########################

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
X_WORLD=23##  24 X_WORLD                Barycenter position along world x axis                     [deg]
Y_WORLD=24#  25 Y_WORLD                Barycenter position along world y axis                     [deg]
ERRX2_WORLD=25#  26 ERRX2_WORLD            Variance of position along X-WORLD (alpha)                 [deg**2]
ERRY2_WORLD=26#  27 ERRY2_WORLD            Variance of position along Y-WORLD (delta)                 [deg**2]


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


Nobjetos=0

listaok=np.array([1.0]*len(objects))

for i in range(len(objects)):
	if objects[i,FLUX_MAX]<=0 or objects[i,FLUX_RADIUS]<=0 or objects[i,FWHM_IMAGE]<=0 or objects[i,FLUX_PSF]<=0:
		listaok[i]=0
	if objects[i,MAG_PSF]==99.0:
		listaok[i]=0
	if objects[i,FLAGS_WEIGHT]==2:
		listaok[i]=0
	if abs(objects[i,SNR_WIN])<5 or abs(objects[i,SNR_WIN])>0.99*10**30:
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
FLUX_SATURATION=min([FLUX_SATURATION_1,FLUX_SATURATION_2])-2*flux_sigma2

FLUX_step=900
FLUX_step_d=10
FLUX_DETECTION_1=0
while abs(FLUX_step_d)>1 and FLUX_step>0: # Pearson r criteria
	FLUX_step=FLUX_step-FLUX_step_d
	r=ajuste_lineal(psf_flux[np.where((max_flux<FLUX_SATURATION) & (max_flux>FLUX_step))],max_flux[np.where((max_flux<FLUX_SATURATION) & (max_flux>FLUX_step))])[4]
	if FLUX_step<0:
		FLUX_step=0
		break
	if r>=0.98:
		continue
	if r<0.98:
		FLUX_DETECTION_1=FLUX_step
		FLUX_step=FLUX_step+FLUX_step_d
		FLUX_step_d=FLUX_step_d/2

FLUX_DETECTION_1=FLUX_step+FLUX_step_d	
FLUX_step=900
FLUX_step_d=10
FLUX_DETECTION_2=0
A0=ajuste_lineal(psf_flux[np.where((max_flux<FLUX_SATURATION) & (max_flux>FLUX_step))],max_flux[np.where((max_flux<FLUX_SATURATION) & (max_flux>FLUX_step))])[0]
while abs(FLUX_step_d)>1 and FLUX_step>0: # Slope criteria
	FLUX_step=FLUX_step-FLUX_step_d
	A=ajuste_lineal(psf_flux[np.where((max_flux<FLUX_SATURATION) & (max_flux>FLUX_step))],max_flux[np.where((max_flux<FLUX_SATURATION) & (max_flux>FLUX_step))])[0]
	if FLUX_step<0:
		FLUX_step=0
		break
	if A>=A0/1.01:
		continue
	if A<A0/1.01:
		FLUX_DETECTION_2=FLUX_step
		FLUX_step=FLUX_step+FLUX_step_d
		FLUX_step_d=FLUX_step_d/2

FLUX_DETECTION_2=FLUX_step+FLUX_step_d	
FLUX_DETECTION=max([FLUX_DETECTION_1,FLUX_DETECTION_2])+2*flux_sigma2

for i in range(len(objects)):
	if objects[i][FLUX_MAX]>FLUX_SATURATION or objects[i][FLUX_MAX]<FLUX_DETECTION:
		listaok[i]=2

plt.figure(figsize=(22.0,7.0))

plt.plot(objects[np.where(listaok==2),FLUX_PSF], objects[np.where(listaok==2),FLUX_MAX],'r.')
plt.plot(objects[np.where(listaok==1),FLUX_PSF], objects[np.where(listaok==1),FLUX_MAX],'g.')
plt.plot(objects[np.where(listaok==-1),FLUX_PSF], objects[np.where(listaok==-1),FLUX_MAX],'b.')
plt.plot(objects[np.where(listaok==2),FLUX_PSF][0], objects[np.where(listaok==2),FLUX_MAX][0],'r.',label='saturated')
plt.plot(objects[np.where(listaok==1),FLUX_PSF][0], objects[np.where(listaok==1),FLUX_MAX][0],'g.',label='Valid')
plt.plot(objects[np.where(listaok==-1),FLUX_PSF][0], objects[np.where(listaok==-1),FLUX_MAX][0],'b.',label='artifacts')

#plt.plot([1*10**4,10**5,max(f1)],ffit([1*10**4,10**5,max(f1)]),'k')
plt.suptitle('FLUX MAX vs FLUX PSF (SNR >=5)')
plt.ylabel('FLUX MAX [count]')
plt.xlabel('FLUX PSF [count]')
plt.yscale('log')
plt.xscale('log')
#plt.legend(framealpha=0.1)

plt.savefig('figures_folder/'+'_flux_selection.pdf')
#plt.savefig('figures_folder/'+fichero[0:len(fichero)-5]+'_flux_selection.pdf')


############################
######### Skymatch #########
############################

objects=objects[np.where(listaok==1)]

SM_flag=objects[:,SPREAD_MODEL]
#for i in range(len(objects)):
#	if objects[i,SPREAD_MODEL]<-0.05:
#		SM_flag[i]=0.0
cl_SDSS=np.array([0.0]*len(objects[:,SPREAD_MODEL]))
final_objects=objects
objects=objects[np.where(SM_flag>-0.05)]

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
	exit()	



#Columns selection
NUMBER_XMATCH=catalog['NUMBER']
mag_sex=catalog['MAG_PSF']
magerr_sex=catalog['MAGERR_PSF']
ellongation=catalog['ELONGATION']
ellipticity=catalog['ELLIPTICITY']
FWHM=catalog['FWHM_WORLD']
SPREAD_VALUE=catalog['SPREAD_MODEL']
limit_detection=25
limit_sat=0
if sdss_key==0:
	q_mode=catalog['q_mode']
	class_sdss=catalog['class']
	for j in range(len(class_sdss)):
		cl_SDSS[np.where(final_objects[:,0]==NUMBER_XMATCH[j])]=class_sdss[j]
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

for k in range(len(lista)):
	lista[k]=lista[k][np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]
if sdss_key==0:
	class_sdss=class_sdss[np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]
	q_mode=q_mode[np.where((mag_sex<0.2) & (np.isnan(pmag)==False) & (pmag<limit_detection) & (pmag>limit_sat) & (e_pmag<0.2))]

mag_sex=lista[0]
magerr_sex=lista[1]
ellongation=lista[2]
ellipticity=lista[3]
FWHM=lista[4]
pmag=lista[5]
e_pmag=lista[6]
SPREAD_VALUE=lista[7]

#Morphology selection

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



print('number of objects after morphology criteria',len(mag_sex))

N_B=len(pmag)
plt.plot(mag_sex,pmag,'b.',label='objects with morphology criteria ('+str(N_B)+')')

if sdss_key==0:
	for k in range(len(lista)):
		try:
			lista[k]=lista[k][np.where((class_sdss==6) & (q_mode=='+'))]
		except:
			lista[k]=lista[k][np.where((class_sdss==6) & (q_mode==1.0))]

	mag_sex=lista[0]
	magerr_sex=lista[1]
	pmag=lista[5]
	e_pmag=lista[6]
	SPREAD_VALUE=lista[7]
	print('number of objects with SDSS class=6 & q mode=1.0 ',len(mag_sex))
	

if len(mag_sex)<6:
	print('Insufficient number of objects for photometric calibration')
	exit()
	
#if abs(max(mag_sex)-min(mag_sex))<2:
#	print('Small magnitude range,photometric calibration may not be correct')

#Photometry calibration


X,Y,Z=sigma_c(X=mag_sex,Y=pmag,n_sigma=2)
N_C=len(X)
plt.plot(X,Y,'r.',label='Calibration objects ('+str(N_C)+')')
plt.plot(mag_sex,Z[1]+Z[0]*mag_sex,'r',label=name_mag+'={0:.3g}'.format(Z[0])+' × MAG PSF + {0:.3g}'.format(Z[1]) + ' (r={0:.3g}'.format(Z[4])+')')
plt.legend()
#plt.ylim(limit_sat,limit_detection+1)

plt.savefig('figures_folder/'+fichero[0:len(fichero)-5]+'_magnitude_calibration.pdf')


###########################################
### Comparasion between both magnitudes ###
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
plt.savefig('figures_folder/'+fichero[0:len(fichero)-5]+'_magnitude_comparation.pdf')



##############
### logout ###
##############

state='rejected'
if Z[4]>=0.98:
	state='calibrated'
images_table=open('logouts_folder/data_table.csv','a')

if sdss_key==0:
	images_table.write(fichero[0:len(fichero)-5] +','+  str(N_C) +','+  str(round(Z[1],3))+','+  str(round(Z[3],3)) +','+  str(round(Z[0],3))+','+  str(round(Z[2],3)) +','+ str(round(Z[4],3)) +','+  state+'\n')
if sdss_key==1:
	images_table.write(fichero[0:len(fichero)-5] +','+  str(N_B) +','+  str(round(Z[1],3))+','+  str(round(Z[3],3)) +','+  str(round(Z[0],3))+','+  str(round(Z[2],3)) +','+ str(round(Z[4],3)) +','+  state+'\n')

images_table.close
###########
if Z[4]>=0.98:
	calibration_mag=Z[1]+Z[0]*final_objects[:,17]
	calibration_mag_error=abs(Z[3])+abs(Z[2]*final_objects[:,17])+abs(Z[0]*final_objects[:,18])

	min_pmag=min(Y)
	max_pmag=max(Y)

	extrapolation_mag=np.array(['A']*len(final_objects))
	for i in range(len(final_objects)):
		if calibration_mag[i]>max_pmag:
			extrapolation_mag[i]='C'
		elif calibration_mag[i]<min_pmag:
			extrapolation_mag[i]='B'
	c1 = fits.Column(name='Image_identifier', array=np.array(len(final_objects[:,0])*['CAHA_CAFOS_BBI_DR1_'+str(CAHA_ID[0])]), format='50A')
	NUMBER_ID=np.arange(1,1+len(final_objects[:,0]),1).astype(np.str)
	DETECTION_ID=np.array(len(NUMBER_ID)*['CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_0000'])
	for j in range(len(NUMBER_ID)):
		DETECTION_ID[j]='CAHA_CAFOS_BBI_DR1_'+CAHA_ID[0]+'_'+'0'*(3-int(np.log10(1+j)))+NUMBER_ID[j]
	c2 = fits.Column(name='Detection_ID', array=DETECTION_ID, format='50A')
	c3 = fits.Column(name='MJD', array=np.array(len(final_objects[:,0])*[str(MJD)]), format='12A')
	c4 = fits.Column(name='SNR_WIN',array=final_objects[:,SNR_WIN], format='E')
	c5 = fits.Column(name='RAJ2000 (deg)',array=np.around(final_objects[:,ALPHA_J2000],5), format='E')
	c6 = fits.Column(name='DEJ2000 (deg)',array=np.around(final_objects[:,DELTA_J2000],5), format='E')
	c7 = fits.Column(name='e_RAJ2000 (arcsec)',array=np.around(3600*final_objects[:,ERRX2_WORLD]**0.5,5), format='E')
	c8 = fits.Column(name='e_DEJ2000 (arcsec)',array=np.around(3600*final_objects[:,ERRY2_WORLD]**0.5,5), format='E')
	RA=Angle(final_objects[:,ALPHA_J2000]* u.deg)
	DEC=Angle(final_objects[:,DELTA_J2000]* u.deg)
	e_RA=Angle(final_objects[:,ERRX2_WORLD]**0.5* u.deg)
	e_DEC=Angle(final_objects[:,ERRY2_WORLD]**0.5* u.deg)
	c9 = fits.Column(name='RAJ2000 (hh:mm:ss)', array=RA.to_string(unit=u.hourangle, sep=(':',':')), format='20A')
	c10 = fits.Column(name='DEJ2000 (dd:mm:ss)', array=DEC.to_string(unit=u.deg, sep=(':',':')), format='20A')
	c11 = fits.Column(name='e_RAJ2000 (hh:mm:ss)', array=e_RA.to_string(unit=u.hourangle, sep=(':',':')), format='20A')
	c12 = fits.Column(name='e_DEJ2000 (dd:mm:ss)', array=e_DEC.to_string(unit=u.deg, sep=(':',':')), format='20A')
	c13 = fits.Column(name='MAG',array=np.around(calibration_mag,3), format='E')
	c14 = fits.Column(name='e_MAG',array=np.around(calibration_mag_error,3), format='E')
	c15 = fits.Column(name='MAG_sex',array=np.around(final_objects[:,MAG_PSF],3), format='E')
	c16 = fits.Column(name='e_MAG_sex',array=np.around(final_objects[:,MAGERR_PSF],3), format='E')
	c17 = fits.Column(name='cl_SDSS',array=cl_SDSS, format='E')
	c18 = fits.Column(name='SPREAD_MODEL',array=np.around(SM_flag,2), format='E')
	c19 = fits.Column(name='flag_calib',array=extrapolation_mag, format='3A')
	c20 = fits.Column(name='Filter',array=np.array([name_filter]*len(final_objects[:,0])), format='10A')
	c21 = fits.Column(name='Elongation',array=np.around(final_objects[:,ELONGATION],2), format='E')
	c22 = fits.Column(name='FWHM (arcsec)',array=np.around(3600*final_objects[:,FWHM_WORLD],2), format='E')
	
	t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22],name='catalog')
	t.writeto('catalogs_folder/'+fichero[0:len(fichero)-5]+'_catalog.fits',overwrite=True)


