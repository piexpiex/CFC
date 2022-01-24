# -*-coding: utf-8 -*-

from read_files import *

def search_name(name_filter):
	filter_name=name_filter
	for j in range(len(filter_name)-4):
		if filter_name[j:j+4]=='sdss' or filter_name[j:j+4]=='SDSS' or filter_name[j:j+4]=='Sdss':
			filter_name=filter_name[j+4:]
	for j in range(len(filter_name)-5):
		if filter_name[j:j+5]=='SLOAN' or filter_name[j:j+5]=='Sloan' or filter_name[j:j+5]=='sloan':
			filter_name=filter_name[j+5:]
			
	filter_name=delete_space(filter_name)
	filter_name=delete_space2(filter_name)
	
	if filter_name=='u' or filter_name=='U':
		filter_name='u'
	elif filter_name=='g' or filter_name=='G':
		filter_name='g'
	elif filter_name=='r' or filter_name=='R':
		filter_name='r'
	elif filter_name=='i' or filter_name=='I':
		filter_name='i'
	elif filter_name=='z' or filter_name=='Z':
		filter_name='z'
	else:
		filter_name='X'
		
	option_A=filter_name
		
	filter_name=name_filter
	
	for j in range(len(filter_name)-4):
		if filter_name[j:j+4]=='sdss' or filter_name[j:j+4]=='SDSS' or filter_name[j:j+4]=='Sdss':
			filter_name=filter_name[:j]
	for j in range(len(filter_name)-5):
		if filter_name[j:j+5]=='SLOAN' or filter_name[j:j+5]=='Sloan' or filter_name[j:j+5]=='sloan':
			filter_name=filter_name[:j]
	
	filter_name=delete_space(filter_name)
	
	if filter_name=='u' or filter_name=='U':
		filter_name='u'
	elif filter_name=='g' or filter_name=='G':
		filter_name='g'
	elif filter_name=='r' or filter_name=='R':
		filter_name='r'
	elif filter_name=='i' or filter_name=='I':
		filter_name='i'
	elif filter_name=='z' or filter_name=='Z':
		filter_name='z'
	else:
		filter_name='X'
	option_B=filter_name
	
	if option_A=='X':
		filter_name=option_B
	else:
		filter_name=option_A
	return(filter_name)
