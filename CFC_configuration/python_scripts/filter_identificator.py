# -*-coding: utf-8 -*-

from read_files import *

def search_name(filter_name):
	for j in range(len(filter_name)-4):
		if filter_name[j:j+4]=='sdss' or filter_name[j:j+4]=='SDSS':
			filter_name=filter_name[j+4:]
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
	return(filter_name)

