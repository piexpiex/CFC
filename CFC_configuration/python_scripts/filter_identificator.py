

def search_name(name):
	filter_name=''
	for j in range(len(name)):
		if name[j]==' ':
			continue
		elif j==0:
			if name[j+1]==' ':
				filter_name=name[j]
		elif j==len(name)-1:
			if name[len(name)-2]==' ':
				filter_name=name[j]
		else:
			if name[j-1]==' ' and name[j+1]==' ':
				filter_name=name[j]		
			else:
				continue
	if filter_name=='' or filter_name==' ':
		filter_name='r'
	if filter_name=='u' or filter_name=='U':
		filter_name='u'
	if filter_name=='g' or filter_name=='G':
		filter_name='g'
	if filter_name=='r' or filter_name=='R':
		filter_name='r'
	if filter_name=='i' or filter_name=='I':
		filter_name='i'
	if filter_name=='z' or filter_name=='Z':
		filter_name='z'
		
	return(filter_name)

