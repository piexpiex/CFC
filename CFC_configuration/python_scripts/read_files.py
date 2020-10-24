def delete_space(A):
	while A[0]==' ' and A!=' ':
		A=A[1:len(A)]
	while A[len(A)-1]==' ' and A!=' ':
		A=A[0:len(A)-1]
	return(A)
def delete_folder_name(A):
	B=''
	key=0
	for j in range(len(A)):
		if A[j]=='/':
			key=j
	B=A[key+1:]
	return(B)
def add_space(A):
	word=''
	for j in range(len(A)):
		if A[j] !=' ':word=word+A[j]
		else:word=word+'_'
	return(word)
def lister(A):
	lista=[]
	word=''
	for j in range(len(A)):
		if A[j] ==',':
			lista.append(word)
			word=''
		else:
			word=word+A[j]
		if j==len(A)-1:
			lista.append(word)
			word=''
	return(lista)
def read_files(fichero,keyword):
	
	CMD = open(fichero)
	
	lista=[]
	run_file=[]

	for linea in CMD:
		if linea[0]==' ':
			continue
		else:
			medidor=0
			cuenta=0
			numero=''
			for k in range(len(linea)):
				if linea[0]==' ' or linea[0]=='#' or linea[0:3]=='rem':
					continue
				if linea[k]!='[' and linea[k]!=']' or medidor==2:
					numero=numero+linea[k]
				if linea[k]=='=':
					lista.append(numero)
					numero=''
					medidor=1
				elif linea[k:k+3]=='rem' or linea[k]=='#' or k==len(linea)-1:
					if medidor==1:
						lista.append(numero[0:len(numero)-1])
						numero=''
						run_file.append(lista)
						medidor=2
				else:
					continue
		
		lista=[]
	
	for j in range(len(run_file)):
		if run_file[j][0][0:len(keyword)]==keyword:
			palabra=lister(delete_space(run_file[j][1][len(keyword)+1,:]))
		
			
	
	return(palabra)

def write_files(fichero,word,value):
	
	file_to_write=[]
	CMD = open(fichero)
	lista=[]
	run_file=[]

	for linea in CMD:
		if linea[0:len(word)]==word:
			post_linea=''
			pre_linea=linea[0:len(linea)]
			for k in range(len(linea)-2):
				if linea[k]=='#':
					post_linea=linea[k:len(linea)]
					pre_linea=linea[0:k]
					break
			for j in range(len(pre_linea)):
				if pre_linea[j]==' ':
					pre_linea=pre_linea[0:j+1]+str(value)+' '
					break
			linea=pre_linea+post_linea #WEIGHT_IMAGE
		elif linea[0:12]=='WEIGHT_IMAGE':
			linea=''
		
		file_to_write.append([linea])
	CMD2 = open(fichero,'w')
	for j in range(len(file_to_write)):
		CMD2.write(file_to_write[j][0])
	return()
