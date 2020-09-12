import numpy as np

def find_sources(X,Y,X2,Y2):
	Numbers=np.array([0]*len(X))
	value= (np.transpose([X] * len(X2)) - np.tile(X2, (len(X), 1)))**2 + (np.transpose([Y] * len(Y2)) - np.tile(Y2, (len(Y), 1)))**2

	for j in range(len(X)):
		Numbers[j]=np.where(value[j]==min(value[j]))[0][0]

	return(Numbers)
