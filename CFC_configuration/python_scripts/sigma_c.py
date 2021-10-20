import matplotlib.pyplot as plt
import numpy as np

def ajuste_lineal(X,Y,W=0):
	W=np.array([W])
	#ajuste lineal por minimos cuadrados de los puntos con coordenadas X e Y con pesos W
	if len(W)==1:
		W=W[0]
		n=len(X)
		x=sum(X)
		y=sum(Y)
		xy=sum(X*Y)
		x_2=sum(X*X)
		div=n*x_2-x*x
		if div==0:
			A=0
			B=0
			d_A=0
			d_B=0
			r=0
			r_2=0
			R_2=0
			sigma=0
		else:
			A=(n*xy-x*y)/div
			#A=(n*xy-x*y)/(n*x_2-x*x)
			B=(y-A*x)/n
			#B=(x_2*y-x*xy)/(n*x_2-x*x) #desarrollado
			rest=(Y-B-A*X)
			sigma=(sum((Y-A*X-B)**2)/(n-2))**0.5
			d_A=sigma*(n/div)**0.5
			#d_A=sigma*(n/(n*x_2-x**2))**0.5
			d_B=d_A*(x_2/n)**0.5
			r=sum((X-np.mean(X))*(Y-np.mean(Y)))/sum((X-np.mean(X))**2)**0.5/sum((Y-np.mean(Y))**2)**0.5
			r_2=r**2
			R_2=1- (sum((rest-n*np.mean(rest))**2)/(n-1))/(sum((Y-n*np.mean(Y))**2)/(n-1))
	else:
		n=len(X)
		N=sum(W)
		x=sum(X*W)
		y=sum(Y*W)
		xy=sum(W*X*Y)
		x_2=sum(W*X*X)
		A=(N*xy-x*y)/(N*x_2-x*x)
		B=(y-A*x)/N
		rest=W*(Y-B-A*X)
		sigma=(sum((Y-A*X-B)**2)/(n-2))**0.5
		d_A=sigma*(N/(N*x_2-x**2))**0.5
		d_B=d_A*(x_2/N)**0.5
		r=sum((W*X-N*np.mean(X))*(W*Y-N*np.mean(Y)))/sum((W*X-N*np.mean(X))**2)**0.5/sum((W*Y-N*np.mean(Y))**2)**0.5
		r_2=r**2
		R_2=1- (sum((rest-N*np.mean(rest))**2)/(n-1))/(sum((W*Y-N*np.mean(Y))**2)/(n-1))
	return(A,B,d_A,d_B,r,r_2,R_2,sigma)
	
def sigma_c(X,Y,idfs,W=0,n_sigma=1,n_iteraciones=1):
	x=X
	y=Y
	identificadores=idfs #orden en el array
	ajuste=ajuste_lineal(x,y)
	for k in range(n_iteraciones):	
		A=ajuste[0]
		B=ajuste[1]
		
		lista=np.where((y<ajuste[1]+ajuste[0]*x+n_sigma*ajuste[7]) & (y>ajuste[1]+ajuste[0]*x-n_sigma*ajuste[7]))
		y_k=y[lista]
		x_k=x[lista]
		NUMBER_XMATCH_k=identificadores[lista]
		ajuste=ajuste_lineal(x_k,y_k)
	return(x_k,y_k,NUMBER_XMATCH_k,ajuste,lista)

