
# import math

import numpy as np

def eckartFit(xdata,ydata):



	return None

def sqrt(x):
	return np.power(x,0.5)

def eckart(x,barrier,asymmetry,length,xts=0.0):

	if isinstance(x,list):
		x = np.array(x)

	reverse_barrier = barrier - asymmetry

	# print(reverse_barrier)

	A = barrier - reverse_barrier

	B = np.power(sqrt(reverse_barrier)+sqrt(barrier),2)

	y = - np.exp(2.0*np.pi*(x-xts)/length)

	V = - A*y*np.power(1.0-y,-1) - B*y*np.power(1.0-y,-2)

	return V
