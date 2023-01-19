
# import math


def find_barrier_stationary_points(xdata,ydata,yerr=None,reduction_order=2,show=False,tshicut=None,tslocut=None):
	import numpy as np
	import scipy.signal as sps

	# find the x-location of the stationary points

	if yerr is None:
		yerr = [None for x in xdata]

	maxima_indices = sps.argrelextrema(np.array(ydata), np.greater)[0]
	minima_indices = sps.argrelextrema(np.array(ydata), np.less)[0]

	max_coords = []
	for i in maxima_indices:
		max_coords.append([xdata[i],ydata[i],yerr[i]])

	min_coords = []
	for i in minima_indices:
		min_coords.append([xdata[i],ydata[i],yerr[i]])

	if reduction_order == 2:

		maxima_ys = []
		for p in max_coords:
			valid = True
			if tshicut is not None and p[0] >= tshicut:
				valid = False
			elif tslocut is not None and p[0] <= tslocut:
				valid = False

			if valid:
				maxima_ys.append(p[1])
			else:
				maxima_ys.append(-1)

		# if tshicut is not None:
		# 	if tslocut is not None:
		# 		maxima_ys = [p[1] for p in max_coords if p[0] <= tshicut and p[0] >= tslocut]
		# 	else:
		# 		maxima_ys = [p[1] for p in max_coords if p[0] <= tshicut]
		# 		# maxima_ys = [p[1] for p in max_coords if p[0] >= tslocut]
		# else:
		# 	maxima_ys = [p[1] for p in max_coords]

		if not maxima_ys:
			import mout
			mout.warningOut("No maximas found!")
			# import matplotlib.pyplot as plt
			# print(min_coords)
			# plt.plot(xdata,ydata)
			# plt.scatter([c[0] for c in min_coords],[c[1] for c in min_coords])
			# plt.show()
			assert len(min_coords) == 1
			return min_coords, []

		highest_maxima_index = maxima_ys.index(max(maxima_ys))

		max_coords = [max_coords[highest_maxima_index]]

		maxima_x = max_coords[0][0]

		min_before = [p for p in min_coords if p[0] < maxima_x]
		min_after = [p for p in min_coords if p[0] > maxima_x]

		min_before = [p for p in min_before if p[1] == min([p[1] for p in min_before])]
		min_after = [p for p in min_after if p[1] == min([p[1] for p in min_after])]

		min_coords = min_before + min_after

	if show:
		import matplotlib.pyplot as plt
		fig,ax = plt.subplots()

		plt.plot(xdata,ydata)
		
		plt.errorbar([p[0] for p in max_coords],[p[1] for p in max_coords],yerr=[p[2] for p in max_coords],fmt="o")
		plt.errorbar([p[0] for p in min_coords],[p[1] for p in min_coords],yerr=[p[2] for p in min_coords],fmt="o")

		plt.show()

	return min_coords,max_coords

# def eckartFit(xdata,ydata):
# 	return None

def sqrt(x):
	import numpy as np
	return np.power(x,0.5)

def eckart(x,barrier,asymmetry,length,xts=0.0):
	import numpy as np

	if isinstance(x,list):
		x = np.array(x)

	reverse_barrier = barrier - asymmetry

	# print(reverse_barrier)

	A = barrier - reverse_barrier

	B = np.power(sqrt(reverse_barrier)+sqrt(barrier),2)

	y = - np.exp(2.0*np.pi*(x-xts)/length)

	V = - A*y*np.power(1.0-y,-1) - B*y*np.power(1.0-y,-2)

	return V
