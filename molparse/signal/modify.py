
def runningAverage(xdata,ydata,averaging_window=5,cutoff_style=None,debug=False):
	import mout

	many = any(isinstance(el,list) for el in ydata)

	if debug: 
		mout.varOut("many",many)
		mout.varOut("len(xdata)",len(xdata))
		mout.varOut("len(ydata)",len(ydata))

	if many:

		result_x = []
		result_y = []

		for dataset in ydata:

			x,y = runningAverage(xdata,dataset,averaging_window=averaging_window,cutoff_style=cutoff_style)
			result_y.append(y)

		result_x = x

		assert len(ydata) == len(result_y)

	else:

		if cutoff_style is not None:
			mout.errorOut("Not supported yet!",fatal=True)

		result_x = []
		result_y = []

		for i,x in enumerate(xdata): 

			result_x.append(x)

			if i == 0:

				y = ydata[0]

			elif i < averaging_window: 

				y = sum(ydata[0:i])/i

			else:

				y = sum(ydata[i-averaging_window:i])/averaging_window

			result_y.append(y)

	return result_x,result_y

def differentiate(xdata,ydata):

	many = any(isinstance(el,list) for el in ydata)

	result_x = []
	result_y = []

	if many:

		for dataset in ydata:

			x,y = differentiate(xdata,dataset)

			result_y.append(y)

		result_x = x

		assert len(ydata) == len(result_y)

	else:

		for i,x in enumerate(xdata):
			result_x.append(xdata[i])

			if i!=0: 

				dy_dx = (ydata[i] - ydata[i-1]) / (xdata[i] - xdata[i-1])

				result_y.append(dy_dx)

			else:

				result_y.append(0.0)

	return result_x,result_y
