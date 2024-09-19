def peakFinder(
    xdata,
    ydata,
    min_width,
    min_height,
    baseline=None,
    threshold=1.0e-2,
    search_coeff=0.1,
):
    import statistics
    import mout
    from .modify import differentiate

    many = any(isinstance(el, list) for el in ydata)

    xdata, dydx = differentiate(xdata, ydata)

    if many:

        peaklist = []

        for data in ydata:
            peaklist.append(
                peakFinder(
                    xdata,
                    data,
                    min_width,
                    min_height,
                    baseline=baseline,
                    threshold=threshold,
                    search_coeff=search_coeff,
                )
            )

        all_peaks = peaklist

    else:

        if baseline is None:
            y_avg = statistics.mean(ydata)
        else:
            y_avg = baseline

        # print(y_avg)

        found = False

        all_peaks = []

        for i, y in enumerate(ydata):

            this_height = y - y_avg

            if this_height > min_height and abs(dydx[i]) < threshold:

                # print(xdata[i],y,this_height,abs(dydx[i]))

                skip = False
                for peak in all_peaks:
                    if peak.in_range(xdata[i]):
                        skip = True
                        break

                if skip:
                    continue

                # backwards search for the start of the peak
                backsearch_i = i
                backsearch_height = this_height

                while backsearch_height > min_height * search_coeff:

                    backsearch_i -= 1
                    if backsearch_i <= 0:
                        break
                    # print(backsearch_i)
                    backsearch_height = ydata[backsearch_i] - y_avg

                mout.varOut("Peak start estimated at", xdata[backsearch_i])

                # forwards search for the start of the peak
                foresearch_i = i
                foresearch_height = this_height

                while foresearch_height > min_height * search_coeff:

                    foresearch_i += 1
                    if foresearch_i >= len(ydata) - 1:
                        break
                    foresearch_height = ydata[foresearch_i] - y_avg

                mout.varOut("Peak end estimated at", xdata[foresearch_i])

                all_peaks.append(Peak(xdata[backsearch_i], xdata[foresearch_i]))

    return all_peaks


class Peak:
    def __init__(self, start_x, end_x):

        self.start_x = start_x
        self.end_x = end_x

    def __str__(self):
        return "Peak[" + str(self.start_x) + ":" + str(self.end_x) + "]"

    def __repr__(self):
        return "Peak[" + str(self.start_x) + ":" + str(self.end_x) + "]"

    def in_range(self, x):
        if x >= self.start_x and x <= self.end_x:
            return True
        else:
            return False

    def gauss_fit(self, xdata, ydata, return_data=False, baseline=0.0, filename=None):
        return gaussian_fit(
            xdata,
            ydata,
            [self.start_x, self.end_x],
            return_data=return_data,
            baseline=baseline,
            filename=filename,
        )


def gaussian_fit(
    xdata,
    ydata,
    window,
    return_data=False,
    index_window=False,
    baseline=0.0,
    filename=None,
):
    import mout
    from scipy.optimize import curve_fit
    from scipy import asarray as ar

    many = any(isinstance(el, list) for el in ydata)

    if many:
        mout.errorOut("gaussian_fit. Unsupported.", fatal=True)

    if index_window:
        x = ar(xdata[window[0] : window[1]])
        y = ar(ydata[window[0] : window[1]])
    else:
        index1 = closest_index(window[0], xdata)
        index2 = closest_index(window[1], xdata)
        x = ar(xdata[index1:index2])
        y = ar(ydata[index1:index2])

    n = len(x)  # the number of data
    mean = sum(x * y) / n
    sigma = sum(y * (x - mean) ** 2) / n

    mean = x[len(x) // 2]
    sigma = (x[-1] - x[0]) / 4

    import mplot

    # print(1.0,mean,sigma)

    global ___BASELINE
    ___BASELINE = baseline

    try:
        popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])
    except:
        mplot.graph2D(x, y, filename=filename, show=False)
        return False

    a, mean, sigma = popt

    # print(a,mean,sigma)

    if return_data:
        gaus_y = []
        for this_x in xdata:
            gaus_y.append(gaus(this_x, *popt))
        if filename is not None:
            mplot.graph2D(x, [y, gaus_y[index1:index2]], filename=filename, show=False)
        return a, mean, sigma, xdata, gaus_y
    else:
        return a, mean, sigma


def gaus(x, a, x0, sigma):
    from scipy import asarray as exp

    return a * exp(-((x - x0) ** 2) / (2 * sigma**2)) + ___BASELINE


# def closest_index(value,xdata):
# 	"""Return the index of the nearest data point in the array"""
# 	if isinstance(value,list):
# 		result = []
# 		for v in value:
# 			result.append(closest_index(v, xdata))
# 		return result
# 	else:
# 		for i,x in enumerate(xdata):
# 			if x > value:
# 				if abs(xdata[i-1]-value) < abs(xdata[i]-value):
# 					return i
# 				else:
# 					return i-1


def closest_index(value, xdata, comparison=None, numpy=False):
    """Return the index of the nearest data point in the array"""
    if numpy:
        if comparison:
            import mout

            mout.warningOut(
                "numpy argument overrides comparison",
                code="amp.signal.peaks.closest_index[1]",
            )
        import numpy as np

        array = np.asarray(xdata)
        return (np.abs(array - value)).argmin()
    elif isinstance(value, list):
        result = []
        for v in value:
            result.append(closest_index(v, xdata, comparison=comparison))
        return result
    else:
        for i, x in enumerate(xdata):
            if (i > 0 and x < xdata[i - 1]) or (i == 0 and x > xdata[i + 1]):
                import mout

                mout.errorOut("Array does not appear to be sorted!")
                return -1
            elif x > value:
                if comparison is not None:
                    if comparison == "<":
                        return i - i
                    else:
                        return i
                else:
                    if abs(xdata[i - 1] - value) < abs(xdata[i] - value):
                        return i
                    else:
                        return i - 1
        return i


def closest_value(
    xvalue, xdata, ydata=None, comparison=None, numpy=False, return_x=False
):
    """Return the yvalue of the data point nearest the given xvalue"""
    if ydata is None:
        ydata = xdata
    if isinstance(xvalue, list):
        x_result = []
        y_result = []
        for x in xvalue:
            xres, yres = closest_value(
                x, xdata, ydata, comparison=comparison, numpy=numpy, return_x=return_x
            )
            x_result = []
            y_result = []
        if return_x:
            return x_result, y_result
        else:
            return y_result
    else:
        index = closest_index(xvalue, xdata, comparison=comparison, numpy=numpy)
        if return_x:
            return xdata[index], ydata[index]
        else:
            return ydata[index]
