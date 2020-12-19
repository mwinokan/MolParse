
import pandas
import mout

def parseDat(filename,num_columns=2,header_rows=1,delimiter=' ',debug=False):

	columns=list(range(num_columns))

	labels=[]

	for i in columns:

		if i==0:
			labels.append("x")
		else:
			labels.append("y"+str(i))

	dataframe = pandas.read_csv(filename,
								skiprows=header_rows,
								delimiter=delimiter,
								usecols=columns,
								names=labels)

	if num_columns > 2:

		if debug: mout.headerOut("many")

		big_y = []
		x = dataframe.iloc[:, 0].values

		for i in range(1,num_columns):

			y = list(dataframe.iloc[:, i].values)

			big_y.append(y)

			if debug: mout.varOut("big_y len",len(big_y))

		return x,big_y

	else:


	# x = dataframe.columns[0].values
	# y = dataframe.columns[1].values

		x = dataframe.iloc[:, 0].values
		y = dataframe.iloc[:, 1].values

		return x,y
