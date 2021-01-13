
import pandas
import mout
import re
import mcol

def csv_strip(filename,output=None,overwrite=False):

	with open(filename,'r') as file:
		lines = file.readlines()
	file.close()

	if output is None and not overwrite:
		mout.errorOut("No output specified",fatal=True)

	if overwrite:
		mout.warningOut("Overwriting "+mcol.file+filename+mcol.warning+"!")
		output=filename

	file = open(output,'w')

	for line in lines:
		line = re.sub(' +',' ',line)
		line = line.strip()
		file.write(line+'\n')

	file.close()

def parseDat(filename,num_columns=2,header_rows=1,delimiter=' ',debug=False,pre_strip=False):

	columns=list(range(num_columns))

	labels=[]

	for i in columns:

		if i==0:
			labels.append("x")
		else:
			labels.append("y"+str(i))

	if pre_strip:
		csv_strip(filename,overwrite=True)

	dataframe = pandas.read_csv(filename,
								skiprows=header_rows,
								delimiter=delimiter,
								usecols=columns,
								names=labels)

	if debug: print(dataframe)

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

		return list(x),list(y)
