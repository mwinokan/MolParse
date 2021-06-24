
import pandas
import mout
import re
import mcol
import os

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

	test_line = re.sub(r'\t', ' ', lines[0])
	test_line = re.sub(' +',' ',test_line)
	test_line = test_line.strip()

	if test_line != lines[0]:
		for line in lines:
			line = re.sub(r'\t', ' ', line)
			line = re.sub(' +',' ',line)
			line = line.strip()
			file.write(line+'\n')
	else:
		mout.warningOut("Skipping already stipped file: "+filename)

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
		# csv_strip(filename,overwrite=True)
		csv_strip(filename,output="__temp__",overwrite=False)
		
		dataframe = pandas.read_csv("__temp__",
									skiprows=header_rows,
									delimiter=delimiter,
									usecols=columns,
									names=labels)

		os.system("rm __temp__")

	else:

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
