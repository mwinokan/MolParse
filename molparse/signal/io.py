
def csv_strip(filename,output=None,overwrite=False,comment_chars='#@!'):
	import re
	import mout
	import mcol

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

			if line[0] in comment_chars:
				continue

			line = re.sub(r'\t', ' ', line)
			line = re.sub(' +',' ',line)
			line = line.strip()
			file.write(line+'\n')
	else:
		mout.warningOut("Skipping already stipped file: "+filename)

	file.close()

def parseDat(filename,num_columns=2,header_rows=1,delimiter=' ',debug=False,pre_strip=False,clean_nan=False,comment_chars='#@!'):
	
	import mout
	import pandas
	import os
	import math

	columns=list(range(num_columns))

	labels=[]

	for i in columns:

		if i==0:
			labels.append("x")
		else:
			labels.append("y"+str(i))

	if pre_strip:
		# csv_strip(filename,overwrite=True)
		csv_strip(filename,output="__temp__",overwrite=False,comment_chars=comment_chars)
		
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
									names=labels,comment='#')

	if debug: print(dataframe)

	if num_columns > 2:

		if debug: mout.headerOut("many")

		big_y = []
		x = dataframe.iloc[:, 0].values

		for i in range(1,num_columns):

			y = list(dataframe.iloc[:, i].values)

			big_y.append(y)

			if debug: mout.varOut("big_y len",len(big_y))

		if clean_nan:
			if num_columns > 2:
				new_x = []
				new_ys = []
				for x,ys in zip(x,big_y):

					# print("DEBUG",x,ys)

					try:
						if math.isnan(float(x)):
							continue
					except ValueError:
						continue
					for y in ys:
						try:
							if math.isnan(float(y)):
								continue
						except ValueError:
							continue
							
					new_x.append(x)
					new_ys.append(ys)
				x = new_x
				big_y = new_ys

		return x,big_y

	elif num_columns == 1:

		x = dataframe.iloc[:, 0].values

		if clean_nan:
			new_x = []
			for this_x in list(x):
				if math.isnan(float(this_x)):
					continue
				new_x.append(this_x)
			return new_x

		return list(x)

	else:

		x = dataframe.iloc[:, 0].values
		y = dataframe.iloc[:, 1].values

		if clean_nan:
			new_x = []
			new_y = []
			for this_x,this_y in zip(list(x),list(y)):
				if math.isnan(float(this_x)):
					continue
				if math.isnan(float(this_y)):
					continue
				new_x.append(this_x)
				new_y.append(this_y)
			return new_x, new_y

		return list(x),list(y)
