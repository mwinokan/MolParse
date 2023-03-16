
import mout
import mcol

def parseXVG(file):
	"""Parse a GROMACS XVG file and return a molparse.XVG object"""

	mout.out(f"Parsing {mcol.file}{file}{mcol.clear}...")
	
	header_buffer = []
	xvg = None
	many = []
	multiname = None

	with open(file) as f:
		for line in f:

			# regular comments
			if line.startswith('#'):
				continue

			elif line.startswith('@'):
				header_buffer.append(line[1:].strip())
				continue

			elif line.startswith('&'):
				if many:
					if xvg.entries:
						many.append(xvg)
						many[-1].title = f"[{len(many)}]"
						xvg = xvg.create_blank_copy()
				else:
					many = [xvg]
					multiname = xvg.title
					many[-1].title = f"[{len(many)}]"
					xvg = xvg.create_blank_copy()
				continue

			elif not xvg:
				xvg = XVG()
				xvg.determine_data_shape(header_buffer,line)
				xvg.parse_data_line(line)

			else:
				xvg.parse_data_line(line)

	if many:
		if xvg.entries:
			many.append(xvg)
			many[-1].title = f"[{len(many)}]"
		return XVGCollection(many,multiname)
	else:
		return xvg
		
class XVG():
	"""Class to store data obtained from an XVG"""
	def __init__(self):
		self.entries = 0
	
	def determine_data_shape(self,header_buffer,demo_line):
		"""Use the header strings and an example data line to construct the data shape"""
		self.title = [line.lstrip("title ") for line in header_buffer if line.startswith("title")]
		self.xlabel = [line.lstrip("xaxis ") for line in header_buffer if line.startswith("xaxis")]
		self.ylabel = [line.lstrip("yaxis ") for line in header_buffer if line.startswith("yaxis")]
		self.type = [line.lstrip("TYPE ") for line in header_buffer if line.startswith("TYPE")]

		if len(self.title) > 1: mout.warningOut("Multiple title headers")
		if len(self.xlabel) > 1: mout.warningOut("Multiple xaxis headers")
		if len(self.ylabel) > 1: mout.warningOut("Multiple yaxis headers")
		if len(self.type) > 1: mout.warningOut("Multiple TYPE headers")

		self.title = self.title[0]
		self.xlabel = self.xlabel[0].lstrip('label "').rstrip('"')
		self.ylabel = self.ylabel[0].lstrip('label "').rstrip('"')
		self.type = self.type[-1]
		
		split_line = demo_line.strip().split()

		self.num_columns = len(split_line)

		if self.type == 'xy' and self.num_columns == 2:
			self.column_labels = ['x','y']
			self.column_types = [float,float]
			self.columns = {'x':[],'y':[]}
		elif self.type == 'xy':
			self.column_labels = ['x'] + [f'y{i+1}' for i in range(self.num_columns-1)]
			self.column_types = [float]*self.num_columns
			self.columns = {key:[] for key in self.column_labels}
		elif self.type == 'xydy':
			self.column_labels = ['x','y','yerr']
			self.column_types = [float,float,float]
			self.columns = {'x':[],'y':[],'yerr':[]}
		else:
			mout.errorOut(f"Unrecognised XVG type: {self.type}, with {self.num_columns} columns",fatal=True)

	def parse_data_line(self,line):
		"""Parse a line of XVG data"""

		split_line = line.strip().split()

		for value,column,data_type in zip(split_line,self.column_labels,self.column_types):
			self.columns[column].append(data_type(value))

		self.entries += 1

	def plot(self,show=False):
		"""Use matplotlib to plot the XVG data"""
		
		import matplotlib.pyplot as plt
		fig,ax = plt.subplots()
		ax.set_xlabel(self.xlabel)
		ax.set_ylabel(self.ylabel)

		if self.type == 'xy' and self.num_columns == 2:
			ax.plot(self.columns['x'],self.columns['y'])

		elif self.type == 'xy':
			for key in self.column_labels[1:]:
				ax.plot(self.columns['x'],self.columns[key],label=key)

		elif self.type == 'xydy' and self.num_columns == 3:
			ax.errorbar(self.columns['x'],self.columns['y'],yerr=self.columns['yerr'])

		else:
			mout.errorOut(f"Unrecognised XVG type (XVG.plot): {self.type}, with {self.num_columns} columns",fatal=True)

		if show:
			plt.show()
		return fig,ax

	def get_go_trace(self,name=None):

		import plotly.graph_objects as go

		if self.type == 'xy' and self.num_columns == 2:
			trace = go.Scatter(x=self.columns['x'],
							   y=self.columns['y'],
							   name=name,
							   mode='lines')
		elif self.type == 'xy':
			traces = []
			for key in self.column_labels[1:]:
				trace = go.Scatter(x=self.columns['x'],
								   y=self.columns[key],
								   mode='lines',
								   name=key)
				traces.append(trace)
			return traces
		elif self.type == 'xydy' and self.num_columns == 3:
			trace = go.Scatter(x=self.columns['x'],
							   y=self.columns['y'],
							   name=name,
							   error_y=dict(type='data',
							   				array=self.columns['yerr'],
							   				visible=True),
							   mode='lines')
		else:
			mout.errorOut(f"Unrecognised XVG type (XVG.plotly): {self.type}, with {self.num_columns} columns",fatal=True)

		return trace

	def plotly(self,show=False):
		"""Use plotly to plot the XVG data"""

		import plotly.graph_objects as go

		fig = go.Figure()

		fig.update_layout(
			title=self.title,
			xaxis_title=self.xlabel,
			yaxis_title=self.ylabel,
			font=dict(family="Helvetica Neue",size=18)
		)

		trace = self.get_go_trace()
		fig.add_trace(trace)

		if show:
			fig.show()
		return fig

	def __repr__(self):
		return f'XVG({self.title=}, {self.type=}, {self.xlabel=}, {self.ylabel=}, {self.entries=}, {self.type=}, {self.column_labels=}, {self.column_types=})'

	def create_blank_copy(self):
		"""Create a copy of this data structure but with no entries"""

		copy = XVG()

		copy.title = self.title
		copy.xlabel = self.xlabel
		copy.ylabel = self.ylabel
		copy.type = self.type
		copy.column_labels = self.column_labels
		copy.column_types = self.column_types

		copy.columns = {}
		for key in self.columns:
			copy.columns[key] = []
		
		copy.num_columns = self.num_columns

		return copy

class XVGCollection():
	"""Class to store many XVG objects"""
	def __init__(self,xvg_list,title):
		self.children = xvg_list
		self.title = title

	def plotly(self,show=False,statistics=False):

		import plotly.graph_objects as go

		fig = go.Figure()

		fig.update_layout(
			title=self.title,
			xaxis_title=self.children[0].xlabel,
			yaxis_title=self.children[0].ylabel,
			font=dict(family="Helvetica Neue",size=18)
		)

		if not statistics:

			for xvg in self.children:

				fig.add_trace(xvg.get_go_trace(name=xvg.title))

		else:

			import numpy as np

			x = self.children[0].columns['x']
			mean = []
			mean_plus_std = []
			mean_minus_std = []
			mins = []
			maxs = []

			# calculate mean curve
			for i in range(self.children[0].entries):
				this_slice = [xvg.columns['y'][i] for xvg in self.children]

				mins.append(min(this_slice))
				maxs.append(max(this_slice))

				mu = np.mean(this_slice)
				std = np.std(this_slice)

				mean.append(mu)
				mean_plus_std.append(mu+std)
				mean_minus_std.append(mu-std)

			fig.add_trace(go.Scatter(x=x,y=maxs,name="max",line=dict(width=0,color='black')))
			fig.add_trace(go.Scatter(x=x,y=mins,name="min",fill='tonexty',line=dict(width=0,color='black'),fillcolor='rgba(0,0,0,0.15)'))
			fig.add_trace(go.Scatter(x=x,y=mean_plus_std,name="mean+std",line=dict(width=0,color='black')))
			fig.add_trace(go.Scatter(x=x,y=mean_minus_std,name="mean-std",fill='tonexty',line=dict(width=0),fillcolor='rgba(0,0,0,0.3)'))
			fig.add_trace(go.Scatter(x=x,y=mean,name="mean",line=dict(color="black",width=4)))

			fig.update_layout(showlegend=False)

		if show:
			fig.show()

		return fig

	def align_ydata(self,function):

		import numpy as np

		for xvg in self.children:

			data = xvg.columns['y']
			ref_value = function(data)
			array = np.array(data)
			array -= ref_value
			xvg.columns['y'] = list(array)


