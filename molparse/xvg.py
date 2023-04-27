
import mout
import mcol

def parseXVG(file,xmin=None,xmax=None,convert_nanometres=True,yscale=1.0,ylabel=None,no_com=False,no_ref=False,ymin=None,ymax=None):
	"""Parse a GROMACS XVG file and return a molparse.XVG object"""

	mout.out(f"Parsing {mcol.file}{file}{mcol.clear}...")
	
	header_buffer = []
	xvg = None
	many = []
	multiname = None
	contains_nan = False

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
				xvg = XVG(ymin=ymin,ymax=ymax)
				if ylabel:
					xvg.ylabel = ylabel
				xvg.determine_data_shape(header_buffer,line,convert_nanometres,xmin=xmin,xmax=xmax,yscale=yscale)

				xvg.title = f'{file}: {xvg.title}'

				status = xvg.parse_data_line(line)
				contains_nan = contains_nan or not status

			else:
				status = xvg.parse_data_line(line)
				contains_nan = contains_nan or not status

	if many:
		for xvg in many:
			if no_com:
				xvg.remove_com_columns()
			if no_ref:
				xvg.remove_ref_columns()
		many.append(xvg)
		many[-1].title = f"[{len(many)}]"
		if contains_nan:
			mout.warningOut("At least one NaN in XVG!")
		return XVGCollection(many,multiname)
	else:
		if no_com:
			xvg.remove_com_columns()
		if no_ref:
			xvg.remove_ref_columns()
		if contains_nan:
			mout.warningOut("At least one NaN in XVG!")
		return xvg
		
def clean_label(string):
	"""Parses Grace label syntax into LaTeX"""

	substring = string
	chunks = []

	while substring:

		for key in _all_codes:
			if substring.startswith(key):
				chunks.append(key)
				substring = substring[len(key):]
				break

		else:
			chunks.append(substring[0])
			substring = substring[1:]

	chunks = [c for c in chunks if c not in _ignorable_codes]

	string = ''
	mode = None

	# sub/superscripts
	for chunk in chunks:

		if chunk == "\\s":
			mode = 'sub'
			string += r'}_{\text{'
		elif chunk == "\\S":
			mode = 'super'
			string += r'}^{\text{'
		elif chunk == "\\N":
			mode = None
			string += r'}}\text{'
		else:
			string += chunk

	string = r'$\text{'+string+'}$'

	return string

class XVG():
	"""Class to store data obtained from an XVG"""
	def __init__(self,ymin=None,ymax=None):
		self.entries = 0
		self._minima_indices = None
		self._maxima_indices = None
		self.ymin = ymin
		self.ymax = ymax

	@property
	def dataframe(self):
		import pandas
		return pandas.DataFrame(data=self.columns)
	
	def determine_data_shape(self,header_buffer,demo_line,convert_nanometres,xmin,xmax,yscale=1.0):
		"""Use the header strings and an example data line to construct the data shape"""
		self.title = [line.lstrip("title ") for line in header_buffer if line.startswith("title")]
		self.xlabel = [line.lstrip("xaxis ") for line in header_buffer if line.startswith("xaxis")]
		self.ylabel = [line.lstrip("yaxis ") for line in header_buffer if line.startswith("yaxis")]
		self.type = [line.lstrip("TYPE ") for line in header_buffer if line.startswith("TYPE")]
		self.series_labels = [line for line in header_buffer if line.startswith("s") and 'legend' in line]

		if self.series_labels:
			self.series_labels = {int(l.split("legend")[0].strip().lstrip("s")):l.split("legend")[1].strip().replace('"','') for l in self.series_labels}

		if len(self.title) > 1: mout.warningOut("Multiple title headers")
		if len(self.xlabel) > 1: mout.warningOut("Multiple xaxis headers")
		if len(self.ylabel) > 1: mout.warningOut("Multiple yaxis headers")
		if len(self.type) > 1: mout.warningOut("Multiple TYPE headers")

		self.xmin = xmin
		self.xmax = xmax

		self.title = self.title[0]
		self.xlabel = self.xlabel[0].lstrip('label "').rstrip('"')
		self.ylabel = self.ylabel[0].lstrip('label "').rstrip('"')

		self.xlabel = clean_label(self.xlabel)
		self.ylabel = clean_label(self.ylabel)

		self.type = self.type[-1]
		
		self.xscale = 1.0
		self.yscale = yscale
		if convert_nanometres and '(nm)' in self.xlabel:
			self.xscale = 10.0
			self.xlabel = self.xlabel.replace('(nm)','(Å)')

		split_line = demo_line.strip().split()

		self.num_columns = len(split_line)

		if self.type == 'xy' and self.num_columns == 2:
			self.column_labels = ['x','y']
			self.column_types = [float,float]
			self.column_scales = [self.xscale,self.yscale]
			self.columns = {'x':[],'y':[]}
		elif self.type == 'xy':
			if self.series_labels:
				self.column_labels = ['x'] + [v for _,v in self.series_labels.items()]
			else:
				self.column_labels = ['x'] + [f'y{i+1}' for i in range(self.num_columns-1)]
			self.column_scales = [self.xscale] + [self.yscale]*(self.num_columns-1)
			self.column_types = [float]*self.num_columns
			self.columns = {key:[] for key in self.column_labels}
		elif self.type == 'xydy':
			self.column_labels = ['x','y','yerr']
			self.column_scales = [self.xscale,self.yscale,self.yscale]
			self.column_types = [float,float,float]
			self.columns = {'x':[],'y':[],'yerr':[]}
		else:
			mout.errorOut(f"Unrecognised XVG type: {self.type}, with {self.num_columns} columns",fatal=True)

	def parse_data_line(self,line):
		"""Parse a line of XVG data"""

		split_line = line.strip().split()

		if len(split_line) < len(self.column_labels):
			mout.warningOut("Skipping incomplete data-line")
			return False

		for value,column,data_type,scale in zip(split_line,self.column_labels,self.column_types,self.column_scales):
			if scale is not None:
				value = data_type(value)*scale
			else:
				value = data_type(value)

			if self.xmin and column == 'x' and float(self.xmin) > value:
				return True
			if self.xmax and column == 'x' and float(self.xmax) < value:
				return True

			self.columns[column].append(value)

		self.entries += 1

		if 'nan' in line:
			return False
		else:
			return True

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

	@property
	def y_columns(self):
		return {k:y for k,y in self.columns.items() if k != 'x'}

	@property
	def summable(self):
		if self.type == 'xy' and self.num_columns > 2:
			return True
		return False
	
	def get_summed_trace(self,name=None,color='black',normalised=False):
		import plotly.graph_objects as go

		if self.type == 'xy' and self.num_columns > 2:

			summed_data = []
			for i in range(self.entries):
				
				value = 0
				for key in self.y_columns:
					value += self.y_columns[key][i]

				summed_data.append(value)

			if normalised:
				maxval = max(summed_data)
				summed_data = [x/maxval for x in summed_data]

			trace = go.Scatter(x=self.columns['x'],
				   y=summed_data,
				   line=dict(color=color),
				   name=self.title+' (Summed)',
				   mode='lines',
				   )

			return trace
		else:
			mout.errorOut(f"Unrecognised XVG type (XVG.get_summed_trace): {self.type}, with {self.num_columns} columns",fatal=True)

	def get_go_trace(self,name=None,color=None,group_from_title=False,histogram=False,normalised=False):
		"""get the plotly.go trace(s) for the data"""

		import plotly.graph_objects as go

		if normalised:
			normalised = 'probability'

		if self.type == 'xy' and self.num_columns == 2:
			if histogram:
				ymin = float(self.ymin or min(self.columns['y']))
				ymax = float(self.ymax or max(self.columns['y']))
				trace = go.Histogram(x=self.columns['y'],
								   	 marker_color=color,
								   	 histnorm=normalised,
								     xbins=dict(start=ymin,end=ymax,size=(ymax-ymin)/float(histogram)),
								   	 name=name)
			else:
				trace = go.Scatter(x=self.columns['x'],
								   y=self.columns['y'],
								   line=dict(color=color),
								   name=name,
								   mode='lines',
								   )
		elif self.type == 'xy':
			traces = []
			for key in self.column_labels[1:]:

				if group_from_title:
					group = self.title

				else:

					group = "default"

					if 'pull' in self.title.lower():
						group = 'RC'

					if 'ref' in key:
						group = 'REF'
					elif 'X' in key or 'Y' in key or 'Z' in key:
						group = 'COM'

				if histogram:
					ymin = float(self.ymin or min(self.columns[key]))
					ymax = float(self.ymax or max(self.columns[key]))
					trace = go.Histogram(x=self.columns[key],
									   marker_color=color,
									   marker_opacity=0.6,
									   name=key,
									   histnorm=normalised,
									   xbins=dict(start=ymin,end=ymax,size=(ymax-ymin)/float(histogram)),
									   legendgroup=group,
									   legendgrouptitle_text=group)
				else:
					trace = go.Scatter(x=self.columns['x'],
									   y=self.columns[key],
									   mode='lines',
									   line=dict(color=color),
									   name=key,
									   legendgroup=group,
									   legendgrouptitle_text=group)

				traces.append(trace)
			return traces
		elif self.type == 'xydy' and self.num_columns == 3:
			if histogram:
				ymin = float(self.ymin or min(self.columns['y']))
				ymax = float(self.ymax or max(self.columns['y']))
				trace = go.Histogram(x=self.columns['y'],
								   name=name,
								   histnorm=normalised,
								   marker_color=color,
								   )
			else:
				trace = go.Scatter(x=self.columns['x'],
								   y=self.columns['y'],
								   name=name,
								   line=dict(color=color),
								   error_y=dict(type='data',
								   				array=self.columns['yerr'],
								   				visible=True),
								   mode='lines')
		else:
			mout.errorOut(f"Unrecognised XVG type (XVG.plotly): {self.type}, with {self.num_columns} columns",fatal=True)

		return trace

	def remove_com_columns(self):
		self.columns = {key:column for key,column in self.columns.items() if not any(['X' in key, 'Y' in key, 'Z' in key])}		
		self.column_labels = [key for key in self.column_labels if not any(['X' in key, 'Y' in key, 'Z' in key])]

	def remove_ref_columns(self):
		self.columns = {key:column for key,column in self.columns.items() if not 'ref' in key}		
		self.column_labels = [key for key in self.column_labels if not 'ref' in key]

	def calculate_stationary_points(self,column='y'):
		"""find any stationary points"""

		import numpy as np
		import scipy.signal as sps

		self.minima_indices = sps.argrelextrema(np.array(self.columns[column]), np.less)[0]
		self.maxima_indices = sps.argrelextrema(np.array(self.columns[column]), np.greater)[0]
		self.stationary_column = column

	@property
	def minima_indices(self):
		"""List of indices pertaining to minima"""
		if self._minima_indices is None:
			self.calculate_stationary_points()
		return self._minima_indices

	@property
	def maxima_indices(self):
		"""List of indices pertaining to maxima"""
		if self._maxima_indices is None:
			self.calculate_stationary_points()
		return self._maxima_indices

	@minima_indices.setter
	def minima_indices(self,a):
		self._minima_indices = a

	@maxima_indices.setter
	def maxima_indices(self,a):
		self._maxima_indices = a

	@property
	def stationary_points(self):
		"""build the stationary points dictionary"""
		dictionary = {'minima':[],'maxima':[]}
		for i in self.minima_indices:
			dictionary['minima'].append(dict(type='min',index=i,x=self.columns['x'][i],y=self.columns[self.stationary_column][i]))
		for i in self.maxima_indices:
			dictionary['maxima'].append(dict(type='max',index=i,x=self.columns['x'][i],y=self.columns[self.stationary_column][i]))
		return dictionary

	def get_stationary_point_trace(self,name=None,column='y',color=None):
		"""find any stationary points and make a plotly.go.Scatter trace"""
		
		import plotly.graph_objects as go

		if self.minima_indices is None or self.maxima_indices is None:
			self.calculate_stationary_points()
		
		indices = list(self.minima_indices) + list(self.maxima_indices)
		xdata = [self.columns['x'][i] for i in indices]
		ydata = [self.columns[column][i] for i in indices]

		return go.Scatter(x=xdata,y=ydata,name=name,mode='markers',marker=dict(color=color))

	def get_closest_value(self,x,column='y',return_x=False):
		from .signal import closest_value
		return closest_value(x,self.columns['x'],self.columns[column],return_x=return_x)

	def plotly(self,show=False,color=None,fig=None,group_from_title=False,histogram=False,summed=False,normalised=False):
		"""Use plotly to plot the XVG data"""

		import plotly.graph_objects as go

		if histogram and normalised:
			mout.errorOut("Not supported!",fatal=True)

		if fig is None:
			fig = go.Figure()

		fig.update_layout(
			title=self.title,
			legend=dict(groupclick="toggleitem"),
			font=dict(family="Helvetica Neue",size=18)
		)

		if histogram:
			fig.update_layout(xaxis_title=self.ylabel,yaxis_title='Count',barmode='overlay')
			fig.update_traces(opacity=0.75)
		else:
			fig.update_layout(xaxis_title=self.xlabel, yaxis_title=self.ylabel)

		if summed:
			trace = self.get_summed_trace(name=self.title,color=color,normalised=normalised)
		else:
			trace = self.get_go_trace(name=self.title,color=color,group_from_title=group_from_title,histogram=histogram,normalised=normalised)

		if isinstance(trace, list):
			for t in trace:
				fig.add_trace(t)
		else:
			fig.add_trace(trace)

		if show:
			fig.show()
		return fig

	def __repr__(self):
		return f'XVG({self.title=}, {self.type=}, {self.xlabel=}, {self.ylabel=}, {self.entries=}, {self.type=}, {self.column_labels=}, {self.column_types=})'

	def create_blank_copy(self):
		"""Create a copy of this data structure but with no entries"""

		copy = XVG(ymin=self.ymin,ymax=self.ymax)

		copy.title = self.title
		copy.xlabel = self.xlabel
		copy.ylabel = self.ylabel
		copy.type = self.type
		copy.column_labels = self.column_labels
		copy.column_types = self.column_types
		copy.column_scales = self.column_scales
		
		copy.xmin = self.xmin
		copy.xmax = self.xmax

		copy.columns = {}
		for key in self.columns:
			copy.columns[key] = []
		
		copy.num_columns = self.num_columns

		return copy

	def align_ydata(self,method,column='y'):
		"""Shift the y-values of all the column to align them according to a function such as min/max or an x-coordinate"""

		import numpy as np
		from .signal import closest_value

		data = self.columns[column]

		if isinstance(method, float):
			ref_value = closest_value(method,self.columns['x'],self.columns[column])
		else:
			ref_value = method(data)

		array = np.array(data)
		array -= ref_value
		self.columns[column] = list(array)

	def smooth(self,column='y',window_length=24,polyorder=3):
		"""Smooth a column of data"""
		from scipy.signal import savgol_filter

		self.columns[column] = list(savgol_filter(self.columns[column], window_length, polyorder))

class XVGCollection():
	"""Class to store many XVG objects"""
	def __init__(self,xvg_list,title):
		self.children = xvg_list
		self.title = title

	def smooth(self,column='y',window_length=24,polyorder=3):
		"""Smooth all the children's data"""
		for xvg in self.children:
			xvg.smooth(column,window_length,polyorder)

	def map_xdata(self,before_low,before_high,after_low,after_high):

		scale = (after_high - after_low)/(before_high - before_low)

		for xvg in self.children:

			data = xvg.columns['x']
			mapped = []

			for x in data:
				mapped.append((x-before_low)*scale + after_low)
			
			xvg.columns['x'] = mapped

	def plotly(self,show=False,fig=None,statistics=False,stationary_points=False,color=None,group_from_title=False,no_layout=False):
		"""Use plotly to plot the XVG collection"""

		import plotly.graph_objects as go

		if fig is None:
			fig = go.Figure()

		if not no_layout:
			fig.update_layout(
				title=self.title,
				xaxis_title=self.children[0].xlabel,
				yaxis_title=self.children[0].ylabel,
				legend=dict(groupclick="toggleitem"),
				font=dict(family="Helvetica Neue",size=18)
			)

		if not statistics:

			for xvg in self.children:

				fig.add_trace(xvg.get_go_trace(name=xvg.title,color=color,group_from_title=group_from_title))

				if stationary_points:
					fig.add_trace(xvg.get_stationary_point_trace(name=xvg.title,color=color))

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

			if not color:
				color = 'rgb(0,0,0)' # default is black

			fig.add_trace(go.Scatter(x=x,y=maxs,name="max",legendgroup=self.title,legendgrouptitle_text=self.title,line=dict(width=0,color=color)))
			fig.add_trace(go.Scatter(x=x,y=mins,name="min",fill='tonexty',legendgroup=self.title,legendgrouptitle_text=self.title,line=dict(width=0,color=color),fillcolor=color.replace('rgb','rgba').replace(')',',0.15)')))
			fig.add_trace(go.Scatter(x=x,y=mean_plus_std,name="mean+std",legendgroup=self.title,legendgrouptitle_text=self.title,line=dict(width=0,color=color)))
			fig.add_trace(go.Scatter(x=x,y=mean_minus_std,name="mean-std",fill='tonexty',legendgroup=self.title,legendgrouptitle_text=self.title,line=dict(width=0,color=color),fillcolor=color.replace('rgb','rgba').replace(')',',0.3)')))
			fig.add_trace(go.Scatter(x=x,y=mean,name="mean",legendgroup=self.title,legendgrouptitle_text=self.title,line=dict(color=color,width=4)))

			fig.update_layout(showlegend=False)

		if show:
			fig.show()

		return fig

	def align_ydata(self,method):
		"""Shift the y-values of all the traces to align them according to a function such as min/max or an x-coordinate"""

		import numpy as np
		from .signal import closest_value

		for xvg in self.children:

			data = xvg.columns['y']

			if isinstance(method, float):
				ref_value = closest_value(method,xvg.columns['x'],xvg.columns['y'])
			else:
				ref_value = method(data)

			array = np.array(data)
			array -= ref_value
			xvg.columns['y'] = list(array)

### Grace Escape Codes:

# ignore for plotly
_font_codes = {
	"\\f{x}": 'switch to font named "x"',
	"\\f{n}": 'switch to font number n',
	"\\f{}": 'return to original font',	
	"\\x": 'switch to Symbol font (same as \\f{Symbol})',
}

# ignore for plotly
_color_codes = {
	"\\R{x}": 'switch to color named "x"',
	"\\R{n}": 'switch to color number n',
	"\\R{}": 'return to original color',
}

# ignore for plotly
_style_codes = {
	"\\u": 'begin underline',
	"\\U": 'stop underline',
	"\\o": 'begin overline',
	"\\O": 'stop overline',
	"\\q": 'make font oblique (same as \\l{0.25})',
	"\\Q": 'undo oblique (same as \\l{-0.25})',
}

# super/subscripting
_script_codes = {
	"\\s": 'begin subscripting (same as \\v{-0.4}\\z{0.71})',
	"\\S": 'begin superscripting (same as \\v{0.6}\\z{0.71})',
	"\\N": 'return to normal style (same as \\v{}\\t{})',
}

_other_codes = {
	"\\#{x}": 'treat "x" (must be of even length) as list of hexadecimal char codes',
	"\\t{xx xy yx yy}": 'apply transformation matrix',
	"\\t{}": 'reset transformation matrix',
	"\\z{x}": 'zoom x times',
	"\\z{}": 'return to original zoom',
	"\\r{x}": 'rotate by x degrees',
	"\\l{x}": 'slant by factor x',
	"\\v{x}": 'shift vertically by x',
	"\\v{}": 'return to unshifted baseline',
	"\\V{x}": 'shift baseline by x',
	"\\V{}": 'reset baseline',
	"\\h{x}": 'horizontal shift by x',
	"\\n": 'new line',
	"\\Fk": 'enable kerning',
	"\\FK": 'disable kerning',
	"\\Fl": 'enable ligatures',
	"\\FL": 'disable ligatures',
	"\\m{n}": 'mark current position as n',
	"\\M{n}": 'return to saved position n',
	"\\dl": 'LtoR substring direction',
	"\\dr": 'RtoL substring direction',
	"\\dL": 'LtoR text advancing',
	"\\dR": 'RtoL text advancing',
	"\\+": 'increase size (same as \\z{1.19} ; 1.19 = sqrt(sqrt(2)))',
	"\\-": 'decrease size (same as \\z{0.84} ; 0.84 = 1/sqrt(sqrt(2)))',
	"\\T{xx xy yx yy}": 'same as \\t{}\\t{xx xy yx yy}',
	"\\Z{x}": 'absolute zoom x times (same as \\z{}\\z{x})',
	"\\\\": 'print \\',
	"\\n": 'switch to font number n (0-9) (deprecated)',
	"\\c": 'begin using upper 128 characters of set (deprecated)',
	"\\C": 'stop using upper 128 characters of set (deprecated)',
}

_ignorable_codes = _font_codes | _color_codes | _style_codes | _other_codes
_all_codes = _font_codes | _style_codes | _color_codes | _script_codes | _other_codes
