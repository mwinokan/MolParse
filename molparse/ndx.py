
from collections import UserDict

def parseNDX(file):

	index_dict = GromacsIndex()

	with open(file) as f:

		group_name = None
		group_indices = None

		for line in f:

			line = line.strip()

			if line.startswith('['):

				if group_name:
					index_dict[group_name] = group_indices
				
				group_name = line.replace('[','').replace(']','').strip()
				group_indices = []

			elif group_name:

				for word in line.split(" "):

					if not word:
						continue

					group_indices.append(int(word.strip()))

			else:
				print(line)

	return index_dict

class GromacsIndex(UserDict):

	def __init__(self):
		super(GromacsIndex, self).__init__()

	@property
	def group_names(self):
		return self.data.keys()

	def shift(self,shift=-1):
		for key in self.data:
			self.data[key] = [i + shift for i in self.data[key]]

	def summary(self):

		import mout
		import mcol

		for group in self.data:
			print(f'{mcol.varType}Group {mcol.arg}{group}{mcol.clear} contains {mcol.result}{len(self.data[group])}{mcol.clear} indices')
