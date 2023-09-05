
from collections import UserDict

def parseNDX(file, verbosity=1):

	import mout
	import mcol

	if verbosity:
		mout.out(f'parsing {mcol.file}{file}{mcol.clear} ... ', end='')

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
			
			if group_name:
				index_dict[group_name] = group_indices

	if verbosity:
		mout.out(f'Done.')

	return index_dict

def writeNDX(file, ndx, chunk_size=15, verbosity=1):

	import mout
	import mcol

	if verbosity:
		mout.out(f'writing {mcol.file}{file}{mcol.clear} ... ', end='')

	with open(file,'wt') as f:

		for group_name in ndx:

			f.write(f'[ {group_name} ]\n')

			indices = ndx[group_name]

			if len(indices) < chunk_size:
				indices = [str(i) for i in indices]
				f.write(f'{" ".join(indices)}\n')

			else:
				for index_chunk in grouper(indices, chunk_size, incomplete='ignore'):
					index_chunk = [str(i) for i in index_chunk]
					f.write(f'{" ".join(index_chunk)}\n')

	if verbosity:
		mout.out(f'Done.')

# https://docs.python.org/3/library/itertools.html#itertools-recipes
def grouper(iterable, n, *, incomplete='fill', fillvalue=None):
    """Collect data into non-overlapping fixed-length chunks or blocks
	   
	* grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
	* grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
	* grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
	   
    """

    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return zip_longest(*args, fillvalue=fillvalue)
    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')

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

	def write(self,file):
		writeNDX(file,self)
