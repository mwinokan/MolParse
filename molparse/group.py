
class AtomGroup():
	"""General class for a group of atoms. Do not construct this object! Let Molparse handle it"""
	def __init__(self,name: str):

		self._name = name

		# GUI state variables
		self._expand = False
		self._show_context = False
		self._context_options = {}

	### PROPERTIES

	# overloaded by child classes
	def fix_names(self):
		pass

	@property
	def name(self):
		return self._name

	@name.setter
	def name(self,name: str):
		self._name = name
		self.fix_names()

	### GUI THINGS

	# open GUI tree viewer
	def tree(self):
		from .tree import tree
		tree(self)

	# expand GUI tree view
	def expand(self):
		self._expand = True
		for child in self.children:
			child._expand = False

	# collapse GUI tree view
	def collapse(self):
		self._expand = False

	# open GUI context window
	def spawn_context(self,line,col):
		self._tree.log(f"trying to spawn context window at {line} {col}")
		self._show_context = True
		self._tree.spawn_context(self,line,col)

	# close GUI context window
	def hide_context(self):
		self._tree.log(f"trying to hide context window")
		self._show_context = False
		self._tree.hide_context_menu()

	### DUNDERS

	def __repr__(self):
		return self.name
	def __str__(self):
		return self.name
	def __len__(self):
		return len(self.children)
