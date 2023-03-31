
class AtomGroup():
	"""General class for a group of atoms. Do not construct this object! Let Molparse handle it"""
	def __init__(self,name: str):

		self._name = name

		# GUI state variables
		self._expand = False
		self._show_context = False
		self._context_options = {}

		from .list import NamedList
		self._atoms = NamedList()

	""" To-Do's:

		- Do not recalculate atomic properties if not necessary

	"""

### PROPERTIES

	@property
	def atoms(self):
		"""Child classes of AtomGroup should overload this method"""
		return self._atoms

	@property
	def children(self):
		"""Child classes of AtomGroup should overload this method"""
		return self.atoms

	@property
	def name(self):
		return self._name

	@name.setter
	def name(self,name: str):
		self._name = name
		self.fix_names()

	@property
	def num_atoms(self):
		return len(self.atoms)

	@property
	def charge(self):
		"""Total charge (float)"""
		charge = 0
		for atom in self.atoms:
			charge += atom.charge
		return charge

	@property
	def atomic_numbers(self):
		"""Get all child Atom atomic numbers (list)"""
		number_list = []
		for atom in self.atoms:
			number_list.append(atom.atomic_number)
		return number_list

	@property
	def positions(self):
		"""Get all child Atom positions (list)"""
		positions_list = []
		for atom in self.atoms:
			positions_list.append(atom.position)
		return positions_list

	@property
	def charges(self):
		"""Get all child Atom charges (list)"""
		charges = []
		for atom in self.atoms:
			charges.append(atom.charge)
		return charges

	@property
	def symbols(self):
		"""Get all child Atom symbols (list)"""
		symbols = []
		for atom in self.atoms:
			symbols.append(atom.symbol)  
		return symbols

	@property
	def masses(self):
		"""Get all child Atom masses (list)"""
		masses = []
		for atom in self.atoms:
			masses.append(atom.mass)  
		return masses

	@property
	def species(self):
		"""Returns species of all child Atoms (list)"""
		species_list = []
		for atom in self.atoms:
			species_list.append(atom.species)
		return ''.join(species_list)

	@property
	def indices(self):
		"""Returns indices of all child Atoms (list)"""
		indices = []
		for atom in self.atoms:
			indices.append(atom.index)
		return indices

	@property
	def bbox(self):
		"""Bounding box of the AtomGroup"""
		import numpy as np
		x = [min([a.np_pos[0] for a in self.atoms]),max([a.np_pos[0] for a in self.atoms])]
		y = [min([a.np_pos[1] for a in self.atoms]),max([a.np_pos[1] for a in self.atoms])]
		z = [min([a.np_pos[2] for a in self.atoms]),max([a.np_pos[2] for a in self.atoms])]
		return [x,y,z]

	@property
	def bbox(self):
		"""Bounding box of the AtomGroup"""
		import numpy as np
		x = [min([a.np_pos[0] for a in self.atoms]),max([a.np_pos[0] for a in self.atoms])]
		y = [min([a.np_pos[1] for a in self.atoms]),max([a.np_pos[1] for a in self.atoms])]
		z = [min([a.np_pos[2] for a in self.atoms]),max([a.np_pos[2] for a in self.atoms])]
		return [x,y,z]

	@property
	def bbox_sides(self):
		"""Bounding box sides of the Atomgroup"""
		import numpy as np
		x = max([a.np_pos[0] for a in self.atoms])-min([a.np_pos[0] for a in self.atoms])
		y = max([a.np_pos[1] for a in self.atoms])-min([a.np_pos[1] for a in self.atoms])
		z = max([a.np_pos[2] for a in self.atoms])-min([a.np_pos[2] for a in self.atoms])
		return [x,y,z]

	@property
	def bbox_norm(self):
		"""Length of bounding box diagonal"""
		import numpy as np
		return np.linalg.norm([x[1]-x[0] for x in self.bbox])

	@property
	def ase_atoms(self):
		"""Construct an equivalent ase.Atoms object"""

		# from .io import write, read
		# write("__temp__.pdb",self,verbosity=0)
		# return read("__temp__.pdb",verbosity=0)
		from ase import Atoms
		return Atoms(symbols=self.symbols,cell=None, pbc=None,positions=self.positions)

	@property
	def covalent_radii(self):
		
		from ase.data import covalent_radii

		radii = []

		for atom in self.atoms:
			radii.append(covalent_radii[atom.atomic_number])

		return radii

### METHODS

	def add_atom(self,atom):
		atom.parent = self
		self._atoms.append(atom)

	def summary(self):
		print(self.name)
		for a in self.atoms:
			print(a,a.position)

	def write(self,filename):
		from .io import write
		write(filename,self)

	# overloaded by child classes
	def fix_names(self):
		pass

	def get_atom_by_index(self,index:int,use_pdb_index:bool=True):
		"""Get Atom by its index"""
		for atom in self.atoms:
			if use_pdb_index:
				if atom.pdb_index == index:
					return atom
			else:
				if atom.index == index:
					return atom

	def centre_of_mass(self,set=None,shift=None):
		"""Calculate centre of mass"""
		return self.CoM(set=set,shift=shift)

	def center(self):
		"""Move the system's CoM to the origin"""
		return self.CoM(set=[0,0,0])

	def CoM(self,set=None,shift=None,verbosity=1):
		"""Calculate or manipulate the system's centre of mass. 

		if not set and not shift: return CoM
		if set: move the System to the CoM
		if shift: move the System by the specified vector"""

		import mout
		import numpy as np

		position_list = self.positions

		centre_of_mass = np.array([sum([pos[0] for pos in position_list])/len(position_list),
			sum([pos[1] for pos in position_list])/len(position_list),
			sum([pos[2] for pos in position_list])/len(position_list)])

		if verbosity > 0: 
			mout.varOut("CoM of "+self.name,
			centre_of_mass,
			unit="Angstroms",precision=4)

		if set is not None:

			assert np.ndim(set) == 1

			try:
				new_com = np.array([set[0],set[1],set[2]])
			except:
				mout.errorOut("Incompatible array input",code="amp.system.CoM.1")

			shift_com = new_com - centre_of_mass

		else:

			shift_com = np.array([0,0,0])

		if shift is not None:

			assert np.ndim(shift) == 1

			try:
				new_shift = np.array([shift[0],shift[1],shift[2]])
			except:
				mout.errorOut("Incompatible array input",code="amp.system.CoM.2")

			shift_com = shift_com + new_shift

		if set is not None or shift is not None:

			for atom in self.atoms:
				atom.position = atom.position + shift_com

		return centre_of_mass

	def guess_bonds(self,scale=1.2):
		"""Guess connectivity using covalent radii (list of pairs)"""
		import numpy as np
		from ase.gui.view import get_bonds
		ase_atoms = self.ase_atoms
		radii = np.array(self.covalent_radii)*(scale/1.5)
		bonds = get_bonds(ase_atoms,radii)
		return [[b[0],b[1]] for b in bonds]

	def guess_names(self,target):
		"""Try and set the atom names of the system by looking for 
		the closest atom of the same species in the target system"""

		import numpy as np

		positions = [b.np_pos for b in target.atoms]
		species = [b.species for b in target.atoms]

		for a in self.atoms:
			a_pos = a.np_pos
			distances = [np.linalg.norm(a_pos - p) if s == a.species else 999 for s,p in zip(species,positions)]
			index = np.argmin(distances)
			b = target.atoms[index]
			a.set_name(b.name,verbosity=0)

	def plot3d(self,extra=[],alpha=1.0,bonds=True):
		"""Render the atoms with plotly graph objects. 
		extra can contain pairs of coordinates to be shown as vectors."""

		if bonds:
			bonds = self.guess_bonds()
			bond_vectors = []
			for a,b in bonds:
				bond_vectors.append([self.atoms[a].np_pos,self.atoms[b].np_pos])
			bonds = bond_vectors
		else:
			bonds = []

		from .go import plot3d
		return plot3d(self.atoms,extra,bonds,alpha)

	def set_coordinates(self,reference):
	    """Set all coordinates according to a reference ase.Atoms object"""
	    if type(reference) is str:
	      from ase.io import read
	      atoms = read(reference)
	    elif isinstance(reference,list):
	      for index,atom in enumerate(self.atoms):
	        atom.position = reference[index]
	      return
	    else:
	      atoms = reference

	    for index,atom in enumerate(self.atoms):
	      atom.position = atoms[index].position

	def rotate(self,angle,vector,center=(0,0,0)):
		"""Rotate the system (see ase.Atoms.rotate)"""
		ase_atoms = self.ase_atoms
		ase_atoms.rotate(angle,vector,center=center)
		self.set_coordinates(ase_atoms)

	def align_by_posmap(self,map):

	    """Align the system to a target by superimposing three shared atoms:

	      a --> A
	      b --> B
	      c --> C

	      (a,b,c) are atoms from this system
	      (A,B,C) are atoms from the target system

	      map should contain Atoms: [[a,b,c],[A,B,C]]

	    """

	    import numpy as np

	    a = map[0][0]
	    b = map[0][1]
	    c = map[0][2]

	    A = map[1][0]
	    B = map[1][1]
	    C = map[1][2]

	    # TRANSLATION

	    displacement = A - a
	    self.translate(displacement)

	    # ROTATION 1

	    d = (b+c)/2
	    D = (B+C)/2
	    self.rotate(d-A.np_pos,D-A.np_pos,center=A.np_pos)

	    # ROTATION 2

	    d = (b+c)/2
	    D = (B+C)/2

	    v_bc = c - b
	    v_BC = C - B

	    def unit_vector(a):
	      return a/np.linalg.norm(a)

	    def remove_parallel_component(vector,reference):
	      unit_reference = unit_vector(reference)
	      projection = np.dot(vector,unit_reference)
	      return vector - projection

	    v_ad = d - a.np_pos
	    v_bc_normal_to_ad = remove_parallel_component(v_bc, v_ad)
	    v_BC_normal_to_ad = remove_parallel_component(v_BC, v_ad)
	    self.rotate(v_bc_normal_to_ad, v_BC_normal_to_ad,center=A.np_pos)

	    # extra stuff for plotly
	    extra = [[A.np_pos,A.np_pos+v_ad]]
	    return extra

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

	@property
	def _context_info(self):

		com = self.CoM(verbosity=0)
		com = f'[{com[0]:.2f} {com[1]:.2f} {com[2]:.2f}]'

		bbox = self.bbox
		bbox = f'[[{bbox[0][0]:.2f} {bbox[0][1]:.2f}] [{bbox[1][0]:.2f} {bbox[1][1]:.2f}] [{bbox[2][0]:.2f} {bbox[2][1]:.2f}]]'

		bbox_sides = self.bbox_sides
		bbox_sides = f'[{bbox_sides[0]:.2f} {bbox_sides[1]:.2f} {bbox_sides[2]:.2f}]'

		items = {
			'Centre of Mass': com,
			'Bounding Box': bbox,
			'Bounding Box Sides': bbox_sides,
			'Bounding Box Diagonal': f'{self.bbox_norm:.2f}',
			'Total Charge': self.charge,
			'Total #Atoms': self.num_atoms,
		}

		from .system import System
		from .group import AtomGroup
		if not isinstance(self,System) and not isinstance(self,AtomGroup):
			file_out = f'{self.__class__.__name__}_{self.name}_{self.index}.pdb'

			clickables = {
				f'Export PDB "{file_out}"': lambda x: x.write(file_out),
			}

		else:
			clickables = {}

		if self.num_atoms == 1:
			items.pop('Bounding Box')
			items.pop('Bounding Box Sides')
			items.pop('Bounding Box Diagonal')

		return items, clickables

### DUNDERS

	def __repr__(self):
		return self.name
	def __str__(self):
		return self.name
	def __len__(self):
		return len(self.children)

	def __getitem__(self,key):

		"""Access child Atom, Residue or Chain by index or name

		group['a0'] : returns the first atom
		group['an10'] : returns the residue with number/resid 10
		group['aH1'] : returns all atoms named 'H1'

		group['r0'] : returns the first residue
		group['rn10'] : returns the residue with number/resid 10
		group['rASN'] : returns all residues named 'ASN'

		group['c0'] : returns the first chain
		group['cA'] : returns all chains named 'A'
		"""

		try:
			
			if key[0] == 'a':
				return self.atoms[key[1:]]

			elif key[0] == 'r':
				return self.residues[key[1:]]

			elif key[0] == 'c':
				return self.chains[key[1:]]

			else:
				import mout
				mout.errorOut("key must start with 'a', 'r', or 'c' to refer to an Atom, Residue, or Chain, respectively. See help(molparse.group.__getitem__)")

		except AttributeError:

			import mout
			lookup = {'a':'atoms','c':'chains','r':'residues'}
			mout.errorOut(f"{type(self)} does not possess {lookup[key[0]]} attribute")
