
class AtomGroup():
	"""General class for a group of atoms. Do not construct this object! Let Molparse handle it"""
	def __init__(self,name: str):

		self._name = name

		# ligand
		self._smiles = None

		# GUI state variables
		self._expand = False
		self._show_context = False
		self._context_options = {}

		self._header_data = []

		from .list import NamedList
		self._atoms = NamedList()

	""" To-Do's:

		- Do not recalculate atomic properties if not necessary

	"""

### FACTORIES

	@classmethod
	def from_any(cls,name,source):

		"""Construct a new AtomGroup from:
			- an mp.NamedList of mp.Atom objects (i.e. mp.System.atoms)
			- an mp.System object
			- an mp.AtomGroup object
			- an mp.Chain object
			- an mp.Residue object
			- a list of ase.Atom objects
			- an ase.Atoms object
			- a list containing any combination of ase.Atom(s), mp.Atom, mp.Chain, mp.Residue, objects
			- a string containing a PDB Block
		"""
		
		import mout

		from ase import Atoms as ase_Atoms
		from ase import Atom as ase_Atom

		from .system import System
		from .chain import Chain
		from .residue import Residue
		from .list import NamedList

		# string (PDB Block)
		if isinstance(source,str):
			return cls.from_pdb_block(source)

		# mp.AtomGroup subclass
		elif issubclass(source.__class__,cls):
			return cls.from_group_subclass(source)

		# mp.AtomGroup
		elif isinstance(source, AtomGroup):
			return cls.from_group(source)

		# NamedList of atoms (i.e. mp.System.atoms):
		elif isinstance(source, NamedList) or isinstance(source, ase_Atoms):
			return cls.from_atoms(name,source)

		# list of various objects
		elif isinstance(source, list):

			from .atom import Atom

			group = cls.__new__(cls)
			group.__init__(name)

			for item in source:

				if issubclass(item.__class__,cls):
					group = cls.from_group_subclass(item,group=group)

				elif isinstance(item, NamedList) or isinstance(item,ase_Atoms):
					group = cls.from_atoms(name,item,group=group)

				elif isinstance(item,Atom):
					group.add_atom(item)

				elif isinstance(item,ase_Atom):
					group.add_atom(Atom(name=item.symbol,position=item.position))

				else:
					mout.errorOut(f"Not supported ({item.__class__=})",fatal=True)

			return group

		else:
			mout.errorOut(f"Not supported ({source.__class__=})",fatal=True)

	@classmethod
	def from_pdb_block(cls,pdb_block):

		from .io import parsePDBAtomLine
		
		name = None
		atom_lines = []

		for line in pdb_block.split('\n'):
			
			if line.startswith('COMPND'):
				name = line.split()[1]
				continue

			if line.startswith('HETATM') or line.startswith('ATOM'):
				atom_lines.append(line)
				continue

		atoms = []
		for i,line in enumerate(atom_lines):
			atom = parsePDBAtomLine(line,0,i,0)
			atoms.append(atom)

		group = cls.from_atoms(name,atoms)

		return group

	@classmethod
	def from_atoms(cls,name,atoms,group=None):

		import mout

		from .list import NamedList
		from ase import Atoms as ase_Atoms
		assert isinstance(atoms, list) or isinstance(atoms, NamedList) or isinstance(atoms, ase_Atoms)

		# create new object
		if group is None:
			group = cls.__new__(cls)
			group.__init__(name)

		from .atom import Atom
		from .residue import Residue
		from ase import Atom as ase_Atom

		for atom in atoms:

			if isinstance(atom, Atom):
				group.add_atom(atom.copy())
			
			elif isinstance(atom, ase_Atom):
				group.add_atom(Atom(name=atom.symbol,position=atom.position))

			elif isinstance(atom, Residue):
				for a in atom.atoms:
					group.add_atom(a.copy())

			else:
				mout.errorOut("item in named list is neither mp.Atom, ase.Atom, mp.Residue",fatal=True)

		return group

	@classmethod
	def from_group_subclass(cls,source,group=None):

		assert issubclass(source.__class__,cls)

		if group is None:
			group = cls.__new__(cls)
			group.__init__(source.name)

		for atom in source.atoms:
			group.add_atom(atom.copy())

		return group

	@classmethod
	def from_group(cls,group):
		assert isinstance(group, AtomGroup)
		return group.copy()

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
		return [atom.atomic_number for atom in self.atoms]

	@property
	def positions(self):
		"""Get all child Atom positions (list)"""
		return [atom.np_pos for atom in self.atoms]

	@property
	def charges(self):
		"""Get all child Atom charges (list)"""
		return [atom.charge for atom in self.atoms]

	@property
	def symbols(self):
		"""Get all child Atom symbols (list)"""
		return [atom.symbol for atom in self.atoms]

	@property
	def present_symbols(self):
		return set(self.symbols)

	@property
	def masses(self):
		"""Get all child Atom masses (list)"""
		masses = []
		for atom in self.atoms:
			masses.append(atom.mass)  
		return masses

	@property
	def species(self):
		"""Returns species of all child Atoms (str)"""
		return ''.join([atom.species for atom in self.atoms])

	@property
	def indices(self):
		"""Returns indices of all child Atoms (list)"""
		return [atom.index for atom in self.atoms]

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
		atoms = self.atoms
		x_coords = [a.np_pos[0] for a in atoms]
		y_coords = [a.np_pos[1] for a in atoms]
		z_coords = [a.np_pos[2] for a in atoms]
		x = [min(x_coords),max(x_coords)]
		y = [min(y_coords),max(y_coords)]
		z = [min(z_coords),max(z_coords)]
		return [x,y,z]

	@property
	def bbox_center(self):
		"""Center of the AtomGroups' bounding box"""
		return self._bbox_center()

	@property
	def bbox_sides(self):
		"""Bounding box sides of the Atomgroup"""
		return self._bbox_sides()

	@property
	def bbox_norm(self):
		"""Length of bounding box diagonal"""
		return self._bbox_norm()

	@property
	def ase_atoms(self):
		"""Construct an equivalent ase.Atoms object"""

		from ase import Atoms
		return Atoms(symbols=self.symbols,cell=None, pbc=None,positions=self.positions)

	@property
	def covalent_radii(self):
		"""Get a list of covalent radii"""
		
		from ase.data import covalent_radii

		radii = []

		for atom in self.atoms:
			radii.append(covalent_radii[atom.atomic_number])

		return radii

	@property
	def contains_alternative_sites(self):
		"""Return true if any child atoms have an alternative site defined"""
		return any([a.alternative_site for a in self.atoms])

	@property
	def pdb_block(self):
		from .io import constructPDBAtomLine
		
		str_buffer = []
		for atom in self.atoms:
			str_buffer.append(constructPDBAtomLine(atom,atom.number,alt_sites=False))
		return ''.join(str_buffer)

	@property
	def pdb_block_with_alt_sites(self):
		from .io import constructPDBAtomLine
		
		str_buffer = []
		for atom in self.atoms:
			str_buffer.append(constructPDBAtomLine(atom,atom.number,alt_sites=True))
		return ''.join(str_buffer)

	@property
	def rdkit_mol(self):
		from .rdkit import mol_from_pdb_block
		return mol_from_pdb_block(self.pdb_block)

	@property
	def smiles(self):
		if self._smiles is None:
			from .rdkit import mol_to_smiles
			self._smiles = mol_to_smiles(self.rdkit_mol)
		return self._smiles
	
### INTERNAL METHODS

	def _bbox_center(self,bbox=None):
		import numpy as np
		if bbox is None:
			bbox = self.bbox
		return np.array([np.mean(bbox[0]),np.mean(bbox[1]),np.mean(bbox[2])])

	def _bbox_norm(self,bbox=None):
		import numpy as np
		if bbox is None:
			bbox = self.bbox
		return np.linalg.norm([x[1]-x[0] for x in bbox])

	def _bbox_sides(self,bbox=None):
		if bbox is None:
			bbox = self.bbox
		return [bbox[0][1]-bbox[0][0],bbox[1][1]-bbox[1][0],bbox[2][1]-bbox[2][0]]

### METHODS

	def add_atom(self,atom):
		import mout
		atom.parent = self
		if atom.chain is None:
			atom.chain = self._atoms[-1].chain
			mout.warningOut('Taking atom chain from last atom in group!')
		if atom.res_number is None:
			if len(self._atoms):
				atom._res_number = self._atoms[-1].res_number
				mout.warningOut('Taking atom res_number from last atom in group!')
		self._atoms.append(atom)

	def summary(self):
		import mcol
		import mout
		mout.header(f'{mcol.varType}AtomGroup{mcol.clear}{mcol.bold}: {mcol.func}{self.name}{mcol.clear}')

		mout.out(f'{mcol.underline}Sy name {"index":>6} {"number":>6} {"res":>4} {"res#":>6} {"chain":<5}')

		for a in self.atoms:
			mout.out(f'{mcol.varName}{a.name:<2}{mcol.clear}',end=' ')
			mout.out(f'{mcol.varName}{a.name:<4}{mcol.clear}',end=' ')
			mout.out(f'{a.index:>6}',end=' ')
			mout.out(f'{a.number:>6}',end=' ')
			mout.out(f'{a.residue:>4}',end=' ')
			mout.out(f'{a.res_number:>6}',end=' ')
			mout.out(f'{a.chain:<5}',end='\n')

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

	def CoM(self,set=None,shift=None,verbosity=0):
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

	def plot(self,ax=None,color=None,center_index=0,show=False,offset=None,padding=1,zeroaxis=None,frame=False,labels=True,textdict={"horizontalalignment":"center","verticalalignment":"center"}):
		"""Render the system with matplotlib"""

		import numpy as np
		from ase.visualize.plot import plot_atoms

		# copy the system so as not to alter it
		copy = self.copy()

		# create new axes if none provided
		if ax is None:
			import matplotlib.pyplot as plt
			fig,ax = plt.subplots()

		# create colours list
		if color is not None:
			color = [color for a in copy.atoms]

		# create the canvas unit size from bounding box
		canvas = np.array(self.bbox_sides)*1.5

		# center the render if needed
		if offset is None:
			offset=[0,0]
			if center_index is not None:
				vec = - copy.atoms[center_index].np_pos + canvas
				copy.CoM(shift=vec,verbosity=0)

		# use ASE to render the atoms
		plot_atoms(copy.ase_atoms,ax,colors=color,offset=offset,bbox=[0,0,canvas[0]*2,canvas[1]*2])

		# do the labels
		if labels:
			for atom in copy.atoms:
				if atom.species != 'H':
					ax.text(atom.position[0],atom.position[1],atom.name,**textdict)

		# crop the plot
		xmax = copy.bbox[0][1] + padding
		xmin = copy.bbox[0][0] - padding
		ymax = copy.bbox[1][1] + padding
		ymin = copy.bbox[1][0] - padding
		ax.set_xlim(xmin,xmax)
		ax.set_ylim(ymin,ymax)

		# render the centering axes
		if zeroaxis is not None:
			ax.axvline(canvas[0],color=zeroaxis)
			ax.axhline(canvas[1],color=zeroaxis)

		if not frame:
			ax.axis('off')

		if show:
			plt.show()

		return ax, copy

	def plot3d(self,extra=[],alpha=1.0,bonds=True,atoms=True,velocity=False,features=False,v_scale=1.0,fig=None,flat=False,show=True,transform=None):
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

		if not isinstance(features,list) and features:
			from .rdkit import features_from_group
			features = features_from_group(self)

		from .go import plot3d
		return plot3d(self.atoms,extra,bonds,alpha,plot_atoms=atoms,features=features,velocity=velocity,v_scale=v_scale,fig=fig,flat=flat,show=show,transform=transform,title=self.name)

	def set_coordinates(self,reference,velocity=False):
		"""Set all coordinates according to a reference ase.Atoms object"""
		if type(reference) is str:
			if velocity:
				from .io import parse
				sys = parse(reference)
				atoms = sys.atoms
			else:
				from ase.io import read
				atoms = read(reference)
		elif isinstance(reference,list):
			if velocity:
				mout.errorOut("Not supported (group.set_coordinates)",fatal=True)
			for index,atom in enumerate(self.atoms):
				atom.position = reference[index]
			return
		else:
			atoms = reference

		for index,atom in enumerate(self.atoms):
			atom.position = atoms[index].position
			if velocity:
				atom.velocity = atoms[index].velocity

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

	def copy(self):
		import copy
		return copy.deepcopy(self)

	def get_nearby(self,candidates,cutoff):

		"""Get a subset of <candidates> that have at least one atom within <cutoff> of any atom in this AtomGroup."""

		import numpy as np

		# precompute properties of self
		self_atoms = self.atoms
		self_bbox = self.bbox
		self_bbox_halfnorm = self._bbox_norm(self_bbox)/2
		self_bbox_center = self._bbox_center(self_bbox)

		def is_intersecting():

			for self_atom in self_atoms:

				d = np.linalg.norm(self_atom.np_pos - group_bbox_center)

				# skip if atom is not within bbox_halfnorm + cutoff of the group
				if d > cutoff + group_bbox_halfnorm:
					continue

				self_atom_np_pos = self_atom.np_pos

				for group_atom in group.atoms:

					d = np.linalg.norm(self_atom_np_pos - group_atom.np_pos)

					# actually test for inter-atomic distance
					if d <= cutoff:
						return True

		# go through list of candidates
		nearby = []
		for group in candidates:

			# precompute properties of self
			group_bbox = group.bbox
			group_bbox_halfnorm = group._bbox_norm(group_bbox)/2
			group_bbox_center = group._bbox_center(group_bbox)

			center_delta = np.linalg.norm(self_bbox_center - group_bbox_center)
			if center_delta == 0.0:
				continue
			min_bbox_separation = center_delta - self_bbox_halfnorm - group_bbox_halfnorm

			# discard any candidate where the boundingboxes are definitely more than cutoff apart
			if min_bbox_separation > cutoff:
				continue

			# check for intersection
			intersects = is_intersecting()

			if intersects:
				nearby.append(group)

		return nearby

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

		if not isinstance(key,str):
			import mout
			mout.error(f'Cannot index an AtomGroup with a non-string key. See help(molparse.AtomGroup.__getitem__)')
			raise IndexError(f'mp.AtomGroup.__getitem__ received non-string key')

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

	# def __del__(self):
	# 	import mout
	# 	mout.debug(f'{self} was deleted')
