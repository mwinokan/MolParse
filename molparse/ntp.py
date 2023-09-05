
from .residue import Residue
from .group import AtomGroup

ALPHABET = {
	"ATP":"A",
	"TTP":"T",
	"GTP":"G",
	"CTP":"C",
}

LONGNAME = {
	"ATP": "Adenine Tri-Phosphate",
	"TTP": "Thymine Tri-Phosphate",
	"GTP": "Guanine Tri-Phosphate",
	"CTP": "Cytosine Tri-Phosphate",
}

class NucleobaseTriPhosphate(Residue):
	"""Class for Nucleobase Tri Phosphate Residue

	These objects should not be created by the user, 
	but constructed automatically when parsing a 
	coordinate file via amp.parsePDB or otherwise"""

	def __init__(self,name: str,index: int=None, number: int=None,chain: str=None,atoms=None):
		super(NucleobaseTriPhosphate, self).__init__(name,index,number,chain,atoms)

		self._backbone = None
		self._nucleobase = None

	@property
	def longname(self):
		return LONGNAME[self.name]
	
	@property
	def letter(self):
		return ALPHABET[self.name]

	@property
	def is_purine(self):
		return self.name in ["GTP","ATP"]

	@property
	def is_pyrimidine(self):
		return self.name in ["CTP","TTP"]
	
	@property
	def type(self):
		"""Fixed typing"""
		return "LIG"

	@property
	def backbone(self):
		"""Get backbone atoms"""
		if not self._backbone:
			self._backbone = self.get_atom(self.bb_names)
		return self._backbone

	@property
	def nucleobase(self):
		"""Get nucleobase atoms"""
		if not self._nucleobase:
			self._nucleobase = self.get_atom(self.nonbb_names)
		return self._nucleobase

	@property
	def bb_names(self):
		return [
			"C1'",
			"H1'",
			"C2'",
			"O2'",
			"H2'",
			"H2''",
			"C3'",
			"H3'",
			"O3'",
			"H3T",
			"HO3'",
			"C4'",
			"H4'",
			"O4'",
			"C5'",
			"H5'",
			"H5''",
			"O5'",
			"PA",
			"O1A",
			"O2A",
			"O3A",
			"PB",
			"O1B",
			"O2B",
			"O3B",
			"PG",
			"O1G",
			"O2G",
			"O3G",
		]

	@property
	def nonbb_names(self):
		"""Get nucleobase atoms"""
		return [n for n in self.atom_names() if n not in self.bb_names]

	def flip(self):

		import mout
		import mcol
		mout.header(f"Flipping {mcol.arg}{self.longname} {mcol.clear}{mcol.bold}({mcol.arg}{self.name_number_str}{mcol.clear}{mcol.bold})")

		nucleobase_group = AtomGroup.from_any(f'{self.name}.nucleobase',self.nucleobase)

		if self.is_purine:
			source = nucleobase_group.atoms['N9'][0]
			target = nucleobase_group.atoms['O6'][0]

		vector = target - source

		nucleobase_group.rotate(180,vector,source.np_pos)

		for a1,a2 in zip(self.nucleobase, nucleobase_group.atoms):
			a1.position = a2.position

	def mutate(self,newname,show=False):
		""" Mutate this nucleic triphosphate to another"""

		import mout
		import mcol
		mout.header(f"Mutating {mcol.arg}{self.longname}{mcol.clear+mcol.bold} --> {mcol.arg}{LONGNAME[newname]}{mcol.clear+mcol.bold} ({mcol.result}{self.letter}{self.number}{ALPHABET[newname]}{mcol.clear+mcol.bold})")

		if self.name == newname:
			mout.warning("Skipping mutation with self == target!")
			return

		NEW_ANALOGUE = {
			"ATP":"DA",
			"TTP":"DT",
			"GTP":"DG",
			"CTP":"DC",
		}

		import os
		amp_path = os.path.dirname(__file__)

		from .io import parse

		if newname in ["ATP","TTP"]:
			ref_sys = parse(f"{amp_path}/ref/AT_FLAT.pdb",verbosity=0)
		elif newname in ["GTP","CTP"]:
			ref_sys = parse(f"{amp_path}/ref/GC_FLAT.pdb",verbosity=0)
		else:
			mout.error("Unsupported target mutation")

		rem = [c for c in ref_sys.chain_names if c not in NEW_ANALOGUE[newname][-1]][0]
		ref_sys.remove_chain(rem,verbosity=0)
		ref_sys.fix_indices()
		ref_res = ref_sys.residues[0]

		if self.is_purine and ref_res.is_purine:

			# get atoms to align by
			names = ["N9","C6","C2"]
			start = ref_res.get_atom(names)
			target = self.get_atom(names)

			# align the system
			extra = ref_sys.align_by_posmap([start,target])

			# delete intersecting atoms
			ref_res.delete_atom("H9",verbosity=0)

		elif self.is_pyrimidine and ref_res.is_pyrimidine:
			
			# get atoms to align by
			names = ["C6","C4","C2"]
			start = ref_res.get_atom(names)
			target = self.get_atom(names)

			# align the system
			extra = ref_sys.align_by_posmap([start,target])

			# shift to connect to the backbone
			ref_res.translate(self.get_atom("N1")-ref_res.get_atom("N1"))

			# delete intersecting atoms
			ref_res.delete_atom("H1",verbosity=0)

		elif self.is_purine and ref_res.is_pyrimidine:
			mout.warning("Purine --> Pyrimidine")

			# get atoms to align by
			start = ref_res.get_atom(["C6","C4","C2"])
			target = self.get_atom(["C4","C6","C2"])

			# align the system
			extra = ref_sys.align_by_posmap([start,target])

			# shift to connect to the backbone
			ref_res.translate(self.get_atom("N9")-ref_res.get_atom("N1"))

			# delete intersecting atoms
			ref_res.delete_atom("H1",verbosity=0)

		elif self.is_pyrimidine and ref_res.is_purine:
			mout.warning("Pyrimidine --> Purine")

			# get atoms to align by
			start = ref_res.get_atom(["C4","C6","C2"])
			target = self.get_atom(["C6","C4","C2"])

			# align the system
			extra = ref_sys.align_by_posmap([start,target])

			# shift to connect to the backbone
			ref_res.translate(self.get_atom("N1")-ref_res.get_atom("N9"))

			# delete intersecting atoms
			ref_res.delete_atom("H9",verbosity=0)
		
		self.name = newname

		if show:
			view_sys = ref_sys.copy()
			for atom in self.atoms:
				view_sys.add_atom(atom)
			view_sys.plot3d(extra,1.0)

		# delete intersecting atoms
		self.delete_atom(self.nonbb_names,verbosity=0)

		# add the reference nucleobase
		for atom in ref_res.atoms:
			self.add_atom(atom)

	def oxidise(self):

		assert self.name == 'GTP'

		import numpy as np
		from .atom import Atom

		import mout
		import mcol

		mout.header(f"Oxidising {mcol.arg}{self.longname}{mcol.header} --> {mcol.result}{mcol.bold}8-Oxo-7,8-dihydroguanine triphosphate (OGTP)")

		# H8 --> O8
		self.get_atom('H8').name = 'O8'

		# +H7
		N7 = self.get_atom('N7')
		N9 = self.get_atom('N9')
		C4 = self.get_atom('C4')
		vec = N7 - (N9 + C4)/2
		vec /= np.linalg.norm(vec)
		atom = Atom("H7")
		atom.position = N7 + vec
		self.add_atom(atom)

		# +H2
		N2 = self.get_atom('N2')
		C2 = self.get_atom('C2')
		N3 = self.get_atom('N3')
		vec = np.cross(N2 - C2, N3 - C2)
		vec /= np.linalg.norm(vec)
		atom = Atom("H2")
		atom.position = C2 + vec
		self.add_atom(atom)

		# +H3
		C6 = self.get_atom('C6')
		vec = N3 - C6
		vec /= np.linalg.norm(vec)
		atom = Atom("H3")
		atom.position = N3 + vec
		self.add_atom(atom)

		self.name = 'OGTP'
