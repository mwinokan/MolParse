
from .residue import Residue, res_type

def alphabet(name):
	match name:
		case "DA": return "A" 
		case "DT": return "T" 
		case "DG": return "G" 
		case "DC": return "C" 
		case "DA3": return "A" 
		case "DT3": return "T" 
		case "DG3": return "G" 
		case "DC3": return "C" 
		case "DA5": return "A" 
		case "DT5": return "T" 
		case "DG5": return "G" 
		case "DC5": return "C" 
		case other:
			import mout
			mout.errorOut(f"Unsupported nucleic acid residue name: {name}")
			return None

def longname(name):
	match name:
		case "DA": return "Adenine" 
		case "DT": return "Thymine" 
		case "DG": return "Guanine" 
		case "DC": return "Cytosine" 
		case "DA3": return "Adenine 3-Terminus" 
		case "DT3": return "Thymine 3-Terminus" 
		case "DG3": return "Guanine 3-Terminus" 
		case "DC3": return "Cytosine 3-Terminus" 
		case "DA5": return "Adenine 5-Terminus" 
		case "DT5": return "Thymine 5-Terminus" 
		case "DG5": return "Guanine 5-Terminus" 
		case "DC5": return "Cytosine 5-Terminus" 
		case other:
			import mout
			mout.errorOut(f"Unsupported nucleic acid residue name: {name}")
			return None

class NucleicAcid(Residue):
	"""Class for Nucleic Acid Residue

	These objects should not be created by the user, 
	but constructed automatically when parsing a 
	coordinate file via amp.parsePDB or otherwise"""

	def __init__(self,name: str,number: int=None,chain: str=None,atoms=None):
		assert res_type(name) == "DNA"
		super(NucleicAcid, self).__init__(name,number,chain,atoms)

		self._backbone = None
		self._nucleobase = None

	@property
	def is_purine(self):
		return self.name in ["DG","DA"]

	@property
	def is_pyrimidine(self):
		return self.name in ["DC","DT"]

	@property
	def letter(self):
		"""Nucleic acid alphabet"""
		return alphabet(self.name)

	@property
	def longname(self):
		"""Nucleic acid chemical name"""
		return longname(self.name)

	@property
	def type(self):
		"""Fixed typing"""
		return "DNA"

	@property
	def backbone(self):
		"""Get backbone atoms
			
			 Y
			 |
		  O==C   HN
		     |   |
	      R--CA--N--X
		     |
		     HA
	
		"""


		if not self._backbone:
			names = ["H5T", "O5'", "C5'",  "H5'", "H5''", 
					 "C4'", "H4'", "O4'",  "C1'", "H1'", 
					 "C2'", "H2'", "H2''", "C3'", "H3'", 
					 "O3'", "P",   "O1P",  "O2P", "H3T"]
			self._backbone = self.get_atom(names)

		return self._backbone

	@property
	def nucleobase(self):
		"""Get nucleobase atoms"""

		if not self._nucleobase:
			bb_names = ["H5T", "O5'", "C5'",  "H5'", "H5''", 
					    "C4'", "H4'", "O4'",  "C1'", "H1'", 
					    "C2'", "H2'", "H2''", "C3'", "H3'", 
					    "O3'", "P",   "O1P",  "O2P", "H3T"]
			names = [n for n in self.atom_names() if n not in bb_names]
			self._nucleobase = self.get_atom(names)

		return self._nucleobase

	def make_5ter(self):
		"""Create a 5-Terminus"""

		# Append 5 to residue name
		if not residue.name.endswith("5"):
			residue.name += "5"

		# Remove HTER/H5T
		residue.delete_atom("HTER")
		residue.delete_atom("H5T")
		residue.delete_atom("HO5'")
		residue.delete_atom("OXT")
		residue.delete_atom("O5T")
		residue.delete_atom("O1P")
		residue.delete_atom("O2P")
		residue.delete_atom("OP1")
		residue.delete_atom("OP2")

		# Rename P->H5T
		atom = residue.get_atom("P")
		if atom is not None:
			atom.name = "H5T"

	def make_3ter(residue):
		"""Create a 3-Terminus"""

		# Append 3 to residue name
		if not residue.name.endswith("3"):
			residue.name += "3"

		# Remove HTER/H3T
		residue.delete_atom("O1P3")
		residue.delete_atom("O2P3")
		residue.delete_atom("O3T")
		residue.delete_atom("H3T")
		residue.delete_atom("HO3'")
		residue.delete_atom("HCAP")

		# Rename P3->H3T
		atom = residue.get_atom("P3")
		if atom is not None:
			atom.name = "H3T"

	def mutate(self,newname):
		""" Mutate this nucleic acid to another"""

		import mout
		import mcol
		mout.headerOut(f"Mutating {mcol.arg}{self.longname}{mcol.clear+mcol.bold} --> {mcol.arg}{longname(newname)}{mcol.clear+mcol.bold} ({mcol.result}{self.letter}{self.number}{alphabet(newname)}{mcol.clear+mcol.bold})")

		if self.name == newname:
			mout.warningOut("Skipping mutation with self == target!")
			return

		import os
		amp_path = os.path.dirname(__file__)

		from .io import parse

		if newname in ["DA","DT"]:
			ref_sys = parse(f"{amp_path}/ref/AT_FLAT.pdb",verbosity=0)
		elif newname in ["DG","DC"]:
			ref_sys = parse(f"{amp_path}/ref/GC_FLAT.pdb",verbosity=0)
		else:
			mout.errorOut("Unsupported target mutation")

		rem = [c for c in ref_sys.chain_names if c not in newname[-1]][0]
		ref_sys.remove_chain(rem,verbosity=0)
		ref_res = ref_sys.residues[0]

		if self.is_purine and ref_res.is_purine:
			print("Purine --> Purine")

			# view_sys = ref_sys.copy()
			# for atom in self.atoms:
			# 	view_sys.add_atom(atom)
			# view_sys.view()

			targets = [self.get_atom(n).np_pos for n in ["N9","C6","C2"]]
			indices = [ref_res.get_atom(n).index for n in ["N9","C6","C2"]]
			ref_sys.align_by_pairs(targets,indices,alt=False)

			view_sys = ref_sys.copy()
			for atom in self.atoms:
				view_sys.add_atom(atom)
			view_sys.view()

		elif self.is_pyrimidine and ref_res.is_pyrimidine:
			print("Pyrimidine --> Pyrimidine")
		elif self.is_purine and ref_res.is_pyrimidine:
			print("Purine --> Pyrimidine")
		elif self.is_pyrimidine and ref_res.is_purine:
			print("Pyrimidine --> Purine")

