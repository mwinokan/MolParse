
from .residue import Residue, res_type

def alphabet(name):
	if "ALA": return "A" 
	elif name == "ARG": return "R" 
	elif name == "ASN":	return "N" 
	elif name == "ASP": return "D" 
	elif name == "CYS": return "C" 
	elif name == "GLN": return "Q" 
	elif name == "GLU": return "E" 
	elif name == "GLY": return "G" 
	elif name == "HIS": return "H" 
	elif name == "HSD": return "H" 
	elif name == "ILE": return "I" 
	elif name == "LEU": return "L" 
	elif name == "LYS": return "K" 
	elif name == "MET": return "M" 
	elif name == "PHE": return "F" 
	elif name == "PRO": return "P" 
	elif name == "SER": return "S" 
	elif name == "THR": return "T" 
	elif name == "TRP": return "W" 
	elif name == "TYR": return "Y" 
	elif name == "VAL": return "V" 
	else:
		import mout
		mout.errorOut(f"Unsupported amino acid residue name: {name}")
		return None

def longname(name):
	if name == "ALA": return "Alanine" 
	elif name == "ARG": return "Arginine" 
	elif name == "ASN":	return "Asparagine" 
	elif name == "ASP": return "Aspartic Acid" 
	elif name == "CYS": return "Cysteine" 
	elif name == "GLN": return "Glutamine" 
	elif name == "GLU": return "Glutamic Acid" 
	elif name == "GLY": return "Glycine" 
	elif name == "HIS": return "Histidine" 
	elif name == "HSD": return "Histidine" 
	elif name == "ILE": return "Isoleucine" 
	elif name == "LEU": return "Leucine" 
	elif name == "LYS": return "Lysine" 
	elif name == "MET": return "Methionine" 
	elif name == "PHE": return "Phenylalanine" 
	elif name == "PRO": return "Proline" 
	elif name == "SER": return "Serine" 
	elif name == "THR": return "Threonine" 
	elif name == "TRP": return "Tryptophan" 
	elif name == "TYR": return "Tyrosine" 
	elif name == "VAL": return "Valine" 
	else:
		import mout
		mout.errorOut(f"Unsupported amino acid residue name: {name}")
		return None

class AminoAcid(Residue):
	"""Class for Amino Acid Residue

	These objects should not be created by the user, 
	but constructed automatically when parsing a 
	coordinate file via amp.parsePDB or otherwise"""

	def __init__(self,name: str,number: int=None,chain: str=None,atoms=None):
		assert res_type(name) == "PRO"
		super(AminoAcid, self).__init__(name,number,chain,atoms)

		self._backbone = None
		self._sidechain = None

	@property
	def letter(self):
		"""Amino acid alphabet"""
		return alphabet(self.name)

	@property
	def longname(self):
		"""Amino acid chemical name"""
		return longname(self.name)

	@property
	def type(self):
		"""Fixed typing"""
		return "PRO"
	
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
			names = ["N","HN","CA","HA","C","O"]
			self._backbone = self.get_atom(names)

		return self._backbone

	@property
	def sidechain(self):
		"""Get sidechain atoms"""

		if not self._sidechain:
			bb_names = ["N","HN","CA","HA","C","O"]
			names = [n for n in self.atom_names() if n not in bb_names]
			self._sidechain = self.get_atom(names)

		return self._sidechain

	@property
	def CA(self):
		"""Get Carbon-Alpha Atom"""
		return self.get_atom("CA")

	@property
	def CB(self):
		"""Get Carbon-Alpha Atom"""
		return self.get_atom("CB")

	@property
	def HA(self):
		"""Get Carbon-Alpha Hydrogen Atom"""
		return self.get_atom("HA")

	def mutate(self,newname):
		"""Mutate this amino acid to another"""
		import mout
		import mcol
		mout.headerOut(f"Mutating {mcol.arg}{self.longname}{mcol.clear+mcol.bold} --> {mcol.arg}{longname(newname)}{mcol.clear+mcol.bold} ({mcol.result}{self.letter}{self.number}{alphabet(newname)}{mcol.clear+mcol.bold})")
		if newname ==  "ALA":
			from .mutate.protein import to_alanine
			to_alanine(self)
		else:
			import mout
			mout.errorOut("Unsupported mutation!")
	