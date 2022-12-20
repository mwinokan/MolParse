
from .residue import Residue, res_type

def alphabet(name):
	match name:
		case "ALA": return "A" 
		case "ARG": return "R" 
		case "ASN":	return "N" 
		case "ASP": return "D" 
		case "CYS": return "C" 
		case "GLN": return "Q" 
		case "GLU": return "E" 
		case "GLY": return "G" 
		case "HIS": return "H" 
		case "HSD": return "H" 
		case "ILE": return "I" 
		case "LEU": return "L" 
		case "LYS": return "K" 
		case "MET": return "M" 
		case "PHE": return "F" 
		case "PRO": return "P" 
		case "SER": return "S" 
		case "THR": return "T" 
		case "TRP": return "W" 
		case "TYR": return "Y" 
		case "VAL": return "V" 
		case other:
			import mout
			mout.errorOut(f"Unsupported amino acid residue name: {name}")
			return None

def longname(name):
	match name:
		case "ALA": return "Alanine" 
		case "ARG": return "Arginine" 
		case "ASN":	return "Asparagine" 
		case "ASP": return "Aspartic Acid" 
		case "CYS": return "Cysteine" 
		case "GLN": return "Glutamine" 
		case "GLU": return "Glutamic Acid" 
		case "GLY": return "Glycine" 
		case "HIS": return "Histidine" 
		case "HSD": return "Histidine" 
		case "ILE": return "Isoleucine" 
		case "LEU": return "Leucine" 
		case "LYS": return "Lysine" 
		case "MET": return "Methionine" 
		case "PHE": return "Phenylalanine" 
		case "PRO": return "Proline" 
		case "SER": return "Serine" 
		case "THR": return "Threonine" 
		case "TRP": return "Tryptophan" 
		case "TYR": return "Tyrosine" 
		case "VAL": return "Valine" 
		case other:
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
	def HA(self):
		"""Get Carbon-Alpha Hydrogen Atom"""
		return self.get_atom("HA")

	def mutate(self,newname):
		"""Mutate this amino acid to another"""
		import mout
		import mcol
		mout.headerOut(f"Mutating {mcol.arg}{self.longname}{mcol.clear+mcol.bold} --> {mcol.arg}{longname(newname)}{mcol.clear+mcol.bold} ({mcol.result}{self.letter}{self.number}{alphabet(newname)}{mcol.clear+mcol.bold})")
		match newname:
			case "ALA":
				from .mutate.protein import to_alanine
				to_alanine(self)
			case other:
				import mout
				mout.errorOut("Unsupported mutation!")
	