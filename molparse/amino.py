
from .residue import Residue, res_type

def alphabet(name):
	try:
		return _letter_from_resname[name]
	except KeyError:
		mout.errorOut(f"Unsupported amino acid residue name: {name}")
		return None

def longname(name):
	try:
		return _longname_from_resname[name]
	except KeyError:
		mout.errorOut(f"Unsupported amino acid residue name: {name}")
		return None

_letter_from_resname = { 
	"ALA": "A",
	"ARG": "R",
	"ASN": "N",
	"ASP": "D",
	"CYS": "C",
	"GLN": "Q",
	"GLU": "E",
	"GLY": "G",
	"HIS": "H",
	"HSD": "H",
	"HSE": "H",
	"HSP": "H",
	"ILE": "I",
	"LEU": "L",
	"LYS": "K",
	"MET": "M",
	"PHE": "F",
	"PRO": "P",
	"SER": "S",
	"THR": "T",
	"TRP": "W",
	"TYR": "Y",
	"VAL": "V"
}

_resname_from_letter = {
	"A":"ALA",
	"R":"ARG",
	"N":"ASN",
	"D":"ASP",
	"C":"CYS",
	"Q":"GLN",
	"E":"GLU",
	"G":"GLY",
	"H":"HIS",
	"H":"HSD",
	"H":"HSE",
	"H":"HSP",
	"I":"ILE",
	"L":"LEU",
	"K":"LYS",
	"M":"MET",
	"F":"PHE",
	"P":"PRO",
	"S":"SER",
	"T":"THR",
	"W":"TRP",
	"Y":"TYR",
	"V":"VAL"
}

_longname_from_resname = {
	"ALA":"Alanine",
	"ARG":"Arginine",
	"ASN":"Asparagine",
	"ASP":"Aspartic Acid",
	"CYS":"Cysteine",
	"GLN":"Glutamine",
	"GLU":"Glutamic Acid",
	"GLY":"Glycine",
	"HIS":"Histidine",
	"HSD":"Histidine (Delta)",
	"HSE":"Histidine (Epsilon)",
	"HSP":"Histidine (Protonated)",
	"ILE":"Isoleucine",
	"LEU":"Leucine",
	"LYS":"Lysine",
	"MET":"Methionine",
	"PHE":"Phenylalanine",
	"PRO":"Proline",
	"SER":"Serine",
	"THR":"Threonine",
	"TRP":"Tryptophan",
	"TYR":"Tyrosine",
	"VAL":"Valine"
}

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
		self._sidechain_names = None

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

	def remove_backbone(self,add_link=True,verbosity=2):
		"""Remove backbone atoms"""

		if add_link:
			if verbosity:
				mout.warningOut("Adding HLNK linker")
			self.CA.set_name('HLNK')

		for atom in self.backbone:
			if not atom:
				continue
			self.delete_atom(atom.name,verbosity=verbosity-1)

	def remove_sidechain(self,verbosity=2):
		"""Remove sidechain atoms"""

		for atom in self.sidechain:
			if not atom:
				continue
			self.delete_atom(atom.name,verbosity=verbosity-1)

	@property
	def sidechain_names(self):
		if not self._sidechain_names:
			bb_names = ["N","HN","CA","HA","C","O"]
			self._sidechain_names = [n for n in self.atom_names() if n not in bb_names]

		return self._sidechain_names
	
	@property
	def sidechain(self):
		"""Get sidechain atoms"""

		if not self._sidechain:
			self._sidechain = self.get_atom(self.sidechain_names)

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

		import os
		import mout
		import mcol
		import glob

		mout.headerOut(f"Mutating {mcol.arg}{self.longname}{mcol.clear+mcol.bold} --> {mcol.arg}{longname(newname)}{mcol.clear+mcol.bold} ({mcol.result}{self.letter}{self.number}{alphabet(newname)}{mcol.clear+mcol.bold})")

		# look for reference files
		amp_path = os.path.dirname(__file__)
		files = glob.glob(f'{amp_path}/ref/???.pdb')
		keys = [os.path.basename(f).replace(".pdb","") for f in files]

		if newname ==  "ALA":
			from .mutate.protein import to_alanine
			to_alanine(self)

		elif newname in keys:

			from .io import parse
			from .mutate.protein import replace_sidechain
			reference = parse(f'{amp_path}/ref/{newname}.pdb').residues[0]
			replace_sidechain(self,reference)

		else:
			import mout
			mout.errorOut("Unsupported mutation!")
	