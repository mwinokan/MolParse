
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

	def __init__(self,name: str,index: int=None,number: int=None, chain: str=None,atoms=None):
		assert res_type(name) == "PRO"
		super(AminoAcid, self).__init__(name,index,number,chain,atoms)

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
			self._backbone = self.get_atom(names,verbosity=0)

		return self._backbone

	def guess_hydrogen_names(self,verbosity=1):
		import mout
		import numpy as np
		for hydrogen in [a for a in self.atoms if a.symbol == 'H']:
			others = [a for a in self.atoms if a.symbol != 'H' and a.index != hydrogen.index]
			other_index = np.argmin([np.linalg.norm(a - hydrogen) for a in others])
			other = others[other_index]
			newname = f'H{other.name}'
			if verbosity:
				mout.warning(f'Renaming: {hydrogen} {hydrogen.number} --> {newname}')
			hydrogen.name = newname

	def remove_backbone(self,add_link=True,copy=False,verbosity=2):
		"""Remove backbone atoms"""

		if copy:
			res = self.copy()
		else:
			res = self

		# res.summary()

		import mout
		if add_link:
			if verbosity:
				import mout
				mout.warningOut("Adding HLNK linker")
			import numpy as np
			vec_CB_CA = res.CA - res.CB
			unit_CB_CA = vec_CB_CA/np.linalg.norm(vec_CB_CA)
			target_length = COVALENT_RADII['H'] + COVALENT_RADII['C']
			res.CA.position = res.CB.np_pos + target_length * unit_CB_CA
			res.CA.set_name('HLNK',verbosity=verbosity-1)

		for atom in res.backbone:
			if atom is None:
				continue
			res.delete_atom(atom.name,verbosity=verbosity-1)

		return res

	def remove_sidechain(self,verbosity=2):
		"""Remove sidechain atoms"""

		for atom in self.sidechain:
			if atom is None:
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
			self._sidechain = self.get_atom(self.sidechain_names,verbosity=0)

		return self._sidechain

	@property
	def CA(self):
		"""Get Carbon-Alpha Atom"""
		return self.get_atom("CA",verbosity=0)

	@property
	def CB(self):
		"""Get Carbon-Alpha Atom"""
		return self.get_atom("CB",verbosity=0)

	@property
	def HA(self):
		"""Get Carbon-Alpha Hydrogen Atom"""
		return self.get_atom("HA",verbosity=0)

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
	
# from ASE
COVALENT_RADII = {
	'X': 0.2, 'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69, 'Bk': 0.2, 'Cf': 0.2, 'Es': 0.2, 'Fm': 0.2, 'Md': 0.2, 'No': 0.2, 'Lr': 0.2, 'Rf': 0.2, 'Db': 0.2, 'Sg': 0.2, 'Bh': 0.2, 'Hs': 0.2, 'Mt': 0.2, 'Ds': 0.2, 'Rg': 0.2, 'Cn': 0.2, 'Nh': 0.2, 'Fl': 0.2, 'Mc': 0.2, 'Lv': 0.2, 'Ts': 0.2, 'Og': 0.2
}