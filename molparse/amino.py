
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
			self._backbone = [a for a in self.get_atom(BB_NAMES,verbosity=0) if a is not None]

		return self._backbone

	def guess_hydrogen_names(self,verbosity=1):
		import mout
		import numpy as np
		for hydrogen in [a for a in self.atoms if a.symbol == 'H']:
			others = [a for a in self.atoms if a.symbol != 'H' and a.index != hydrogen.index]

			if len(others) == 0:
				mout.var('this_residue',self.name_number_str)
				self.plot3d()
				raise Exception('No other atoms in residue when guessing hydrogen names')
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
			target_length = 1.07 # from ase covalent radii
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
			self._sidechain_names = [n for n in self.atom_names() if n not in BB_NAMES]

		return self._sidechain_names
	
	@property
	def sidechain(self):
		"""Get sidechain atoms"""

		if not self._sidechain:
			self._sidechain = [a for a in self.get_atom(self.sidechain_names,verbosity=0) if a is not None]

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

	@property
	def interaction_sites(self):
		
		import numpy as np
		from collections import namedtuple

		from .sites import Site
		
		sites = []

		for interaction in BB_INTERACTION_SITES + INTERACTION_SITES[self.name]:

			atoms = interaction['atoms']

			if len(atoms) == 1:
				atoms = [self.get_atom(atoms[0])]
				position = atoms[0].np_pos
			else:
				atoms = self.get_atom(atoms)
				position = np.mean([a.np_pos for a in atoms],axis=0).round(3)

			# check if an equivalent site exists

			matches = [s for s in sites if np.linalg.norm(s.position - position) < 0.1]
			if len(matches) == 1:
				matches[0].types.append(interaction['type'])
				
			elif len(matches) == 0:
				sidechain = atoms[0].name in self.sidechain_names
				site = Site([interaction['type']],atoms,position,sidechain,self.name,self.number,self.chain)
				sites.append(site)

			else:
				raise Exception("Should never have multiple sites with the same position")

		return sites

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

	@property
	def features(self):

		import numpy as np
		from .rdkit import Feature

		features = []

		for f_dict in BB_FEATURES + FEATURES[self.name]:

			atoms = f_dict['atoms']

			if len(atoms) == 1:
				atoms = [self.get_atom(atoms[0])]
				if atoms[0] is None:
					self.plot3d()
					mout.error(f'No atom {f_dict["atoms"][0]}',code='AminoAcid.features.0',fatal=True)
				position = atoms[0].np_pos
			else:
				atoms = self.get_atom(atoms)
				position = np.mean([a.np_pos for a in atoms],axis=0).round(3)

			feature = Feature(
				family=f_dict['family'], 
				atoms=atoms,
				position=position,
				sidechain=None,
				res_name=self.name,
				res_number=self.number,
				res_chain=self.chain,
			)

			# print(f_dict)
			features.append(feature)

		return features

"""

	Collate interactions per site

	Which interactions are shared among all amino acid backbones?

	Auto-checking rules:

		+ Hydrogen donor infers water donor
		+ Hydrogen acceptor infers water acceptor
		+ Most carbons should be able to form hydrophobic interactions
		+ Can hydrogen donors also be acceptors?


"""

RES_NAMES = [
	"ALA",
	"ARG",
	"ASN",
	"ASP",
	"CYS",
	"GLN",
	"GLU",
	"GLY",
	"HIS",
	"HSD",
	"HSE",
	"HSP",
	"ILE",
	"LEU",
	"LYS",
	"MET",
	"PHE",
	"PRO",
	"SER",
	"THR",
	"TRP",
	"TYR",
	"VAL",
]

BB_NAMES = ["N", "HN", "CA", "HA", "C", "O"]

BB_INTERACTION_SITES = [
	{"type": "water_donor", "atoms": ["N"], "source": "collated from ZV PLIP"},
	{"type": "hydrogen_donor", "atoms": ["N"], "source": "collated from ZV PLIP"},
	{"type": "water_acceptor", "atoms": ["O"], "source": "collated from ZV PLIP"},
	{"type": "hydrogen_acceptor", "atoms": ["O"], "source": "collated from ZV PLIP"},
]

BB_FEATURES = [
	dict(family="Donor", atoms=["N"], source="collated from ZV PLIP"),
	dict(family="Acceptor", atoms=["O"], source="collated from ZV PLIP"),

	dict(family="Donor", atoms=["O"], source="unsure"),
	dict(family="Acceptor", atoms=["N"], source="unsure"),
]

FEATURES = {
	"ALA": [
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
	],

	"ASN": [
		dict(family="Acceptor", atoms=["ND2"], source="MPro PLIP"),
		dict(family="Acceptor", atoms=["OD1"], source="unsure"),
		dict(family="Donor", atoms=["OD1"], source="MPro PLIP"),
		dict(family="Hydrophobe", atoms=["CB"], source="MPro PLIP"),
	],

	"ASP": [
		dict(family="Acceptor", atoms=["OD1"], source="unsure"),
		dict(family="Donor", atoms=["OD1"], source="unsure"),
		dict(family="Acceptor", atoms=["OD2"], source="unsure"),
		dict(family="NegIonizable", atoms=["OD1","OD2"], source="unsure"),
	],

	"ASP": [
		dict(family="Acceptor", atoms=["OD1"], source="ZV PLIP"),
		dict(family="Donor", atoms=["OD1"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["OD2"], source="ZV PLIP"),
		dict(family="NegIonizable", atoms=["OD1","OD2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
	],

	"ARG": [
		dict(family="Donor", atoms=["NH1"], source="unsure"),
		dict(family="Donor", atoms=["NH2"], source="unsure"),
		dict(family="Donor", atoms=["NE"], source="unsure"),
		dict(family="PosIonizable", atoms=["NH1"], source="unsure"),
		dict(family="PosIonizable", atoms=["NH2"], source="unsure"),
	],

	"CYS": [
		dict(family="Hydrophobe", atoms=["CB"], source="unsure"),
		dict(family="PosIonizable", atoms=["SG"], source="unsure"),
		dict(family="Donor", atoms=["SG"], source="unsure"),
	],
	
	"GLN": [
		dict(family="Hydrophobe", atoms=["CG"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["NE2"], source="inferred (ASN)"),
		dict(family="Acceptor", atoms=["OE1"], source="inferred (ASN)"),
		dict(family="Donor", atoms=["OE1"], source="inferred (ASN)"),
		dict(family="Hydrophobe", atoms=["CB"], source="inferred (ASN)"),
	],
	
	"GLU": [
		dict(family="Acceptor", atoms=["OE1"], source="ZV PLIP"),
		dict(family="Donor", atoms=["OE1"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["OE2"], source="ZV PLIP"),
		dict(family="NegIonizable", atoms=["OE1","OE2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CG"], source="ZV PLIP"),
	],

	"GLY": [
	],

	"HIS": [
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Aromatic", atoms=["CG", "ND1", "CE1", "NE2", "CD2"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["ND1"], source="unsure"),
		dict(family="Acceptor", atoms=["NE2"], source="unsure"),
		dict(family="Donor", atoms=["ND1"], source="unsure"),
		dict(family="Donor", atoms=["NE2"], source="unsure"),
	],

	"HSE": [
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Aromatic", atoms=["CG", "ND1", "CE1", "NE2", "CD2"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["ND1"], source="unsure"),
		dict(family="Acceptor", atoms=["NE2"], source="unsure"),
		dict(family="Donor", atoms=["NE2"], source="unsure"),
	],

	"HSD": [
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Aromatic", atoms=["CG", "ND1", "CE1", "NE2", "CD2"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["ND1"], source="unsure"),
		dict(family="Acceptor", atoms=["NE2"], source="unsure"),
		dict(family="Donor", atoms=["ND1"], source="unsure"),
	],

	"HSP": [
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Aromatic", atoms=["CG", "ND1", "CE1", "NE2", "CD2"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["ND1"], source="unsure"),
		dict(family="Acceptor", atoms=["NE2"], source="unsure"),
		dict(family="Donor", atoms=["ND1"], source="unsure"),
		dict(family="Donor", atoms=["NE2"], source="unsure"),
	],

	"ILE": [
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CD1"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CG1"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CG2"], source="ZV PLIP"),
	],

	"LEU": [
		dict(family="Hydrophobe", atoms=["CB"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CD1"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CD2"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CG"], source="unsure"),
	],

	"LYS": [
		dict(family="Hydrophobe", atoms=["CB"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CG"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CD"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CE"], source="unsure"),
		dict(family="PosIonizable", atoms=["NZ"], source="unsure"),
		dict(family="Donor", atoms=["NZ"], source="unsure"),
	],

	"LYS": [
		dict(family="Hydrophobe", atoms=["CB"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CG"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CD"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CE"], source="unsure"),
		dict(family="PosIonizable", atoms=["NZ"], source="unsure"),
		dict(family="Donor", atoms=["NZ"], source="unsure"),
	],

	"MET": [
		dict(family="Hydrophobe", atoms=["CE"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CG"], source="unsure"),
		dict(family="PosIonizable", atoms=["SD"], source="unsure"),
	],

	"PHE": [
		dict(family="Hydrophobe", atoms=["CB"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CD1"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CD2"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CE1"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CE2"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CZ"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CG","CD1","CD2","CE1","CE2","CZ"], source="unsure"),
		dict(family="Aromatic", atoms=["CG","CD1","CD2","CE1","CE2","CZ"], source="unsure"),
	],

	"PRO": [
		dict(family="Hydrophobe", atoms=["CG"], source="ZV PLIP"),
		dict(family="Aromatic", atoms=["N","CA","CB","CG","CD"], source="unsure"),
	],

	"SER": [
		dict(family="Donor", atoms=["OG"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["OG"], source="ZV PLIP"),
	],

	"THR": [
		dict(family="Donor", atoms=["OG1"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["OG1"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CG2"], source="ZV PLIP"),
	],

	"TRP": [
		dict(family="Hydrophobe", atoms=["CZ2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CH2"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CZ3"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CE3"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Aromatic", atoms=["CG","CD1","NE1","CE2","CD2"], source="unsure"),
		dict(family="Aromatic", atoms=["CE2","CD2","CZ2","CZ3","CE3","CH2"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CE2","CD2","CZ2","CZ3","CE3","CH2"], source="unsure"),
	],
	
	"TYR": [
		dict(family="Aromatic", atoms=["CG", "CD1", "CE1", "CZ", "CE2", "CD2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CG", "CD1", "CE1", "CZ", "CE2", "CD2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CD1"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CD2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CE1"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CE2"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CG"], source="ZV PLIP"),
		dict(family="Hydrophobe", atoms=["CB"], source="ZV PLIP"),
		dict(family="Donor", atoms=["OH"], source="ZV PLIP"),
		dict(family="Acceptor", atoms=["OH"], source="ZV PLIP"),
	],

	"VAL": [
		dict(family="Hydrophobe", atoms=["CB"], source="unsure"),
		dict(family="Hydrophobe", atoms=["CG1"], source="ZV PLIP (inferred)"),
		dict(family="Hydrophobe", atoms=["CG2"], source="ZV PLIP"),
	],
}

INTERACTION_SITES = {

	"ALA": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP"},
	], 

	"ASN": [
		# {"type": "hydrophobic", "atoms": ["CB"], "source": "MPro PLIP"},
		# {"type": "hydrophobic", "atoms": ["CG"], "source": "MPro PLIP"},
		# {"type": "water_acceptor", "atoms": ["ND2"], "source": "MPro PLIP"},
		# {"type": "hydrogen_donor", "atoms": ["ND2"], "source": "unsure"},
	],
	
	"ASP": [
		{"type": "hydrogen_acceptor", "atoms": ["OD2"], "source": "ZV PLIP"},
		{"type": "water_acceptor", "atoms": ["OD2"], "source": "ZV PLIP (inferred)"},
		{"type": "hydrogen_acceptor", "atoms": ["OD1"], "source": "ZV PLIP"},
		{"type": "hydrogen_acceptor", "atoms": ["CG"], "source": "ZV PLIP"},
		{"type": "water_acceptor", "atoms": ["CG"], "source": "ZV PLIP (inferred)"},
		{"type": "salt_bridge", "atoms": ["OD1", "OD2"], "source": "ZV PLIP"},
		{"type": "water_acceptor", "atoms": ["OD1"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP"},
		{"type": "hydrogen_donor", "atoms": ["OD1"], "source": "ZV PLIP"},
		{"type": "water_donor", "atoms": ["OD1"], "source": "ZV PLIP (inferred)"},
	], 

	"ARG": [
		{"type": "hydrogen_donor", "atoms": ["NH1"], "source": "unsure"},
		{"type": "hydrogen_donor", "atoms": ["NH2"], "source": "unsure"},
		{"type": "water_donor", "atoms": ["NH1"], "source": "unsure"},
		{"type": "water_donor", "atoms": ["NH2"], "source": "unsure"},
		{"type": "salt_bridge", "atoms": ["NH1"], "source": "unsure"},
		{"type": "salt_bridge", "atoms": ["NH2"], "source": "unsure"},
		{"type": "pi_cation", "atoms": ["NH1"], "source": "unsure"},
		{"type": "pi_cation", "atoms": ["NH2"], "source": "unsure"},
	],

	"CYS": [
		{"type": "pi_cation", "atoms": ["SG"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CB"], "source": "unsure"},
	],

	"GLN": [
		# {"type": "hydrophobic", "atoms": ["CB"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CG"], "source": "ZV PLIP"},
		# {"type": "hydrogen_donor", "atoms": ["NE2"], "source": "unsure"},
		# {"type": "hydrogen_acceptor", "atoms": ["OE1"], "source": "unsure"},
		# {"type": "halogen_acceptor", "atoms": ["O2"], "source": "MPro PLIP"},
	],

	"GLU": [
		{"type": "hydrogen_donor", "atoms": ["OE1"], "source": "ZV PLIP"},
		{"type": "water_donor", "atoms": ["OE1"], "source": "ZV PLIP (inferred)"},
		# {"type": "hydrogen_donor", "atoms": ["CA"], "source": "MPro PLIP"},
		{"type": "hydrogen_acceptor", "atoms": ["OE1"], "source": "unsure"},
		{"type": "water_acceptor", "atoms": ["OE1"], "source": "unsure"},
		# {"type": "hydrogen_donor", "atoms": ["N"], "source": "MPro PLIP"},
		# {"type": "hydrophobic", "atoms": ["CB"], "source": "MPro PLIP"},
	], 

	"GLY": [
		# {"type": "water_donor", "atoms": ["O"]},
	], 

	"HIS": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP"},
		{"type": "pi_stacking", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP"},
		{"type": "pi_cation", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
	], 

	"HSE": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP (inferred)"},
		{"type": "pi_stacking", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
		{"type": "pi_cation", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
	], 

	"HSD": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP (inferred)"},
		{"type": "pi_stacking", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
		{"type": "pi_cation", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
	], 

	"HSP": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP (inferred)"},
		{"type": "pi_stacking", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
		{"type": "pi_cation", "atoms": ["CG", "ND1", "CE1", "NE2", "CD2"], "source": "ZV PLIP (inferred)"},
	], 

	"ILE": [
		{"type": "hydrophobic", "atoms": ["CD1"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP"},
	], 

	"LEU": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CG"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CD1"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CD2"], "source": "unsure"},
	],

	"LYS": [
		{"type": "hydrogen_donor", "atoms": ["NZ"], "source": "unsure"},
		{"type": "water_donor", "atoms": ["NZ"], "source": "unsure"},
		{"type": "salt_bridge", "atoms": ["NZ"], "source": "unsure"},
		{"type": "pi_cation", "atoms": ["NZ"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CB"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CG"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CD"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CE"], "source": "unsure"},
	],
	
	"MET": [
		{"type": "pi_cation", "atoms": ["SD"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CE"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CG"], "source": "unsure"},
		# {"type": "hydrophobic", "atoms": ["CB"], "source": "MPro PLIP"},
	],

	"PHE": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CG"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CD1"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CD2"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CE1"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CE2"], "source": "unsure"},
		{"type": "hydrophobic", "atoms": ["CZ"], "source": "unsure"},
		{"type": "pi_stacking", "atoms": ["CG","CD1","CD2","CE1","CE2","CZ"], "source": "unsure"},
		{"type": "pi_cation", "atoms": ["CG","CD1","CD2","CE1","CE2","CZ"], "source": "unsure"},
	],

	"PRO": [
		{"type": "hydrophobic", "atoms": ["CG"], "source": "ZV PLIP"},
	], 

	"SER": [
		{"type": "hydrogen_acceptor", "atoms": ["OG"], "source": "ZV PLIP"},
		{"type": "water_acceptor", "atoms": ["OG"], "source": "ZV PLIP (inferred)"},
		{"type": "hydrogen_donor", "atoms": ["OG"], "source": "ZV PLIP"},
		{"type": "water_donor", "atoms": ["OG"], "source": "ZV PLIP (inferred)"},
		# {"type": "hydrogen_donor", "atoms": ["CB"], "source": "MPro PLIP"},
	], 

	"THR": [
		{"type": "hydrophobic", "atoms": ["CG2"], "source": "ZV PLIP"},
		{"type": "hydrogen_donor", "atoms": ["OG1"], "source": "ZV PLIP"},
		{"type": "water_donor", "atoms": ["OG1"], "source": "ZV PLIP (inferred)"},
		{"type": "hydrogen_acceptor", "atoms": ["OG1"], "source": "ZV PLIP"},
		{"type": "water_acceptor", "atoms": ["OG1"], "source": "ZV PLIP"},
		# {"type": "water_acceptor", "atoms": ["CB"]},
	], 

	"TRP": [
		{"type": "hydrophobic", "atoms": ["CZ2"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP"},
		{"type": "pi_stacking", "atoms": ["CG", "CD1", "CD2", "NE1", "CE2"], "source": "ZV PLIP"},
		{"type": "pi_cation", "atoms": ["CG", "CD1", "CD2", "NE1", "CE2"], "source": "ZV PLIP (inferred)"},
	], 
	
	"TYR": [
		{"type": "hydrophobic", "atoms": ["CB"], "source": "ZV PLIP"},
		{"type": "pi_stacking", "atoms": ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CE1"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CD1"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CE2"], "source": "ZV PLIP"},
		{"type": "hydrophobic", "atoms": ["CG"], "source": "ZV PLIP"},
		{"type": "pi_cation", "atoms": ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"], "source": "ZV PLIP"},
	], 
	
	"VAL": [
		{"type": "hydrophobic", "atoms": ["CG2"], "source": "ZV PLIP"},
	], 
}
