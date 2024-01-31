
import mout
from rdkit import Chem

# class Molecule(Chem.rdchem.Mol):
# 	def __init__(self,*args,**kwargs):
# 		super(ClassName, self).__init__(*args,**kwargs)

def mol_from_pdb_block(pdb_block):
	return Chem.rdmolfiles.MolFromPDBBlock(pdb_block)

def mol_to_smiles(mol):
	return Chem.MolToSmiles(mol)

def mol_from_smiles(smiles):
	return Chem.MolFromSmiles(smiles)

def mol_to_pdb_block(mol):
	return Chem.MolToPDBBlock(mol)

def mol_to_AtomGroup(mol):
	from ..group import AtomGroup
	group = AtomGroup.from_pdb_block(mol_to_pdb_block(mol))
	if hasattr(mol, '_Name'):
		group.name = mol._Name
	return group

def mol_from_AtomGroup(group):
	return group.mol

def protonate(mol,embed=True,align=True,verbosity=1):
	mol_prot = Chem.AddHs(mol)
	if embed and align:
		if verbosity:
			mout.warning('May lose exact pose',code="mp.rdkit.protonate.1")
		mol = Chem.AllChem.ConstrainedEmbed(mol_prot, mol)
	elif embed:
		if verbosity:
			mout.warning('Disregarding original coordinates',code="mp.rdkit.protonate.2")
		ps = Chem.AllChem.ETKDGv3()
		Chem.AllChem.EmbedMolecule(mol_prot,ps)
	return mol_prot
