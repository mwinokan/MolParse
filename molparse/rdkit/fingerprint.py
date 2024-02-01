
from .mol import mol_from_smiles
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

def get_similar(mols, query, index=False, threshold=1.0):

	if isinstance(query, str):
		query = mol_from_smiles(query)

	kwargs = dict(minPath=1, maxPath=7, fpSize=2048, bitsPerHash=2, useHs=True, tgtDensity=0.0, minSize=128)	

	fp = FingerprintMols.FingerprintMol(query, **kwargs)
	fps = [FingerprintMols.FingerprintMol(m, **kwargs) for m in mols]

	scores = DataStructs.BulkTanimotoSimilarity(fp, fps)

	if index:
		return [i for i,s in enumerate(scores) if s >= threshold]
	else:
		return [m for m,s in zip(mols, scores) if s >= threshold]
