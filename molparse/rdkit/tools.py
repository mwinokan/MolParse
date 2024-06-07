
from mlog import setup_logger
logger = setup_logger('MolParse')

def sdf_align(sdf_path: str, source_protein_path: str, target_protein_path: str, mol_col='ROMol', name_col='ID', backbone_only=True, heavy_atoms_only=True):

	from ..io import parse
	from numpy import array, double
	from rdkit.Chem import PandasTools, rdMolTransforms, Mol
	from .mol import mol_from_smiles, mol_to_smiles, mol_to_pdb_block

	# parse the reference files
	source = parse(source_protein_path)
	target = parse(target_protein_path)

	# align the proteins
	source_center, target_center, rotation = source.align_to(
							target, 
							protein_only=True, 
							return_transformations=True, 
							backbone_only=backbone_only, 
							heavy_atoms_only=heavy_atoms_only
						)

	# create the transformation matrices

	translation_1 = array([
	    [1, 0, 0, -source_center[0]],
	    [0, 1, 0, -source_center[1]],
	    [0, 0, 1, -source_center[2]],
	    [0, 0, 0, 1],
	], dtype=double)

	translation_2 = array([
	    [1, 0, 0, target_center[0]],
	    [0, 1, 0, target_center[1]],
	    [0, 0, 1, target_center[2]],
	    [0, 0, 0, 1],
	], dtype=double)

	rotation = array([
		[*rotation[0], 0],
		[*rotation[1], 0],
		[*rotation[2], 0],
		[0, 0, 0, 1],
	], dtype=double)

	df = PandasTools.LoadSDF(sdf_path)

	for i,row in df.iterrows():
		if row[name_col] == 'ver_1.2':
			logger.warning('Skipping Fragalysis header molecule')
			continue

		mol = row[mol_col]

		# translate to origin
		rdMolTransforms.TransformConformer(mol.GetConformer(0), translation_1)
		
		# apply rotation
		rdMolTransforms.TransformConformer(mol.GetConformer(0), rotation)

		# translate to target
		rdMolTransforms.TransformConformer(mol.GetConformer(0), translation_2)

	out_path = sdf_path.replace('.sdf', '_aligned.sdf')

	logger.writing(out_path)
	PandasTools.WriteSDF(df, out_path, mol_col, name_col, list(df.columns))
