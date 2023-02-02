
def to_alanine(residue):
	"""Remove the amino acid's sidechain and replace with a methyl group"""
	from ..amino import AminoAcid
	assert isinstance(residue,AminoAcid)

	import mout
	from .functional import methylate
	
	# get the original amino-acids CB position
	pos_CB = residue.CB.np_pos

	# delete sidechain
	names = [a.name for a in residue.sidechain]
	residue.delete_atom(names,verbosity=1)

	# add the methyl group
	methylate(residue,residue.CA,pos_CB)

	# rename the residue
	residue.name = "ALA"
	