
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
	
def replace_sidechain(residue,reference):
	"""Remove the amino acid's sidechain and replace with one from the reference residue"""

	import mout
	from ..amino import AminoAcid
	assert isinstance(residue,AminoAcid)
	assert isinstance(reference,AminoAcid)

	# common sidechain atoms
	common = list(set(residue.sidechain_names).intersection(reference.sidechain_names))
	nonH_common = [n for n in common if not n.startswith("H")]
	mout.out(f"{residue} and {reference} share {len(nonH_common)} non-hydrogen side chain atoms")

	if len(nonH_common) > 1:

		# align using these atoms
		align = ['CA'] + nonH_common[:3]

		# align the reference to the target
		res_map = [residue.get_atom(n) for n in align]
		ref_map = [reference.get_atom(n) for n in align]
		reference.align_by_posmap([ref_map,res_map])

		# replace the sidechain
		residue.remove_sidechain(verbosity=0)
		for atom in reference.sidechain:
			residue.add_atom(atom)

		# rename the residue
		residue.name = reference.name

	else:

		mout.errorOut("Unsupported (residues do not share enough sidechain atoms)")
