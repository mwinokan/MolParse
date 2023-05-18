
def get_link_atom(site,reference_atom,l=1.0):
	"""Get a hydrogen link atom positioned at l Å away from site towards reference_atom"""
	import numpy as np
	from ..atom import Atom

	direction = reference_atom - site
	new_position = site + l*(direction/np.linalg.norm(direction))

	atom = Atom("H",position=list(new_position),chain=site.chain,residue=site.residue,res_number=site.res_number,res_index=site.res_index,occupancy=1.0,temp_factor=0.0)
	# atom = Atom("H",position=list(new_position))

	return atom

