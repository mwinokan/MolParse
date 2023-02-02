
def methylate(residue,site,position=None):

	import numpy as np
	from ..atom import Atom

	if position is None:
		position = site.np_pos + np.array([1.538,0,0])

	l0_CB_HB = 1.111 # Angstrom
	a0_HB_CB_CS = 110.10 # degrees

	# place CB
	CB = Atom("CB",position=position,residue=residue.name,chain=residue.chain,res_number=residue.number)
	residue.add_atom(CB)

	# distances
	dist_Z = np.sin(np.pi/180*(a0_HB_CB_CS - 90))*l0_CB_HB
	r_XY = np.cos(np.pi/180*(a0_HB_CB_CS - 90))*l0_CB_HB

	# vectors
	vec_CA_CB = CB - site
	vec_X = np.cross(vec_CA_CB,[1,0,0])
	vec_X /= np.linalg.norm(vec_X)
	vec_Y = np.cross(vec_CA_CB,vec_X)
	vec_Y /= np.linalg.norm(vec_Y)

	# place HB1	
	dist_X = r_XY*np.sin(0.0)
	dist_Y = r_XY*np.cos(0.0)
	pos_HB1 = CB + vec_CA_CB*dist_Z + vec_X*dist_X + dist_Y*vec_Y
	HB1 = Atom("HB1",position=pos_HB1,residue=residue.name,chain=residue.chain,res_number=residue.number)
	residue.add_atom(HB1)

	# place HB2
	dist_X = r_XY*np.sin(2/3*np.pi)
	dist_Y = r_XY*np.cos(2/3*np.pi)
	pos_HB2 = CB + vec_CA_CB*dist_Z + vec_X*dist_X + dist_Y*vec_Y
	HB2 = Atom("HB2",position=pos_HB2,residue=residue.name,chain=residue.chain,res_number=residue.number)
	residue.add_atom(HB2)

	# place HB3
	dist_X = r_XY*np.sin(4/3*np.pi)
	dist_Y = r_XY*np.cos(4/3*np.pi)
	pos_HB3 = CB + vec_CA_CB*dist_Z + vec_X*dist_X + dist_Y*vec_Y
	HB3 = Atom("HB3",position=pos_HB3,residue=residue.name,chain=residue.chain,res_number=residue.number)
	residue.add_atom(HB3)
