
def interpolate(input1,input2,frames=10,verbosity=2,smooth=False,indices=None):
	import mout
	from ase import aseatoms
	from .system import System

	if verbosity > 0:
		mout.headerOut("Interpolating...")
	if verbosity > 1:
		mout.varOut("input1",str(input1))
		mout.varOut("input2",str(input2))
		mout.varOut("frames",frames)

	if isinstance(input1,aseatoms):
		mout.out("Input1 = ase.Atoms")
		mout.errorOut("Not supported yet",fatal=True)
	elif isinstance(input1,System):
		mout.out("Input1 = amp.system")
	else:
		mout.out(str(type(input1)))
		mout.errorOut("Not supported yet",fatal=True)

	system = input1.copy()

	system_array=[]

	for i in range(frames): 

		system.name = "Interpolation Frame "+str(i)

		if indices is not None:
			for index in indices:
				atom = system.atoms[index]

				if smooth:
					atom.position = smooth_interpolate(input1.atoms[index].np_pos,
													   input2.atoms[index].np_pos,
													   frames,i)
				else:
					atom.position = simple_interpolate(input1.atoms[index].np_pos,
													   input2.atoms[index].np_pos,
													   frames,i)

		else:
			for index,atom in enumerate(system.atoms):

				if smooth:
					atom.position = smooth_interpolate(input1.atoms[index].np_pos,
													   input2.atoms[index].np_pos,
													   frames,i)
				else:
					atom.position = simple_interpolate(input1.atoms[index].np_pos,
													   input2.atoms[index].np_pos,
													   frames,i)
			
		system_array.append(system.copy())

	return system_array

def simple_interpolate(start,end,frames,i):
	return start+i*(end-start)/(frames-1)

def smooth_interpolate(start,end,frames,i):
	import math
	angle = simple_interpolate(0.0,math.pi,frames,i)
	return start + 0.5*(1-math.cos(angle))*(end-start)

def auto_rotate(atoms):

	atoms = atoms.copy()

	positions = atoms.get_positions()

	positions = [
				[ p[0] for p in positions],
				[ p[1] for p in positions],
				[ p[2] for p in positions],
				]

	import numpy as np
	from numpy.linalg import svd
	
	# fit a plane to the atomic positions
	points = np.reshape(positions, (np.shape(positions)[0], -1)) # Collapse trialing dimensions
	assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
	central_point = points.mean(axis=1)
	x = points - central_point[:,np.newaxis]
	M = np.dot(x, x.T) # Could also use np.cov(x) here.
	normal_vector = svd(M)[0][:,-1]

	# rotate the normal onto the Z-axis
	atoms.rotate(normal_vector,'z')

	return atoms
