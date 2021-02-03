
import mout
import mcol
import ase

from .system import System

def interpolate(input1,input2,filename,frames=10,verbosity=2):

	if verbosity > 0:
		mout.headerOut("Interpolating...")
	if verbosity > 1:
		mout.varOut("input1",str(input1))
		mout.varOut("input2",str(input2))
		mout.varOut("filename",filename)
		mout.varOut("frames",frames)

	if isinstance(input1,ase.Atoms):
		mout.out("Input1 = ase.atoms")
		mout.errorOut("Not supported yet",fatal=True)
	elif isinstance(input1,System):
		mout.out("Input1 = amp.system")
	else:
		mout.out(str(type(input1)))
		mout.errorOut("Not supported yet",fatal=True)

	start=0.0
	end=1.0

	system = input1.copy()

	system_array=[]

	for i in range(frames): 

		system.name = "Interpolation Frame "+str(i)

		for index,atom in enumerate(system.atoms):


			atom.position = simple_interpolate(input1.atoms[index].np_pos,
											   input2.atoms[index].np_pos,
											   frames,i)
			
			if index == 20:
				print(simple_interpolate(input1.atoms[index].np_pos,
											   input2.atoms[index].np_pos,
											   frames,i))

		system_array.append(system.copy())

	print(system_array)

	return system_array

		# print(simple_interpolate(start,end,frames,i))

def simple_interpolate(start,end,frames,i):
	return start+i*(end-start)/(frames-1)

# def vector_interpolate(start,end,frames,i):
# 	return start+i*(end-start)/(frames-1)
