
def view(atoms,**kwargs):
	from .system import System
	if isinstance(atoms,System):
		atoms = atoms.ase_atoms
	from ase import visualize
	visualize.view(atoms,**kwargs)
	