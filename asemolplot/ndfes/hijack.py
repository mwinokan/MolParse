
import numpy as np
from ase.calculators.calculator import Calculator
from ase import Atoms
import mout
import mcol

class FakeAtoms(Atoms):

	def __init__(self,rcs=None,**kwargs):
		mout.debugHeader("FakeAtoms.__init__()")

		if rcs is not None:

			symbols = len(rcs)*'H'

			positions = []

			for r in rcs:
				positions.append((r,0.0,0.0))

			Atoms.__init__(self,symbols,positions=positions,**kwargs)
		
		else:
			
			Atoms.__init__(self,**kwargs)

	def set_constraints(self,locut=0.5,hicut=3.0,k=20):

		del self.constraints

		from ase.constraints import Hookean

		constraints = []

		for i in range(self.__len__()):
			constraints.append(Hookean(a1=i,a2=(1,0,0,-hicut),k=k))
			constraints.append(Hookean(a1=i,a2=(-1,0,0,locut),k=k))

			# constraints.append(Hookean(a1=i,a2=(-1,0,0,-hicut),k=k))
			# constraints.append(Hookean(a1=i,a2=(1,0,0,locut),k=k))

		self.set_constraint(constraints)

	def positions_from_rcs(self,rcs):
		positions = []
		for r in rcs:
			positions.append((r,0.0,0.0))
		self.positions = positions

	@property
	def rcs(self):
		return [p[0] for p in self.positions]	

	def set_positions(self, newpositions, apply_constraint=False):
		newpositions = np.array([np.array([p[0],0.0,0.0]) for p in newpositions])
		self.set_array('positions', newpositions, shape=(3,))

class FakeSurfaceCalculator(Calculator):

	"""
		
		Calculator to hijack BFGS, NEB and other algorithms.

		A multidimensional free energy surface which is described by several 1D reaction coordinates, is described by a proxy Atoms object.

	"""

	implemented_properties = ['energy','forces']

	# default_parameters = {'delta': 0.01}

	nolabel = True

	def __init__(self, pmf, delta=0.001, suppress_warnings=False, **kwargs):
		mout.debugHeader(f"FakeSurfaceCalculator.__init__(delta={delta})")
		self._pmf = pmf
		self.delta = delta
		self._history = []
		self._warnings = not suppress_warnings
		Calculator.__init__(self, **kwargs)
	
	def pmf(self,arg):
		from ase.units import kcal,mol
		return self._pmf(arg)*(kcal/mol)

	@property
	def history(self):
		return self._history

	@history.setter
	def history(self,arg):
		self._history = arg

	def calculate(self,atoms=None,properties=['energy','forces'],system_changes=['positions','numbers','cell','pbc']):
		Calculator.calculate(self, atoms, properties, system_changes)

		self.energy = 0.0

		# print("calculate:",self.atoms.get_positions())

		xs = np.array([p[0] for p in self.atoms.get_positions()])

		energy = self.pmf([xs])
		mout.debugOut(f"xs: {xs}, pmf(xs): {energy}")

		# if energy > 40:
		# 	mout.warningOut("Unsampled region! Going back one step!")
		# 	xs = self.history[-1]
		# 	energy = self.pmf([xs])
		# 	mout.debugOut(f"xs: {xs}, pmf(xs): {energy}")
		
		self.history.append(xs)

		forces = []

		# partial derivatives
		for i in range(len(self.atoms)):

			h = np.array([0.0 for i in range(len(self.atoms))])
			h[i] = self.delta
			y_h = self.pmf([xs+h])
			y = self.pmf([xs])

			# print(f"xs+h: {xs+h}")
			# print(f"pmf(xs+h): {y_h}")

			if self._warnings and y > 40:
				mout.warningOut("Unsampled region!")
				mout.out(f"{mcol.warning}xs: {xs}, pmf(xs): {y}")

			dx = (y_h - y)/self.delta
			# print(dx)
			forces.append(-dx[0])

		mout.debugOut(f"forces: {forces}")

		forces = np.array([np.array([f,0.0,0.0]) for f in forces])

		# print(forces)

		### convert forces to correct units

		self.results['energy'] = energy # / 100.0 ### now in eV
		self.results['forces'] = forces# / 100.0
