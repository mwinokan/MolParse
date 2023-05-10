
import numpy as np
from ase.calculators.calculator import Calculator
from ase import Atoms
import mout
import mcol

class NaNEncounteredInPmfError(Exception):
	pass

class FakeAtoms(Atoms):

	"""This class creates an object resembling ase.Atoms that encodes any number (N) of reaction coordinates into the x-positions of ase.Atom objects. Use this class together with FakeSurfaceCalculator to perform optimisations inside an N-dimensional PMF.

	Initialise as:

	FakeAtoms(rcs=[...],name=...)

	- rcs is a list of floats representing a vector in reaction-coordinate space
	- name is an optional string naming this point in RC-space

	"""

	def __init__(self,rcs=None,name=None,**kwargs):
		mout.debugOut(f"FakeAtoms.__init__({rcs=},{name=},...)")

		self.name = name

		if rcs is not None:

			symbols = len(rcs)*'H'

			positions = []

			for r in rcs:
				positions.append((r,0.0,0.0))

			Atoms.__init__(self,symbols,positions=positions,**kwargs)
		
		else:
			
			Atoms.__init__(self,**kwargs)

	def set_constraints(self,locut=0.5,hicut=3.0,k=20,additional=None):

		del self.constraints

		from ase.constraints import Hookean

		constraints = []

		for i in range(self.__len__()):
			constraints.append(Hookean(a1=i,a2=(1,0,0,-hicut),k=k))
			constraints.append(Hookean(a1=i,a2=(-1,0,0,locut),k=k))

			# constraints.append(Hookean(a1=i,a2=(-1,0,0,-hicut),k=k))
			# constraints.append(Hookean(a1=i,a2=(1,0,0,locut),k=k))

		if additional is not None:
			constraints.append(additional)

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
	
		Calculator to hijack ASE optimisation algorithms.

		Initialise as:

		FakeSurfaceCalculator(pmf=function())

		pmf is a N-dimensional free energy surface. Should be callable with the same number of arguments as reaction coordinates described in the FakeAtoms object.

		Optimisation example:

		atoms = FakeAtoms(rcs=[1.0,1.0])

		def pmf(x,y):
			...

		calc = FakeSurfaceCalculator(pmf=pmf)

		atoms.set_calculator(calc)

		from ase.optimize import BFGS

		dyn = BFGS(atoms=atoms)

		dyn.run()

	"""

	implemented_properties = ['energy','forces']

	# default_parameters = {'delta': 0.01}

	nolabel = True

	def __init__(self, pmf, delta=0.001, suppress_warnings=False, restoring_force = 10, **kwargs):
		mout.debugOut(f"FakeSurfaceCalculator.__init__(delta={delta})")
		self._pmf = pmf
		self.delta = delta
		self._history = []
		self._last_ok = None
		self._restoring_force = restoring_force
		self._warnings = not suppress_warnings
		Calculator.__init__(self, **kwargs)
	
	def pmf(self,arg):
		from ase.units import kcal,mol
		e = self._pmf(arg)

		if np.isnan(e):
			mout.errorOut(f"NaN detected in PMF! {arg=}")
			raise NaNEncounteredInPmfError

		return self._pmf(arg)

	@property
	def history(self):
		return self._history

	@history.setter
	def history(self,arg):
		self._history = arg

	def calculate(self,atoms=None,properties=['energy','forces'],system_changes=['positions','numbers','cell','pbc']):
		Calculator.calculate(self, atoms, properties, system_changes)

		self.energy = 0.0

		xs = np.array([p[0] for p in self.atoms.get_positions()])

		energy = self.pmf([xs])
		mout.debugOut(f"xs: {xs}, pmf(xs): {energy}")
		
		self.history.append(xs)

		forces = []

		y = self.pmf([xs])

		self._last_ok = [x for x in xs]
	
		# partial derivatives
		for i in range(len(self.atoms)):

			h = np.array([0.0 for i in range(len(self.atoms))])
			h[i] = self.delta
			y_h = self.pmf([xs+h])

			dx = (y_h - y)/self.delta
			forces.append(-dx[0])

		mout.debugOut(f"forces: {forces}")

		forces = np.array([np.array([f,0.0,0.0]) for f in forces])

		### convert forces to correct units

		self.results['energy'] = energy
		self.results['forces'] = forces
