
from .hijack import FakeAtoms, FakeSurfaceCalculator
import mout
import mcol

from ase.optimize import BFGS
from ase.units import kcal,mol
from ase.calculators.singlepoint import SinglePointCalculator

class NDFES(object):

	def __init__(self, outdir, n_coords):
		mout.debugHeader("NDFES.__init__()")
		self._pmf = None
		self._meta = None
		self._calc = None
		self.outdir = outdir
		self.checkpoint = f'{outdir}/checkpoint'
		self.metafile = f'{outdir}/metafile'
		self.n_coords = n_coords
		super(NDFES, self).__init__()
	
	def run(self,method="vfep",binlength=0.2,temperature=300,script="eureka.sh"):
		mout.debugHeader("NDFES.run()")

		mout.varOut("metafile",self.metafile,valCol=mcol.file)
		mout.varOut("checkpoint",self.checkpoint,valCol=mcol.file)
		mout.varOut("method",f'{method}',valCol=mcol.arg)
		mout.varOut("temperature",f'{temperature}',unit="K",valCol=mcol.arg)
		mout.varOut("script",f'{script}',valCol=mcol.file)

		import os
		amp_path = os.environ['MWAMPPATH']
		script = f'{amp_path}/molparse/ndfes/{script}'

		length_str = ""
		for i in range(self.n_coords):
			length_str += f"-w {binlength} "

		mout.headerOut("Running ndfes...")
		os.system(f'{script} --{method} -c {self.checkpoint} {length_str} {self.metafile} -t {temperature} > {self.outdir}/ndfes.log')

	@property
	def pmf(self):
		if self._pmf is None:
			self.generate_pmf_object()
		return self._pmf

	@property
	def meta(self):
		if self._meta is None:
			self.generate_pmf_object()
		return self._meta

	@property
	def calc(self):
		return self._calc
	
	@calc.setter
	def calc(self,arg):
		self._calc = arg

	def energy_along_path(self,rcs):
		energies = []
		for xs in rcs:
			e = self.pmf([xs])
			# if e > 1000:
			# 	energies.append(None)
			# else:
			# 	energies.append(e*(kcal/mol))
			energies.append(e*(kcal/mol))
		return energies

	def generate_pmf_object(self):
		mout.debugHeader("NDFES.generate_pmf_object()")
		try:
			import ndfes
		except:
			mout.errorOut("Missing ndfes library. Try sourcing $MWSHPATH/load_ndfes.sh",fatal=True)

		meta = ndfes.Metafile(self.metafile)
		pmf = ndfes.interpolator([self.checkpoint],meta)

		self._meta = meta
		self._pmf = pmf

	def plot_2d_fes(self,dims):
		mout.debugHeader("NDFES.plot_2d_fes()")
		try:
			import ndfes
		except:
			mout.errorOut("Missing ndfes library. Try sourcing $MWSHPATH/load_ndfes.sh",fatal=True)

		mout.varOut("checkpoint",f'{self.checkpoint}',valCol=mcol.file)

		ndgrid = ndfes.NDUniformRegularGrid()
		for d in dims:
			ndgrid.add_dim_using_delta(d[0],d[1],d[2])

		pts = ndgrid.get_list_of_tuples()

		vals = self.pmf(pts)

		# Print the PMF to a file
		fh = open(f'{self.outdir}/2dfes.dat',"w")
		for pt,val in zip(pts,vals):
		    # if pt[0] == -1.0: fh.write("\n")

			# collapse onto 2d

			if val > 1000:
				continue

			fh.write("%15.5f %15.5f %20.10e\n"%(pt[0]-pt[1],pt[2]-pt[3],val))
		
		fh.close()

	def plot_window_pmf(self,i1,i2):
		import matplotlib.pyplot as plt
		import numpy as np
		mout.debugHeader("NDFES.plot_window_pmf()")

		try:
			import ndfes
		except:
			mout.errorOut("Missing ndfes library. Try: "+mcol.file+"source $MWSHPATH/load_ndfes.sh",fatal=True)

		cs = self.meta.get_centers()

		# reactant free energy
		self.reactant = cs[i1]
		mout.varOut(f"Reactant Free Energy (window {i1})",self.pmf([self.reactant])[0],unit="kcal/mol")

		# product free energy
		self.product = cs[i2]
		mout.varOut(f"Product Free Energy (window {i2})",self.pmf([self.product])[0],unit="kcal/mol")

		cs = cs[0:5]

		# reshape into 1d

		vals = self.pmf(cs)

		x = []
		y = []

		mout.headerOut("Window Centres & Energies:")

		for c,v in zip(cs,vals):
			x.append(c[0]-c[1]+c[2]-c[3])
			print(c,v)
			if v < 1000:
				y.append(v*(kcal/mol))
			else:
				y.append(None)

		fig,ax = plt.subplots()

		plt.scatter(x,y)

		plt.savefig(f"{self.outdir}/window_pmf.png")

		plt.close()

	def find_minima(self,start,fmax=0.01,optimizer=BFGS,suppress_warnings=False):
		mout.debugHeader("NDFES.find_minima()")

		# mout.varOut("start",list(start))

		atoms = FakeAtoms(start)

		# print("1",atoms.get_positions())

		self.calc = FakeSurfaceCalculator(self.pmf,suppress_warnings=suppress_warnings)

		atoms.calc = self.calc
		
		# print("2",atoms.get_positions())

		# from ase.optimize import BFGS

		dyn = optimizer(atoms)
		
		# print("3",atoms.get_positions())

		dyn.run(fmax=fmax)

		mout.varOut("Post-BFGS RC's",atoms.rcs)

	def constraint_pair(self,index,locut=0.5,hicut=3.0,k=20):
		from ase.constraints import Hookean
		print(index,locut,hicut)
		if locut > hicut:
			c1 = Hookean(a1=index,a2=(1,0,0,-locut),k=k)
			c2 = Hookean(a1=index,a2=(-1,0,0,hicut),k=k)
		else:
			c1 = Hookean(a1=index,a2=(1,0,0,-hicut),k=k)
			c2 = Hookean(a1=index,a2=(-1,0,0,locut),k=k)
		return [c1,c2]

	def find_path(self,start,final,n_images=10,fmax=0.01,traj=None,optimizer=BFGS,suppress_warnings=True,constraints=None,constrain_final=False):

		assert n_images > 2

		if traj is None:
			traj = f"{self.outdir}/neb_rcs.traj"

		from ase.neb import NEB
		from ase.optimize import MDMin
		if constrain_final:
			from ase.constraints import FixAtoms

		start = FakeAtoms(start)
		final = FakeAtoms(final)

		mout.varOut("NEB Start",start.rcs)
		mout.varOut("NEB Final",final.rcs)

		images = [start]
		for i in range(n_images-2):
			image = start.copy()
			image.calc = FakeSurfaceCalculator(self.pmf,suppress_warnings=suppress_warnings)
			image.set_constraints(constraints[0],constraints[1],constraints[2])
			images.append(image)
		images.append(final)

		neb = NEB(images)
		neb.interpolate()

		import matplotlib.pyplot as plt

		fig,ax = plt.subplots()

		interpolated_rcs = [[p[0] for p in a.get_positions()] for a in images]

		fig,ax = plt.subplots()
		plt.plot(self.energy_along_path(interpolated_rcs))
		plt.savefig(f"{self.outdir}/neb_interpolated.png")
		plt.close()

		dyn = optimizer(neb, trajectory=traj)
		dyn.run(fmax=fmax)

def pullx2rco(infile,output,n_coords,skip=1):
	"""convert a pullx to an rco file"""

	import mout
	import os

	mout.out(f"{infile} --> {os.path.basename(output)}")

	f_input = open(infile,'r')
	f_output = open(output,'w')

	mins = [None for i in range(n_coords)]
	maxs = [None for i in range(n_coords)]

	for line in f_input:

		if line[0] in ["#","@"]:
			count = 0
			continue

		if count%skip != 0:
			count += 1
			continue

		# the new RCO should contain
		# time, N*(restraint)

		split_line=line.split()

		time = float(split_line[0])
		RCs = []

		for i in range(n_coords):
			x = float(split_line[1+8*i])*10
			RCs.append(x)
			if mins[i] is None or x < mins[i]:
				mins[i] = x
			if maxs[i] is None or x > maxs[i]:
				maxs[i] = x

		f_output.write(f'{time:.4f}')
		for x in RCs:
			f_output.write(f' {x:.6f}')

		f_output.write("\n")

		count += 1

	# return reference information
	window_centres = []
	for i in range(n_coords):
		window_centres.append(float(split_line[2+8*i]))

	f_input.close()
	f_output.close()

	return window_centres, mins, maxs

def prepare_metafile(outdir,winkeys,n_coords,force,skip=1):
	import mout
	import mcol
	import os
	
	mout.debugHeader("prepare_metafile()")

	os.system(f'mkdir -p {outdir}')

	metafile = open(f'{outdir}/metafile','w')

	mout.varOut("metafile",f'{outdir}/metafile',valCol=mcol.file)
	mout.varOut("force",f'{force}',unit="kcal/mol/Length^2",valCol=mcol.arg)

	all_mins = [None for i in range(n_coords)]
	all_maxs = [None for i in range(n_coords)]

	mout.headerOut("Creating amber-style RCO's")
	for key in winkeys:

		# create the 'amber-style' RCO's
		centres, mins, maxs = pullx2rco(f'{key}/pullx.xvg',f'{outdir}/{key}.rco',n_coords,skip=skip)

		for i,(this_min,this_max) in enumerate(zip(mins,maxs)):
			if all_mins[i] is None or this_min < all_mins[i]:
				all_mins[i] = this_min
			if all_maxs[i] is None or this_max > all_maxs[i]:
				all_maxs[i] = this_max

		# add to the metafile
		# the file should contain
		# RCO_path, N*(RN_window_centre, RN_force)

		metafile.write(f'{outdir}/{key}.rco')
		for x in centres:
			metafile.write(f' {x*10:.6f} {force:.2f}')
		metafile.write('\n')

	for i in range(n_coords):
		mout.varOut(f"RC{i+1} range: ",f'({all_mins[i]} - {all_maxs[i]})',unit="Å")
## use scipy to find minima and maxima in the 4D PMF