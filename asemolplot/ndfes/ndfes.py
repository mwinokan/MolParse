
def pullx2rco(infile,output,n_coords):
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
			continue

		# the new RCO should contain
		# time, N*(restraint)

		split_line=line.split()

		time = float(split_line[0])
		RCs = []

		for i in range(n_coords):
			x = float(split_line[1+8*i])
			RCs.append(x)
			if mins[i] is None or x < mins[i]:
				mins[i] = x
			if maxs[i] is None or x > maxs[i]:
				maxs[i] = x

		f_output.write(f'{time:.4f}')
		for x in RCs:
			f_output.write(f' {x:.6f}')

		f_output.write("\n")

	# return reference information
	window_centres = []
	for i in range(n_coords):
		window_centres.append(float(split_line[2+8*i]))

	f_input.close()
	f_output.close()

	return window_centres, mins, maxs

def prepare_metafile(outdir,winkeys,n_coords,force):
	import mout
	import mcol
	import os
	
	mout.debugOut("prepare_metafile()")

	os.system(f'mkdir -p {outdir}')

	metafile = open(f'{outdir}/metafile','w')

	mout.varOut("metafile",f'{outdir}/metafile',valCol=mcol.file)
	mout.varOut("force",f'{force}',unit="kcal/mol/Length^2",valCol=mcol.arg)

	all_mins = [None for i in range(n_coords)]
	all_maxs = [None for i in range(n_coords)]

	mout.headerOut("Creating amber-style RCO's")
	for key in winkeys:

		# create the 'amber-style' RCO's
		centres, mins, maxs = pullx2rco(f'{key}/pullx.xvg',f'{outdir}/{key}.rco',n_coords)

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
			metafile.write(f' {x:.6f} {force:.2f}')
		metafile.write('\n')

	for i in range(n_coords):
		mout.varOut(f"RC{i+1} range: ",f'({all_mins[i]} - {all_maxs[i]})',unit="nm")

def run_ndfes(outdir,n_coords,method="vfep",binlength=0.2,temperature=300,script="eureka.sh"):
	import os
	import mout
	import mcol

	mout.varOut("metafile",f'{outdir}/metafile',valCol=mcol.file)
	mout.varOut("checkpoint",f'{outdir}/checkpoint',valCol=mcol.file)
	mout.varOut("method",f'{method}',valCol=mcol.arg)
	mout.varOut("temperature",f'{temperature}',unit="K",valCol=mcol.arg)
	mout.varOut("script",f'{script}',valCol=mcol.file)

	amp_path = os.environ['MWAMPPATH']
	script = f'{amp_path}/asemolplot/ndfes/{script}'

	length_str = ""
	for i in range(n_coords):
		length_str += f"-w {binlength} "

	mout.headerOut("Running ndfes...")
	os.system(f'{script} --{method} -c {outdir}/checkpoint {length_str} {outdir}/metafile -t {temperature} > {outdir}/ndfes.log')

def plot_2d_fes(outdir,dims):
	try:
		import ndfes
	except:
		mout.errorOut("Missing ndfes library. Try sourcing $MWSHPATH/load_ndfes.sh",fatal=True)
	import mout
	import mcol

	mout.varOut("checkpoint",f'{outdir}/checkpoint',valCol=mcol.file)

	ndgrid = ndfes.NDUniformRegularGrid()
	for d in dims:
		ndgrid.add_dim_using_delta(d[0],d[1],d[2])

	pts = ndgrid.get_list_of_tuples()

	meta,pmf = generate_pmf_object(f'{outdir}/checkpoint',f'{outdir}/metafile')

	vals = pmf(pts)

	# Print the PMF to a file
	fh = open(f'{outdir}/2dfes.dat',"w")
	for pt,val in zip(pts,vals):
	    # if pt[0] == -1.0: fh.write("\n")

		# collapse onto 2d

		if val > 1000:
			continue

		fh.write("%15.5f %15.5f %20.10e\n"%(pt[0]-pt[1],pt[2]-pt[3],val))
	
	fh.close()

def generate_pmf_object(checkpoint,metafile):
	try:
		import ndfes
	except:
		mout.errorOut("Missing ndfes library. Try sourcing $MWSHPATH/load_ndfes.sh",fatal=True)

	meta = ndfes.Metafile(metafile)
	pmf = ndfes.interpolator([checkpoint],meta)

	return meta,pmf

def plot_window_pmf(outdir,i1,i2):
	import mout
	import mcol
	import matplotlib.pyplot as plt
	import numpy as np

	meta,pmf = generate_pmf_object(f'{outdir}/checkpoint',f'{outdir}/metafile')

	cs = meta.get_centers()

	# reactant free energy
	print(pmf([cs[i1]]))

	# product free energy
	print(pmf([cs[i2]]))

	cs = cs[0:5]

	# reshape into 1d

	vals = pmf(cs)

	x = []
	y = []

	for c,v in zip(cs,vals):
		x.append(c[0]-c[1]+c[2]-c[3])
		print(c)
		if v < 1000:
			y.append(v)
		else:
			y.append(None)

	fig,ax = plt.subplots()

	plt.scatter(x,y)

	plt.savefig(f"{outdir}/window_pmf.png")

	plt.close()
