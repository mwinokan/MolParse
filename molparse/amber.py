
"""Experimental and undocumented. Speak to Max."""

from .restraint import Restraint

def write_amber_restraints(restraint_list,filename_prefix,zfill=2,filename_suffix=".RST"):

	import mcol
	import mout
	import os

	os.system("mkdir -p "+os.path.dirname(filename_prefix))

	assert isinstance(restraint_list,list)

	# get all the restraint values
	all_values = []
	for res in restraint_list:
		assert isinstance(res,Restraint)
		values = res.values()
		assert values is not None
		all_values.append(values)
		assert len(values) == len(all_values[0])

	num_windows = len(all_values[0])

	for i in range(num_windows):

		i_str = str(i).zfill(zfill)

		filename = filename_prefix+i_str+filename_suffix
		
		# print(i)
		# print(filename)

		rst_buffer = "# "+filename+"\n"

		for restraint in restraint_list:
			rst_buffer += restraint.amber_block(i)

		# print(rst_buffer)

		out_rst = open(filename,"w")
		out_rst.write(rst_buffer)
		out_rst.close()
		mout.out("File written to "+mcol.file+filename)

	import pickle

	pickle_file = os.path.dirname(filename_prefix)+"/amp_res.pkl"
	with open(pickle_file, 'wb') as file:
		pickle.dump(restraint_list, file)
	mout.out("Restraints dumped to to "+mcol.file+pickle_file)

		# exit()

def umbrella_plotter(filenames,bins=20,subdir=None,show_level=2):

	import os
	import mplot

	from . import signal

	big_ydata = []
	labels=[]

	if subdir is not None:
		os.system("mkdir -p "+subdir)

	for file in filenames:

		print(file)

		xdata,ydata = signal.parseDat(file,num_columns=2,pre_strip=True)

		label=os.path.basename(file).replace(".rco",'')
		labels.append(label)

		plotfile=None
		if subdir is not None:
			plotfile=subdir+"/"+label+".png"

		if show_level > 1:
			show=True
		else:
			show=False

		mplot.hist1D(ydata,show=show,xlab="Reaction Coordinate",ylab="Frequency",bins=bins,title=label,filename=plotfile)

		big_ydata.append(ydata)

	if show_level > 0:
		show=True
	else:
		show=False

	plotfile=None
	if subdir is not None:
		plotfile=subdir+"/allwindows.png"

	mplot.graph2D(xdata,big_ydata,show=show,filename=plotfile,ytitles=labels,xlab="MD Step",ylab="Reaction Coordinate")

def umbrella_helper_2dist(atoms,weights,coord_range,num_windows,force_constant,harmonic_width,subdir=None,samples=1000,graph=False):
	
	import mcol
	import mout
	import os
	import mplot

	assert len(atoms) == 4
	assert len(weights) == 2

	mout.headerOut("Atoms")

	for i,atom in enumerate(atoms):

		mout.varOut("Atom #"+str(i+1),[atom.name,atom.residue,atom.pdb_index],integer=True,list_length=False)

	mout.headerOut("Windows")
	mout.varOut("#Windows",num_windows,valCol=mcol.arg)

	centres=[]
	
	for i in range(num_windows):

		centres.append(coord_range[0]+i*(coord_range[1]-coord_range[0])/(num_windows-1))

		mout.varOut("Window #"+str(i+1)+" Centre",centres[-1],precision=4)

	mout.headerOut("Amber Restraint File:")

	if subdir is not None:
		os.system("mkdir -p "+subdir)

	restraints = []

	for i,centre in enumerate(centres):

		end = "\n"

		rst_buffer = "# window_"+str(i+1)+".RST"+end
		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[0].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[2].pdb_index)+","
		rst_buffer += str(atoms[3].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[0])+","
		rst_buffer += str(weights[1])+","+end

		r2 = centre
		r3 = centre

		r1 = centre - harmonic_width
		r4 = centre + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end+end

		restraints.append({'centre':centre,
						   'r1':r1,
						   'r2':r2,
						   'r3':r3,
						   'r4':r4,
						   'rk2':rk2,
						   'rk3':rk3})

		if subdir is not None:
		    out_rst = open(subdir+"/window_"+str(i+1).zfill(2)+".RST","w")
		    out_rst.write(rst_buffer)
		    out_rst.close()
		    mout.out("File written to "+mcol.file+subdir+"/window_"+str(i+1).zfill(2)+".RST")
	
	xdata = []
	big_ydata = []

	for i in range(samples):

		# print(i,coord_range[0]+i*(coord_range[1]-coord_range[0])/(samples-1))
		x = coord_range[0]+(i*1.4/(samples-1)-0.2)*(coord_range[1]-coord_range[0])
		xdata.append(x)

	for restraint in restraints:

		ydata = []

		r1 = restraint['r1']
		r2 = restraint['r2']
		r3 = restraint['r3']
		r4 = restraint['r4']
		rk2 = restraint['rk2']
		rk3 = restraint['rk3']
		
		for x in xdata:

			y = restraint_potential(x,r1,r2,r3,r4,rk2,rk3)

			ydata.append(y)

		big_ydata.append(ydata)

	# if graph:
	if subdir is not None:
		graphfile=subdir+"/allwindows.png"
	else:
		graphfile=None

	if subdir is not None and graph is not None:
		mplot.graph2D(xdata,big_ydata,show=graph,ymax=120,ymin=-5,filename=graphfile)

	if subdir is not None:

		pot_buffer = "# restraint potentials"+end

		pot_buffer += str(x) + " "

		for restraint in restraints:
			r1 = restraint['r1']
			r2 = restraint['r2']
			r3 = restraint['r3']
			r4 = restraint['r4']
			rk2 = restraint['rk2']
			rk3 = restraint['rk3']
			y = restraint_potential(x,r1,r2,r3,r4,rk2,rk3)
			pot_buffer += str(y) + " "

		out_dat = open(subdir+"/allwindows.dat","w")
		out_dat.write(pot_buffer)
		out_dat.close()

def umb_rst_2prot(atoms,weights,coord_range,num_windows,force_constant,harmonic_width,subdir=None,samples=1000,graph=False,adiab_windows=True,fix_angle=None,fix_third=None):
	
	import mcol
	import mout
	import os
	import mplot

	assert len(atoms) == 6
	assert len(weights) == 2

	# 2x (donor-hydrogen...acceptor)

	mout.headerOut("Atoms")

	for i,atom in enumerate(atoms):

		mout.varOut("Atom #"+str(i+1),[atom.name,atom.residue,atom.pdb_index],integer=True,list_length=False)

	mout.headerOut("Windows")
	mout.varOut("#Windows",num_windows,valCol=mcol.arg)

	centres=[]
	
	for i in range(num_windows):

		centres.append(coord_range[0]+i*(coord_range[1]-coord_range[0])/(num_windows-1))

		mout.varOut("Window #"+str(i+1)+" Centre",centres[-1],precision=4)

	mout.headerOut("Amber Restraint File:")

	if subdir is not None:
		os.system("mkdir -p "+subdir)

	restraints = []

	### RC:
	### 		dist(donor-hydrogen)	  [1.0->2.0]
	### 	  - dist(hydrogen...acceptor) [2.0->1.0]
	### 	  =			   				  [-1.0->1.0]

	mout.array2file("centres.dat",centres)

	if fix_angle is not None:
		mout.warningOut("Angle restraints are active!")
	if fix_third is not None:
		mout.warningOut("Third bond restraints are active!")

	for i,centre in enumerate(centres):

		end = "\n"

		rst_buffer = "# window_"+str(i+1)+".RST"+end

		# First reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[0].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[2].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[0])+","
		rst_buffer += str(-weights[0])+","+end

		r2 = centre
		r3 = centre

		r1 = centre - harmonic_width
		r4 = centre + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end

		# Second reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[3].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[5].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[1])+","
		rst_buffer += str(-weights[1])+","+end

		r2 = centre
		r3 = centre

		r1 = centre - harmonic_width
		r4 = centre + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end

		restraints.append({'centre':centre,
						   'r1':r1,
						   'r2':r2,
						   'r3':r3,
						   'r4':r4,
						   'rk2':rk2,
						   'rk3':rk3})

		###### ANGLE FIX ######

		# fix_angle should be a list with [atom,atom,atom,angle]

		if fix_angle is not None:

			# try to catch some errors
			assert isinstance(fix_angle,list)
			assert len(fix_angle) == 4

			# construct the restraint

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(fix_angle[0].pdb_index)+","
			rst_buffer += str(fix_angle[1].pdb_index)+","
			rst_buffer += str(fix_angle[2].pdb_index)+","+end

			r2 = fix_angle[-1]
			r3 = 180

			r1 = r2 - 20
			r4 = r2 + 20

			rk2 = force_constant
			rk3 = force_constant

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

		###### THIRD BOND FIX ######

		# fix_third should be a list with [atom,atom,atom]
		# [donor,H,acceptor]

		if fix_third is not None:

			# try to catch some errors
			assert isinstance(fix_third,list)
			assert len(fix_third) == 3

			# construct the restraints

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(fix_third[0].pdb_index)+","
			rst_buffer += str(fix_third[1].pdb_index)+","+end

			r2 = 0.9
			r3 = 1.5

			r1 = 0.0
			r4 = 2.4

			rk2 = force_constant
			rk3 = force_constant

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(fix_third[1].pdb_index)+","
			rst_buffer += str(fix_third[2].pdb_index)+","+end

			r2 = 1.4
			r3 = 2.4

			r1 = 0.5
			r4 = 3.3

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

		rst_buffer += end

		if subdir is not None:
		    out_rst = open(subdir+"/window_"+str(i+1).zfill(2)+".RST","w")
		    out_rst.write(rst_buffer)
		    out_rst.close()
		    mout.out("File written to "+mcol.file+subdir+"/window_"+str(i+1).zfill(2)+".RST")

	if adiab_windows:
		
		for i,centre in enumerate(centres):

			end = "\n"

			rst_buffer = "# adiab_"+str(i+1)+".RST"+end

			# First reaction coord

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(atoms[0].pdb_index)+","
			rst_buffer += str(atoms[1].pdb_index)+","
			rst_buffer += str(atoms[1].pdb_index)+","
			rst_buffer += str(atoms[2].pdb_index)+","+end

			rst_buffer += "rstwt="
			rst_buffer += str(weights[0])+","
			rst_buffer += str(-weights[0])+","+end

			r2 = centre
			r3 = centre

			r1 = centre - harmonic_width*10
			r4 = centre + harmonic_width*10

			rk2 = force_constant*50
			rk3 = force_constant*50

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

			# Second reaction coord

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(atoms[3].pdb_index)+","
			rst_buffer += str(atoms[4].pdb_index)+","
			rst_buffer += str(atoms[4].pdb_index)+","
			rst_buffer += str(atoms[5].pdb_index)+","+end

			rst_buffer += "rstwt="
			rst_buffer += str(weights[1])+","
			rst_buffer += str(-weights[1])+","+end

			r2 = centre
			r3 = centre

			r1 = centre - harmonic_width*10
			r4 = centre + harmonic_width*10

			rk2 = force_constant*50
			rk3 = force_constant*50

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end+end

			if subdir is not None:
			    out_rst = open(subdir+"/adiab_"+str(i+1).zfill(2)+".RST","w")
			    out_rst.write(rst_buffer)
			    out_rst.close()
			    mout.out("File written to "+mcol.file+subdir+"/adiab_"+str(i+1).zfill(2)+".RST")
	
	xdata = []
	big_ydata = []

	for i in range(samples):

		# print(i,coord_range[0]+i*(coord_range[1]-coord_range[0])/(samples-1))
		x = coord_range[0]+(i*1.4/(samples-1)-0.2)*(coord_range[1]-coord_range[0])
		xdata.append(x)

	for restraint in restraints:

		ydata = []

		r1 = restraint['r1']
		r2 = restraint['r2']
		r3 = restraint['r3']
		r4 = restraint['r4']
		rk2 = restraint['rk2']
		rk3 = restraint['rk3']
		
		for x in xdata:

			y = restraint_potential(x,r1,r2,r3,r4,rk2,rk3)

			ydata.append(y)

		big_ydata.append(ydata)

	# if graph:
	if subdir is not None:
		graphfile=subdir+"/allwindows.png"
	else:
		graphfile=None

	if subdir is not None and graph is not None:
		mplot.graph2D(xdata,big_ydata,show=graph,ymax=120,ymin=-5,filename=graphfile)

	if subdir is not None:

		pot_buffer = "# restraint potentials"+end

		pot_buffer += str(x) + " "

		for restraint in restraints:
			r1 = restraint['r1']
			r2 = restraint['r2']
			r3 = restraint['r3']
			r4 = restraint['r4']
			rk2 = restraint['rk2']
			rk3 = restraint['rk3']
			y = restraint_potential(x,r1,r2,r3,r4,rk2,rk3)
			pot_buffer += str(y) + " "

		out_dat = open(subdir+"/allwindows.dat","w")
		out_dat.write(pot_buffer)
		out_dat.close()

def umb_rst_2prot_new(atoms,weights,coord_range,num_windows,force_constant,harmonic_width,subdir=None,samples=1000,graph=False,adiab_windows=True,fix_angle=None,fix_third=None,add_len=None):

	import mcol
	import mout
	import os

	assert len(atoms) == 6
	assert len(weights) == 2

	# 2x (donor-hydrogen...acceptor)

	mout.headerOut("Atoms")

	for i,atom in enumerate(atoms):

		mout.varOut("Atom #"+str(i+1),[atom.name,atom.residue,atom.pdb_index],integer=True,list_length=False)

	mout.headerOut("Windows")
	mout.varOut("#Windows",num_windows,valCol=mcol.arg)

	print(coord_range)

	if add_len is not None:
		coord_range[0] -= 0.5*add_len
		coord_range[1] += 0.5*add_len
		coord_range[2] -= 0.5*add_len
		coord_range[3] += 0.5*add_len
	print(coord_range)

	centres=[]
	
	for i in range(num_windows):

		centre_1 = coord_range[0]+i*(coord_range[1]-coord_range[0])/(num_windows-1)
		centre_2 = coord_range[2]+i*(coord_range[3]-coord_range[2])/(num_windows-1)

		centres.append([centre_1,centre_2])

		mout.varOut("Window #"+str(i+1)+" Centres",centres[-1],precision=4)

	mout.headerOut("Amber Restraint File:")

	if subdir is not None:
		os.system("mkdir -p "+subdir)

	restraints = []

	### RC:
	### 		dist(donor-hydrogen)	  [1.0->2.0]
	### 	  - dist(hydrogen...acceptor) [2.0->1.0]
	### 	  =			   				  [-1.0->1.0]

	mout.array2file("centres.dat",centres)

	if fix_angle is not None:
		mout.warningOut("Angle restraints are active!")
	if fix_third is not None:
		mout.warningOut("Third bond restraints are active!")

	for i,centre in enumerate(centres):

		end = "\n"

		rst_buffer = "# window_"+str(i+1)+".RST"+end

		# First reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[0].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[2].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[0])+","
		rst_buffer += str(-weights[0])+","+end

		r2 = centre[0]
		r3 = centre[0]

		r1 = centre[0] - harmonic_width
		r4 = centre[0] + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end

		# Second reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[3].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[5].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[1])+","
		rst_buffer += str(-weights[1])+","+end

		r2 = centre[1]
		r3 = centre[1]

		r1 = centre[1] - harmonic_width
		r4 = centre[1] + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end

		# restraints.append({'centre':centre,
		# 				   'r1':r1,
		# 				   'r2':r2,
		# 				   'r3':r3,
		# 				   'r4':r4,
		# 				   'rk2':rk2,
		# 				   'rk3':rk3})

		###### ANGLE FIX ######

		# fix_angle should be a list with [atom,atom,atom,angle]

		if fix_angle is not None:

			# try to catch some errors
			assert isinstance(fix_angle,list)
			assert len(fix_angle) == 4

			# construct the restraint

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(fix_angle[0].pdb_index)+","
			rst_buffer += str(fix_angle[1].pdb_index)+","
			rst_buffer += str(fix_angle[2].pdb_index)+","+end

			r2 = fix_angle[-1]
			r3 = 180

			r1 = r2 - 20
			r4 = r2 + 20

			rk2 = force_constant
			rk3 = force_constant

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

		###### THIRD BOND FIX ######

		# fix_third should be a list with [atom,atom,atom]
		# [donor,H,acceptor]

		if fix_third is not None:

			# try to catch some errors
			assert isinstance(fix_third,list)
			assert len(fix_third) == 3

			# construct the restraints

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(fix_third[0].pdb_index)+","
			rst_buffer += str(fix_third[1].pdb_index)+","+end

			if add_len is not None:
				assert isinstance(add_len,float) or isinstance(add_len,int)
				this_add = add_len
			else:
				this_add = 0.0

			r2 = 0.9
			r3 = 1.5

			r1 = 0.0
			r4 = 2.4

			rk2 = force_constant
			rk3 = force_constant

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(fix_third[1].pdb_index)+","
			rst_buffer += str(fix_third[2].pdb_index)+","+end

			r2 = 1.4 + this_add
			r3 = 2.4 + this_add

			r1 = 0.5 + this_add
			r4 = 3.3 + this_add

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

		if add_len is not None:

			# Additional restraints to extend hydrogen bonds

			assert isinstance(add_len,float) or isinstance(add_len,int)

			# construct the restraints

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(atoms[0].pdb_index)+","
			rst_buffer += str(atoms[2].pdb_index)+","+end

			r2 = 2.5+add_len
			r3 = 3.0+add_len

			r1 = r2 - 1.0
			r4 = r3 + 1.0

			rk2 = force_constant
			rk3 = force_constant

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(atoms[3].pdb_index)+","
			rst_buffer += str(atoms[5].pdb_index)+","+end

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

		rst_buffer += end

		if subdir is not None:
		    out_rst = open(subdir+"/window_"+str(i+1).zfill(2)+".RST","w")
		    out_rst.write(rst_buffer)
		    out_rst.close()
		    mout.out("File written to "+mcol.file+subdir+"/window_"+str(i+1).zfill(2)+".RST")

def umb_rst_2prot_1rc(atoms,weights,coord_range,num_windows,force_constant,harmonic_width,subdir=None,samples=1000,graph=False,adiab_windows=False,fix_length=False):

	import mcol
	import mout
	import os

	assert len(atoms) == 6
	assert len(weights) == 2

	# 2x (donor-hydrogen...acceptor)

	mout.headerOut("Atoms")

	for i,atom in enumerate(atoms):

		mout.varOut("Atom #"+str(i+1),[atom.name,atom.residue,atom.pdb_index],integer=True,list_length=False)

	mout.headerOut("Windows")
	mout.varOut("#Windows",num_windows,valCol=mcol.arg)

	centres=[]
	
	for i in range(num_windows):

		centres.append(coord_range[0]+i*(coord_range[1]-coord_range[0])/(num_windows-1))

		mout.varOut("Window #"+str(i+1)+" Centre",centres[-1],precision=4)

	mout.headerOut("Amber Restraint File:")

	if subdir is not None:
		os.system("mkdir -p "+subdir)

	restraints = []

	### RC:
	### 		dist(donor-hydrogen)	  [1.0->2.0]
	### 	  - dist(hydrogen...acceptor) [2.0->1.0]
	### 	  =			   				  [-1.0->1.0]

	for i,centre in enumerate(centres):

		end = "\n"

		rst_buffer = "# window_"+str(i+1)+".RST"+end

		# First reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[0].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[2].pdb_index)+","
		rst_buffer += str(atoms[3].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[5].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[0])+","
		rst_buffer += str(-weights[0])+","
		rst_buffer += str(weights[1])+","
		rst_buffer += str(-weights[1])+","+end

		r2 = centre
		r3 = centre

		r1 = centre - harmonic_width
		r4 = centre + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end

		# Output only dummy-restraints

		rk2 = 0.0
		rk3 = 0.0

		rst_buffer += "&rst"+end
		rst_buffer += "iat="
		rst_buffer += str(atoms[0].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[2].pdb_index)+","+end
		rst_buffer += "rstwt="
		rst_buffer += str(weights[0])+","
		rst_buffer += str(-weights[0])+","+end
		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		rst_buffer += "/"+end

		# Second reaction coord

		rst_buffer += "&rst"+end
		rst_buffer += "iat="
		rst_buffer += str(atoms[3].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[5].pdb_index)+","+end
		rst_buffer += "rstwt="
		rst_buffer += str(weights[1])+","
		rst_buffer += str(-weights[1])+","+end
		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		rst_buffer += "/"+end

		# Restrain donor-acceptor H-bond lengths

		if fix_length:

			r1 = 0.0
			r2 = 2.0
			r3 = 4.0
			r4 = 6.0

			rk2 = force_constant * 2
			rk3 = force_constant * 2

			rst_buffer += "&rst"+end
			rst_buffer += "iat="
			rst_buffer += str(atoms[0].pdb_index)+","
			rst_buffer += str(atoms[2].pdb_index)+","+end
			rst_buffer += "rstwt="
			rst_buffer += str(weights[1])+","+end
			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			rst_buffer += "/"+end

			rst_buffer += "&rst"+end
			rst_buffer += "iat="
			rst_buffer += str(atoms[3].pdb_index)+","
			rst_buffer += str(atoms[5].pdb_index)+","+end
			rst_buffer += "rstwt="
			rst_buffer += str(weights[1])+","+end
			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			rst_buffer += "/"+end

		rst_buffer += end

		if subdir is not None:
		    out_rst = open(subdir+"/window_"+str(i+1).zfill(2)+".RST","w")
		    out_rst.write(rst_buffer)
		    out_rst.close()
		    mout.out("File written to "+mcol.file+subdir+"/window_"+str(i+1).zfill(2)+".RST")

	if adiab_windows:
		mout.warningOut("Adiabatic windows not implemented")

def umb_rst_2prot_asy(atoms,weights,coord_range,num_windows,force_constant,harmonic_width,subdir=None,samples=1000,graph=False,adiab_windows=True,reverse=False):

	import mcol
	import mout
	import os

	num_windows = num_windows//2+1

	assert len(atoms) == 6
	assert len(weights) == 2

	# 2x (donor-hydrogen...acceptor)

	mout.headerOut("Atoms")

	for i,atom in enumerate(atoms):

		mout.varOut("Atom #"+str(i+1),[atom.name,atom.residue,atom.pdb_index],integer=True,list_length=False)

	mout.headerOut("Windows")
	mout.varOut("#Windows",num_windows,valCol=mcol.arg)

	centres=[]
	
	for i in range(num_windows):

		centres.append(coord_range[0]+i*(coord_range[1]-coord_range[0])/(num_windows-1))

		mout.varOut("Window #"+str(i+1)+" Centre",centres[-1],precision=4)

	centres_2d=[]
	
	if not reverse:
		for centre in centres:
			centres_2d.append([centre,centres[0]])
		for centre in centres:
			if centre == centres[0]: continue
			centres_2d.append([centres[-1],centre])
	else:
		for centre in centres:
			centres_2d.append([centres[0],centre])
		for centre in centres:
			if centre == centres[0]: continue
			centres_2d.append([centre,centres[-1]])

	mout.headerOut("Amber Restraint File:")

	# mout.varOut("centres",centres_2d,precision=2,sf=False)

	if subdir is not None:
		os.system("mkdir -p "+subdir)

	restraints = []

	### RC:
	### 		dist(donor-hydrogen)	  [1.0->2.0]
	### 	  - dist(hydrogen...acceptor) [2.0->1.0]
	### 	  =			   				  [-1.0->1.0]

	for i,centre_pair in enumerate(centres_2d):

		end = "\n"

		rst_buffer = "# window_"+str(i+1)+".RST"+end

		# First reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[0].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[1].pdb_index)+","
		rst_buffer += str(atoms[2].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[0])+","
		rst_buffer += str(-weights[0])+","+end

		r2 = centre_pair[0]
		r3 = centre_pair[0]

		r1 = centre_pair[0] - harmonic_width
		r4 = centre_pair[0] + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end

		# Second reaction coord

		rst_buffer += "&rst"+end

		rst_buffer += "iat="
		rst_buffer += str(atoms[3].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[4].pdb_index)+","
		rst_buffer += str(atoms[5].pdb_index)+","+end

		rst_buffer += "rstwt="
		rst_buffer += str(weights[1])+","
		rst_buffer += str(-weights[1])+","+end

		r2 = centre_pair[1]
		r3 = centre_pair[1]

		r1 = centre_pair[1] - harmonic_width
		r4 = centre_pair[1] + harmonic_width

		rk2 = force_constant
		rk3 = force_constant

		rst_buffer += "r1="+str(r1)+","+end
		rst_buffer += "r2="+str(r2)+","+end
		rst_buffer += "r3="+str(r3)+","+end
		rst_buffer += "r4="+str(r4)+","+end
		rst_buffer += "rk2="+str(rk2)+","+end
		rst_buffer += "rk3="+str(rk3)+","+end
		
		rst_buffer += "/"+end+end

		restraints.append({'centre_pair':centre_pair,
						   'r1':r1,
						   'r2':r2,
						   'r3':r3,
						   'r4':r4,
						   'rk2':rk2,
						   'rk3':rk3})

		if subdir is not None:
		    out_rst = open(subdir+"/window_"+str(i+1).zfill(2)+".RST","w")
		    out_rst.write(rst_buffer)
		    out_rst.close()
		    mout.out("File written to "+mcol.file+subdir+"/window_"+str(i+1).zfill(2)+".RST")

	if adiab_windows:
		
		for i,centre_pair in enumerate(centres_2d):

			end = "\n"

			rst_buffer = "# adiab_"+str(i+1)+".RST"+end

			# First reaction coord

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(atoms[0].pdb_index)+","
			rst_buffer += str(atoms[1].pdb_index)+","
			rst_buffer += str(atoms[1].pdb_index)+","
			rst_buffer += str(atoms[2].pdb_index)+","+end

			rst_buffer += "rstwt="
			rst_buffer += str(weights[0])+","
			rst_buffer += str(-weights[0])+","+end

			r2 = centre_pair[0]
			r3 = centre_pair[0]

			r1 = centre_pair[0] - harmonic_width*10
			r4 = centre_pair[0] + harmonic_width*10

			rk2 = force_constant*50
			rk3 = force_constant*50

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end

			# Second reaction coord

			rst_buffer += "&rst"+end

			rst_buffer += "iat="
			rst_buffer += str(atoms[3].pdb_index)+","
			rst_buffer += str(atoms[4].pdb_index)+","
			rst_buffer += str(atoms[4].pdb_index)+","
			rst_buffer += str(atoms[5].pdb_index)+","+end

			rst_buffer += "rstwt="
			rst_buffer += str(weights[1])+","
			rst_buffer += str(-weights[1])+","+end

			r2 = centre_pair[1]
			r3 = centre_pair[1]

			r1 = centre_pair[1] - harmonic_width*10
			r4 = centre_pair[1] + harmonic_width*10

			rk2 = force_constant*50
			rk3 = force_constant*50

			rst_buffer += "r1="+str(r1)+","+end
			rst_buffer += "r2="+str(r2)+","+end
			rst_buffer += "r3="+str(r3)+","+end
			rst_buffer += "r4="+str(r4)+","+end
			rst_buffer += "rk2="+str(rk2)+","+end
			rst_buffer += "rk3="+str(rk3)+","+end
			
			rst_buffer += "/"+end+end

			if subdir is not None:
			    out_rst = open(subdir+"/adiab_"+str(i+1).zfill(2)+".RST","w")
			    out_rst.write(rst_buffer)
			    out_rst.close()
			    mout.out("File written to "+mcol.file+subdir+"/adiab_"+str(i+1).zfill(2)+".RST")
	
	# import mplot

	# xdata = []
	# big_ydata = []

	# for i in range(samples):

	# 	# print(i,coord_range[0]+i*(coord_range[1]-coord_range[0])/(samples-1))
	# 	x = coord_range[0]+(i*1.4/(samples-1)-0.2)*(coord_range[1]-coord_range[0])
	# 	xdata.append(x)

	# for restraint in restraints:

	# 	ydata = []

	# 	r1 = restraint['r1']
	# 	r2 = restraint['r2']
	# 	r3 = restraint['r3']
	# 	r4 = restraint['r4']
	# 	rk2 = restraint['rk2']
	# 	rk3 = restraint['rk3']
		
	# 	for x in xdata:

	# 		y = restraint_potential(x,r1,r2,r3,r4,rk2,rk3)

	# 		ydata.append(y)

	# 	big_ydata.append(ydata)

	# # if graph:
	# if subdir is not None:
	# 	graphfile=subdir+"/allwindows.png"
	# else:
	# 	graphfile=None

	# if subdir is not None and graph is not None:
	# 	mplot.graph2D(xdata,big_ydata,show=graph,ymax=120,ymin=-5,filename=graphfile)

	# if subdir is not None:

	# 	pot_buffer = "# restraint potentials"+end

	# 	pot_buffer += str(x) + " "

	# 	for restraint in restraints:
	# 		r1 = restraint['r1']
	# 		r2 = restraint['r2']
	# 		r3 = restraint['r3']
	# 		r4 = restraint['r4']
	# 		rk2 = restraint['rk2']
	# 		rk3 = restraint['rk3']
	# 		y = restraint_potential(x,r1,r2,r3,r4,rk2,rk3)
	# 		pot_buffer += str(y) + " "

	# 	out_dat = open(subdir+"/allwindows.dat","w")
	# 	out_dat.write(pot_buffer)
	# 	out_dat.close()

def restraint_potential(x,r1,r2,r3,r4,rk2,rk3):

	if x < r1:
		return lin1(x,r1,r2,rk2)

	if (x >= r1 and x <= r2):
		return par2(x,r2,rk2)

	if (x >= r3 and x <= r4):
		return par3(x,r3,rk3)

	if x > r4:
		return lin4(x,r3,r4,rk3)

	return 0.0

def lin1(x,r1,r2,rk2):
	return 2*rk2*(r1-r2)*(x-r1)+par2(r1,r2,rk2)

def lin4(x,r3,r4,rk3):
	return 2*rk3*(r4-r3)*(x-r4)+par3(r4,r3,rk3)

def par2(x,r2,rk2):
	return pow(rk2*(x-r2),2)

def par3(x,r3,rk3):
	return pow(rk3*(x-r3),2)

def prep4amber(system):

	import mcol
	import mout

	# Fix histidine residue names
	if "HSD" in str(system.residues):
		i = system.rename_residues("HSD","HID",verbosity=0)
		mout.out("Fixed names of "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"HSD")

	###### Fix names common to many animo acids

	# HN -> H
	aa_hnfix_residues = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
	  					  "GLY", "HID", "ILE", "LEU", "LYS", "MET", "PHE", 
	  					  "SER", "THR", "TRP", "TYR", "VAL"]

	i=0
	for resname in aa_hnfix_residues:
		i += system.rename_atoms("HN","H",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HN->H"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")

	# HBN->HBN+1
	aa_hbfix_residues = [ "ARG", "ASN", "ASP", "CYS", "GLN", "GLU","HID", "LEU",
						  "LYS", "MET", "PHE", "PRO", "SER", "TRP", "TYR"]

	i=0
	for resname in aa_hbfix_residues:
		system.rename_atoms("HB2","HB3",res_filter=resname,verbosity=0)
		i += system.rename_atoms("HB1","HB2",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HBn->HBn+1"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")

	# HGN->HGN+1
	i=0
	for resname in [ "ALA", "ARG", "GLN", "GLU", "LYS", "MET", "PRO" ]:
		system.rename_atoms("HG2","HG3",res_filter=resname,verbosity=0)
		i += system.rename_atoms("HG1","HG2",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HGn->HGn+1"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")

	# HDN->HDN+1
	i=0
	for resname in [ "ARG", "LYS", "PRO"]:
		system.rename_atoms("HD2","HD3",res_filter=resname,verbosity=0)
		i += system.rename_atoms("HD1","HD2",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HDn->HDn+1"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")

	###### Fix atom names in specific amino acids

	# Fix LEU atom names
	if "LEU" in str(system.residues):
		system.rename_atoms("HT1","H1",res_filter="LEU",verbosity=0)
		system.rename_atoms("HT2","H2",res_filter="LEU",verbosity=0)
		i = system.rename_atoms("HT3","H3",res_filter="LEU",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"LEU")

	# Fix SER atom names
	if "SER" in str(system.residues):
		i = system.rename_atoms("HG1","HG",res_filter="SER",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"SER")

	# Fix LYS atom names
	if "LYS" in str(system.residues):
		system.rename_atoms("HE2","HE3",res_filter="LYS",verbosity=0)
		i = system.rename_atoms("HE1","HE2",res_filter="LYS",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"LYS")

	# Fix GLY atom names
	if "GLY" in str(system.residues):
		system.rename_atoms("HA2","HA3",res_filter="GLY",verbosity=0)
		i = system.rename_atoms("HA1","HA2",res_filter="GLY",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"GLY")

	# Fix ILE atom names
	if "ILE" in str(system.residues):
		system.rename_atoms("HG12","HG13",res_filter="ILE",verbosity=0)
		system.rename_atoms("HG11","HG12",res_filter="ILE",verbosity=0)
		# system.rename_atoms("2HG1","HG13",res_filter="ILE",verbosity=1)
		# system.rename_atoms("1HG1","HG12",res_filter="ILE",verbosity=1)
		system.rename_atoms("CD","CD1",res_filter="ILE",verbosity=0)
		system.rename_atoms("HD1","HD11",res_filter="ILE",verbosity=0)
		system.rename_atoms("HD2","HD12",res_filter="ILE",verbosity=0)
		i = system.rename_atoms("HD3","HD13",res_filter="ILE",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"ILE")

	# Fix CYS atom names
	if "CYS" in str(system.residues):
		i = system.rename_atoms("HG1","HG",res_filter="CYS",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"CYS")

	# Fix ARG atom names
	if "ARG" in str(system.residues):
		system.rename_atoms("OT1","O",res_filter="ARG",verbosity=0)
		i = system.rename_atoms("OT2","OXT",res_filter="ARG",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"ARG")

	# Fix thymine atom names
	if "DT" in str(system.residues):
		system.rename_atoms("H51","H71",res_filter="DT",verbosity=0)
		system.rename_atoms("H52","H72",res_filter="DT",verbosity=0)
		i = system.rename_atoms("H53","H73",res_filter="DT",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"DT")

	# HCAP
	i=0
	for resname in [ "DC3", "DA3", "DT3", "DG3"]:
		i += system.rename_atoms("HCAP","HO3'",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HCAP->HO3'"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of 3TER")

	# Fix ATP atom names
	if "ATP" in str(system.residues):
		system.rename_atoms("O5'","O5*",res_filter="ATP",verbosity=0)
		system.rename_atoms("C5'","C5*",res_filter="ATP",verbosity=0)
		system.rename_atoms("C4'","C4*",res_filter="ATP",verbosity=0)
		system.rename_atoms("O4'","O4*",res_filter="ATP",verbosity=0)
		system.rename_atoms("C1'","C1*",res_filter="ATP",verbosity=0)
		system.rename_atoms("C3'","C3*",res_filter="ATP",verbosity=0)
		system.rename_atoms("O3'","O3*",res_filter="ATP",verbosity=0)
		system.rename_atoms("C2'","C2*",res_filter="ATP",verbosity=0)
		system.rename_atoms("O2'","O2*",res_filter="ATP",verbosity=0)
		system.rename_atoms("H4'","H40",res_filter="ATP",verbosity=0)
		system.rename_atoms("H1'","H10",res_filter="ATP",verbosity=0)
		system.rename_atoms("H5'","H50",res_filter="ATP",verbosity=0)
		system.rename_atoms("H5''","H51",res_filter="ATP",verbosity=0)
		system.rename_atoms("H8","H80",res_filter="ATP",verbosity=0)
		system.rename_atoms("H61","H60",res_filter="ATP",verbosity=0)
		system.rename_atoms("H62","H61",res_filter="ATP",verbosity=0)
		system.rename_atoms("H2''","H20",res_filter="ATP",verbosity=0)
		i = system.rename_atoms("H3T","H30",res_filter="ATP",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"ATP")

def amber2charmm(system):

	import mcol
	import mout

	# Fix histidine residue names
	if "HSD" in str(system.residues):
		i = system.rename_residues("HSD","HID",verbosity=0)
		mout.out("Fixed names of "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"HSD")

	###### Fix names common to many animo acids

	# H -> HN
	aa_hnfix_residues = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
	  					  "GLY", "HID", "ILE", "LEU", "LYS", "MET", "PHE", 
	  					  "SER", "THR", "TRP", "TYR", "VAL"]

	i=0
	for resname in aa_hnfix_residues:
		i += system.rename_atoms("H","HN",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"H->HN"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")


	# HBN->HBN-1
	aa_hbfix_residues = [ "ARG", "ASN", "ASP", "CYS", "GLN", "GLU","HID", "LEU",
						  "LYS", "MET", "PHE", "PRO", "SER", "TRP", "TYR"]

	i=0
	for resname in aa_hbfix_residues:
		i += system.rename_atoms("HB2","HB1",res_filter=resname,verbosity=0)
		system.rename_atoms("HB3","HB2",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HBn->HBn-1"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")

	# HGN->HGN-1
	i=0
	for resname in [ "ALA", "ARG", "GLN", "GLU", "LYS", "MET", "PRO" ]:
		i += system.rename_atoms("HG2","HG1",res_filter=resname,verbosity=0)
		system.rename_atoms("HG3","HG2",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HGn->HGn-1"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")
	
	# HDN->HDN-1
	i=0
	for resname in [ "ARG", "LYS", "PRO"]:
		i += system.rename_atoms("HD2","HD1",res_filter=resname,verbosity=0)
		system.rename_atoms("HD3","HD2",res_filter=resname,verbosity=0)
	mout.out(mcol.arg+"HDn->HDn-1"+mcol.clear+" rename performed in "+
			 mcol.result+str(i)+mcol.clear+" instances of amino acids")

	###### Fix atom names in specific amino acids

	# Fix LEU atom names
	if "LEU" in str(system.residues):
		system.rename_atoms("H1","HT1",res_filter="LEU",verbosity=0)
		system.rename_atoms("H2","HT2",res_filter="LEU",verbosity=0)
		i = system.rename_atoms("H3","HT3",res_filter="LEU",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"LEU")

	# Fix SER atom names
	if "SER" in str(system.residues):
		i = system.rename_atoms("HG","HG1",res_filter="SER",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"SER")

	# Fix LYS atom names
	if "LYS" in str(system.residues):
		i = system.rename_atoms("HE2","HE1",res_filter="LYS",verbosity=0)
		system.rename_atoms("HE3","HE2",res_filter="LYS",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"LYS")

	# Fix GLY atom names
	if "GLY" in str(system.residues):
		i = system.rename_atoms("HA2","HA1",res_filter="GLY",verbosity=0)
		system.rename_atoms("HA3","HA2",res_filter="GLY",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"GLY")

	# Fix ILE atom names
	if "ILE" in str(system.residues):
		system.rename_atoms("HG12","HG11",res_filter="ILE",verbosity=0)
		system.rename_atoms("HG13","HG12",res_filter="ILE",verbosity=0)
		# system.rename_atoms("2HG1","HG13",res_filter="ILE",verbosity=1)
		# system.rename_atoms("1HG1","HG12",res_filter="ILE",verbosity=1)
		system.rename_atoms("CD1","CD",res_filter="ILE",verbosity=0)
		system.rename_atoms("HD11","HD1",res_filter="ILE",verbosity=0)
		system.rename_atoms("HD12","HD2",res_filter="ILE",verbosity=0)
		i = system.rename_atoms("HD13","HD3",res_filter="ILE",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"ILE")

	# Fix CYS atom names
	if "CYS" in str(system.residues):
		i = system.rename_atoms("HG","HG1",res_filter="CYS",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"CYS")

	# Fix ARG atom names
	if "ARG" in str(system.residues):
		system.rename_atoms("OT1","O",res_filter="ARG",verbosity=0)
		i = system.rename_atoms("OXT","OT2",res_filter="ARG",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"ARG")

	# Fix thymine atom names
	if "DT" in str(system.residues):
		system.rename_atoms("H71","H51",res_filter="DT",verbosity=0)
		system.rename_atoms("H72","H52",res_filter="DT",verbosity=0)
		i = system.rename_atoms("H73","H53",res_filter="DT",verbosity=0)
		mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"DT")

	# # Fix HIS atom names
	# if "HIS" in str(system.residues):
	# 	i = system.rename_atoms("HD2","HD1",res_filter="HIS",verbosity=0)
	# 	mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"HIS")

	# # Fix ATP atom names
	# if "ATP" in str(system.residues):
	# 	system.rename_atoms("O5'","O5*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("C5'","C5*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("C4'","C4*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("O4'","O4*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("C1'","C1*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("C3'","C3*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("O3'","O3*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("C2'","C2*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("O2'","O2*",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H4'","H40",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H1'","H10",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H5'","H50",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H5''","H51",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H8","H80",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H61","H60",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H62","H61",res_filter="ATP",verbosity=0)
	# 	system.rename_atoms("H2''","H20",res_filter="ATP",verbosity=0)
	# 	i = system.rename_atoms("H3T","H30",res_filter="ATP",verbosity=0)
	# 	mout.out("Fixed atom names in "+mcol.result+str(i)+mcol.clear+" instances of "+mcol.arg+"ATP")

def parseRST(filename,key,convert_to_float=True):

	import re

	file = open(filename, "r")

	output = []

	for line in file:
		result = re.search(key+'=(.*),', line)
		if result is not None:
			value = result.group(1)
			if convert_to_float: value = float(value)
			output.append(value)

			# print(result.group(1))

	return output
