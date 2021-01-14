
import mcol
import mout
import os
import mplot

from . import signal


def umbrella_plotter(filenames,bins=20,subdir=None,show_level=2):

	# xdata = []
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

		# mplot.graph2D(xdata,ydata)
		mplot.hist1D(ydata,show=show,xlab="Reaction Coordinate",ylab="Frequency",bins=bins,title=label,filename=plotfile)

		big_ydata.append(ydata)

	# mout.varOut("x",xdata)
	# mout.varOut("y",big_ydata)

	if show_level > 0:
		show=True
	else:
		show=False

	plotfile=None
	if subdir is not None:
		plotfile=subdir+"/allwindows.png"

	mplot.graph2D(xdata,big_ydata,show=show,filename=plotfile,ytitles=labels,xlab="MD Step",ylab="Reaction Coordinate")

def umbrella_helper_2dist(atoms,weights,coord_range,num_windows,force_constant,harmonic_width,subdir=None,samples=1000,graph=False):

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
	
	import mplot

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

	if graph:
		mplot.graph2D(xdata,big_ydata,ymax=120,ymin=-5)

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

