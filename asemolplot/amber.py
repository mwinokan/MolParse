
import mcol
import mout

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

def dnacleanup(system):

DG H2
DG N6
DG H61
DG H62
DC N9
DC N2
DC H21
DC H22
DC H1
DC O6
DC N7
DC C8
DC H8
DG H6
DG O2
DG H3
DG O4
DG C5M
DG H71
DG H72
DG H73
DC N9
DC N2
DC H21
DC H22
DC H1
DC O6
DC N7
DC C8
DC H8
DG H6
DG H5
DG O2
DG N4
DG H41
DG H42
DC N9
DC N7
DC C8
DC H8
DC H2
DC N6
DC H61
DC H62
DC H3
DC O4
DC C5M
DC H71
DC H72
DC H73
DG  H6
DG  O2
DG  H3
DG  O4
DG  C5M
DG  H71
DG  H72
DG  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
DC3 H3
DC3 O4
DC3 C5M
DC3 H71
DC3 H72
DC3 H73
DC  N9
DC  N7
DC  C8
DC  H8
DC  H2
DC  N6
DC  H61
DC  H62
DG  H2
DG  N6
DG  H61
DG  H62
DG  H6
DG  O2
DG  H3
DG  O4
DG  C5M
DG  H71
DG  H72
DG  H73
DC  N9
DC  N2
DC  H21
DC  H22
DC  H1
DC  O6
DC  N7
DC  C8
DC  H8
DG  H6
DG  H5
DG  O2
DG  N4
DG  H41
DG  H42
DC  N9
DC  N7
DC  C8
DC  H8
DC  H2
DC  N6
DC  H61
DC  H62
DG  H6
DG  H5
DG  O2
DG  N4
DG  H41
DG  H42
DC  H3
DC  O4
DC  C5M
DC  H71
DC  H72
DC  H73
