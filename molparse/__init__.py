"""

Molecular Parsing

A python package for parsing, modifying, and analysis of molecular structure files. 

It is recommended to import as the acronym:

	import molparse as mp

Key features:

* Parsing coordinate files into System objects:

	system = mp.parsePDB(filestr)
	system = mp.parseGRO(filestr)

* Viewing the molecular hierarchy of a system:

	mp.tree(system)

* Parsing of delimited text files:

	see help(mp.signal.parseDat)

* Writing System or ase.Atoms objects or lists thereof:

	mp.write(filestr,object)

* An object oriented approach to molecular system hierarchies:

	mp.System
	    |__ mp.Chain
	            |__ mp.Residue
	                    |__ mp.Atom

	See also:

	help(mp.System), help(mp.Chain), etc.

* Wrappers for common ASE methods (with mp class support):

	mp.read --> ase.io.read
	mp.write --> ase.io.write
	mp.view --> ase.visualize.view

* Functions to render images and animations of ASE Atoms:

	with matplotlib:
		mp.makeImage, 
		mp.makeImages, 
		mp.makeAnimation

	with PoVRay:
		mp.makePovImage, 
		mp.makePovImages, 
		mp.makePovAnimation

* Methods relating to DNA:

	see help(mp.dna)

* Finding the stationary points of a double well potential

	see help(mp.tunnel.find_barrier_stationary_points)

* Graphing of common molecular properties:

	see help(mp.graphing)

* Analysis of common molecular properties:

	see help(mp.analysis)

* Tunnel

"""

# from .version import import_checks
# import_checks()

# Custom IO
from .io import read
from .io import write
from .io import parse
from .io import parsePDB
from .io import writePDB
from .io import parseGRO
from .io import parseXYZ
from .io import modifyPDB

from .signal import parseDat

# CLI-Viewer
# from .tree import tree

# Manipulation
from .manipulate import interpolate
from .manipulate import auto_rotate

# Custom Classes
from .group import AtomGroup
from .system import System
from .chain import Chain
from .residue import Residue
from .atom import Atom
from .restraint import Restraint
from .network import Network

# ase.visualize wrappers
from .gui import view
from .go import plot3d

# ase standard rendering
from .plot import makeImage
from .plot import makeImages
from .plot import makeAnimation

# PoV-Ray rendering
from .povplot import loadPov
from .povplot import makePovImage
from .povplot import makePovImages
from .povplot import makePovAnimation
from .povplot import crop

# Analysis
from .analysis import bondLengthStats
from .analysis import bondAngleStats
from .analysis import getCentreOfMass
from .analysis import getRMSD
from .analysis import getAngle
from .analysis import getDisplacement
from .analysis import getBondLabel
from .analysis import getAtomLabel

# Structure Comparison
# from .compare import compareSystems
from .compare import euclid_dist

# Graphing
from .graphing import showFigs
from .graphing import graphEnergy
from .graphing import graphForces
from .graphing import graphDisplacement
from .graphing import graphBondLength
from .graphing import graphBondAngle
from .graphing import graphBondVibSpec

# Console IO
from .console_io import printEnergy # needs documentation

# XVG Parsing
from .xvg import parseXVG

# NDX Parsing
from .ndx import parseNDX

# Amber specific
from .amber import prep4amber
from .amber import amber2charmm
from .amber import parseRST
from .amber import umbrella_helper_2dist
from .amber import umb_rst_2prot
from .amber import umb_rst_2prot_new
from .amber import umb_rst_2prot_1rc
from .amber import umb_rst_2prot_asy
from .amber import umbrella_plotter
from .amber import write_amber_restraints

# Umbrella Integration
# from .umbrella_integrate import umbrella_integrate

# File conversion
from .convert import pdb2traj # needs documentation
from .convert import gro2traj # needs documentation
from .convert import xyz2traj # needs documentation

# import subpackages

from . import styles
from . import units
# from . import dl_poly
from . import dna
from . import signal
from . import tunnel
from . import mutate
from . import monte
from . import rdkit

from .amino import alphabet as amino_alphabet
from .amino import longname as amino_longname

from . import hijack
