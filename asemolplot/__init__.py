"""

ASE Molecular Plotting


Extensions to ASE for tasks involving molecular simulations

It is recommended to import as the acronym:

	import asemolplot as amp

Key features:

* Parsing coordinate files into System objects:

	system = amp.parsePDB(filestr)
	system = amp.parseGRO(filestr)

* Parsing of delimited text files:

	see help(amp.signal.parseDat)

* Writing System or ase.Atoms objects or lists thereof:

	amp.write(filestr,object)

* An object oriented approach to molecular system hierarchies:

	amp.System
	    |__ amp.Chain
	            |__ amp.Residue
	                    |__ amp.Atom

	See also:

	help(amp.System), help(amp.Chain), etc.

* Wrappers for common ASE methods (with AMP class support):

	amp.read --> ase.io.read
	amp.write --> ase.io.write
	amp.view --> ase.visualize.view

* Functions to render images and animations of ASE Atoms:

	with matplotlib:
		amp.makeImage, 
		amp.makeImages, 
		amp.makeAnimation

	with PoVRay:
		amp.makePovImage, 
		amp.makePovImages, 
		amp.makePovAnimation

* Methods relating to DNA:

	see help(amp.dna)

* Finding the stationary points of a double well potential

	see help(amp.tunnel.find_barrier_stationary_points)

* Graphing of common molecular properties:

	see help(amp.graphing)

* Analysis of common molecular properties:

	see help(amp.analysis)

* Tunnel

"""

# Custom IO
from .io import read
from .io import write
from .io import parsePDB
from .io import parseGRO

# Manipulation
from .manipulate import interpolate
from .manipulate import auto_rotate

# Custom Classes
from .system import System
from .chain import Chain
from .residue import Residue
from .atom import Atom
from .restraint import Restraint

# ase.visualize wrappers
from .gui import view

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
from .compare import compareSystems
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

# from . import ndfes
