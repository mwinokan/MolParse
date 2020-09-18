
from ase.calculators.general import Calculator
from ase.calculators.calculator import compare_atoms, PropertyNotImplementedError

from .dl_poly_py import SPE

from asemolplot import write, read

import mout
import mcol

import os

class DL_POLY(Calculator):

  def __init__(self,field,control,workdir="DL_POLY",verbosity=1,debug=False):

    if debug: mout.debugOut("DL_POLY.__init__")

    self.field = field
    self.control = control
    self.workdir = workdir
    self.verbosity = verbosity
    self.debug = debug
    self.energy_zero = None
    self.forces = None
    self.old_atoms = None
    self.run_count = 0

    self.results = {'energy':self.energy_zero,'forces':self.forces}

    self.implemented_properties=['energy','forces']

    self.clean_directory(full=True,verbosity=verbosity-2)

  def update(self,atoms):
    if self.debug: mout.debugOut("DL_POLY.update")

    # Check if an update is required:

    if self.old_atoms is None: # No calculation has been run yet
      self.calculate(atoms)
    else:                      # Previous calculation data exists

      # Get the names of changed properties
      changed_properties = compare_atoms(atoms,self.old_atoms)

      # If there are changes, run DL_POLY SPE
      if len(changed_properties) > 0:
        if self.verbosity > 1:
          mout.varOut("Changed Properties",changed_properties)
        self.calculate(atoms)

  def set_field(self,field):
    if self.debug: mout.debugOut("DL_POLY.set_field")
    self.field = field

  def set_control(self,control):
    if self.debug: mout.debugOut("DL_POLY.set_control")
    self.control = control

  def clean_directory(self,full=False,verbosity=1):
    if self.debug: mout.debugOut("DL_POLY.clean_directory")

    if full:
      if verbosity > 0:
        os.system("rm -rfv "+self.workdir)
      else:
        os.system("rm -rf "+self.workdir)
    else:
      if verbosity > 0:
        os.system("rm -fv "+self.workdir+"/REVCON")
        os.system("rm -fv "+self.workdir+"/OUTPUT")
        os.system("rm -fv "+self.workdir+"/HISTORY")
        os.system("rm -fv "+self.workdir+"/STATIS")
      else:
        os.system("rm -fv "+self.workdir+"/REVCON")
        os.system("rm -fv "+self.workdir+"/OUTPUT")
        os.system("rm -fv "+self.workdir+"/HISTORY")
        os.system("rm -fv "+self.workdir+"/STATIS")

  def calculate(self,atoms):
    if self.debug: mout.debugOut("DL_POLY.calculate")

    # Check pre-requisites
    if self.field is None:
      mout.errorOut("No FIELD specified")
    if self.control is None:
      mout.errorOut("No CONTROL specified")

    if self.run_count > 0:
      os.system("cp "+self.workdir+"/auto.CONFIG "+self.workdir+"/backup"+str(self.run_count)+".CONFIG")

    # Write new system configuration in DL_POLY format
    write("auto.CONFIG",atoms,format="dlp4",verbosity=self.verbosity-1)

    # Get rid of old outputs
    self.clean_directory(full=False,verbosity=self.verbosity-2)

    # Call DL_POLY
    e_tot = SPE(control=self.control,
                config="auto.CONFIG",
                field=self.field,
                workdir=self.workdir,
                verbosity=self.verbosity-1)

    # Store the atoms on which the calculation was run
    self.old_atoms = atoms.copy()

    # Get the trajectory
    history = read(self.workdir+"/HISTORY",format="dlp-history",verbosity=self.verbosity-1)
    
    # Store the calculation results
    self.energy_zero = e_tot
    self.forces = history.get_forces()
    self.results = {'energy':self.energy_zero,'forces':self.forces}

    # Iterate the run_count
    self.run_count += 1

  def get_forces(self,atoms):
    if self.debug: mout.debugOut("DL_POLY.get_forces")
    self.update(atoms)
    return self.forces

  def get_stress(self,atoms):
    raise PropertyNotImplementedError

  def get_dipole(self,atoms):
    raise PropertyNotImplementedError

  def get_charges(self,atoms):
    raise PropertyNotImplementedError

  def get_magmom(self,atoms):
    raise PropertyNotImplementedError

  def get_magmoms(self,atoms):
    raise PropertyNotImplementedError

  def reset(self):
    self.energy_zero = None
    self.forces = None
    self.old_atoms = None
