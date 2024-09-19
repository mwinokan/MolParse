from ase.calculators.general import Calculator


class DL_POLY(Calculator):
    import mout
    import mcol
    from ase.calculators.calculator import compare_atoms, PropertyNotImplementedError
    import os

    def __init__(
        self,
        field,
        control,
        workdir="DL_POLY",
        name="DL_POLY",
        output="OUTPUT",
        verbosity=1,
        debug=False,
        num_procs=8,
        backup=False,
    ):

        if debug:
            mout.debugOut(name + ".__init__")

        # Manage arguments
        self.name = name
        self.workdir = workdir
        self.field = field
        self.control = control
        self.output = output
        self.debug = debug
        if debug:
            verbosity = 99
        self.verbosity = verbosity
        self.energy_zero = None
        self.energy_free = None
        self.forces = None
        self.old_atoms = None
        self.run_count = 0
        self.num_procs = num_procs
        self.backup = backup

        self.field_path = self.workdir + "/" + self.field
        self.control_path = self.workdir + "/" + self.control
        self.output_path = self.workdir + "/" + self.output
        self.revcon_path = self.workdir + "/REVCON"
        self.history_path = self.workdir + "/HISTORY"
        self.statis_path = self.workdir + "/STATIS"
        self.auto_config = "auto.CONFIG"
        self.auto_config_path = self.workdir + "/" + self.auto_config

        self.results = {"energy": self.energy_zero, "forces": self.forces}

        self.implemented_properties = ["energy", "forces"]

        if self.verbosity > 1:
            mout.varOut("IO: Workdir", self.workdir, valCol=mcol.file)
            mout.varOut("IO: Field", self.field, valCol=mcol.file)
            mout.varOut("IO: Control", self.control, valCol=mcol.file)
            mout.varOut("IO: Output", self.output, valCol=mcol.file)

        self.clean_directory(full=True, verbosity=verbosity - 2)

    def update(self, atoms):
        if self.debug:
            mout.debugOut(self.name + ".update")

        # Check if an update is required:

        if self.old_atoms is None:  # No calculation has been run yet

            if self.verbosity > 1:
                mout.out("Calculation needed (first run)")

            self.calculate(atoms)

        else:  # Previous calculation data exists

            # Get the names of changed properties
            changed_properties = compare_atoms(atoms, self.old_atoms)

            # If there are changes, run DL_POLY SPE
            if len(changed_properties) > 0:
                if self.verbosity > 1:
                    mout.out("Calculation needed (atoms have changed)")
                    mout.varOut("Changed Properties", changed_properties)
                self.calculate(atoms)
            else:
                if self.verbosity > 1:
                    mout.out("No calculation needed (atoms are unchanged)")

    def clean_directory(self, full=False, verbosity=1):
        if self.debug:
            mout.debugOut(self.name + ".clean_directory")

        if full:  # Remove the work directory
            if verbosity > 0:
                os.system("rm -rfv " + self.workdir)
            else:
                os.system("rm -rf " + self.workdir)
        else:
            if verbosity > 0:  # Only remove output files
                os.system("rm -f " + self.revcon_path)
                os.system("rm -f " + self.output_path)
                os.system("rm -f " + self.history_path)
                os.system("rm -f " + self.statis_path)
            else:
                os.system("rm -fv " + self.revcon_path)
                os.system("rm -fv " + self.output_path)
                os.system("rm -fv " + self.history_path)
                os.system("rm -fv " + self.statis_path)

    def calculate(self, atoms):
        from molparse import write, read
        from .dl_poly_py import SPE

        if self.debug:
            mout.debugOut(self.name + ".calculate")

        # Check pre-requisites
        if self.field is None:
            mout.errorOut("No FIELD specified")
        if self.control is None:
            mout.errorOut("No CONTROL specified")

        if self.run_count > 0:
            if self.backup:
                os.system(
                    "cp "
                    + self.auto_config_path
                    + " "
                    + self.workdir
                    + "/backup"
                    + str(self.run_count)
                    + ".CONFIG"
                )
        else:
            os.system("mkdir -p " + self.workdir)

        # Write new system configuration in DL_POLY format
        write(self.auto_config, atoms, format="dlp4", verbosity=self.verbosity - 1)

        # Get rid of old outputs
        self.clean_directory(full=False, verbosity=self.verbosity - 2)

        # Call DL_POLY
        e_tot = SPE(
            control=self.control,
            config=self.auto_config,
            field=self.field,
            workdir=self.workdir,
            output=self.output,
            num_procs=self.num_procs,
            verbosity=self.verbosity - 1,
        )

        # Clean up files in root directory
        os.system("rm " + self.auto_config)
        os.system("rm OUTPUT")

        # Store the atoms on which the calculation was run
        self.old_atoms = atoms.copy()

        # Get the trajectory
        history = read(
            self.history_path, format="dlp-history", verbosity=self.verbosity - 1
        )

        # Store the calculation results
        self.energy_zero = e_tot
        self.energy_free = e_tot
        self.forces = history.get_forces()
        self.results = {"energy": self.energy_zero, "forces": self.forces}

        # Iterate the run_count
        self.run_count += 1

    def get_forces(self, atoms):
        if self.debug:
            mout.debugOut(self.name + ".get_forces")
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        raise PropertyNotImplementedError

    def get_dipole(self, atoms):
        raise PropertyNotImplementedError

    def get_charges(self, atoms):
        raise PropertyNotImplementedError

    def get_magmom(self, atoms):
        raise PropertyNotImplementedError

    def get_magmoms(self, atoms):
        raise PropertyNotImplementedError

    def reset(self):
        self.energy_zero = None
        self.forces = None
        self.old_atoms = None
        self.energy_free = None
        self.run_count = 0
        self.results = {"energy": self.energy_zero, "forces": self.forces}

    # def set_field(self,field):
    #   if self.debug: mout.debugOut(self.name+".set_field")
    #   self.field = field

    # def set_control(self,control):
    #   if self.debug: mout.debugOut(self.name+".set_control")
    #   self.control = control
