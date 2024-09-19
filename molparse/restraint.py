class Restraint:
    def __init__(self, atoms, type=None, weights=None, force_constant=None):

        self.atoms = atoms
        self.range = None
        self.force_constant = force_constant
        self.weights = weights
        self._values = None

        if type == None:
            self.type = self.determine_type()

        self.name = self.determine_name()

    def determine_type(self):
        import mout

        """ Restraint types:
                Two-Atom: Distance							Implemented
                Three-Atom: Angle							Implemented
                Four-Atom: Torsional / J-Coupling
                Four-Atom: Generalised Distance (Transfer)	Implemented
                Five-Atom: Plane-Point
                Six-Atom: Generalised Distance
                Eight-Atom: Plane-Plane						Implemented
        """

        if len(self.atoms) == 2:
            return "distance"

        elif len(self.atoms) == 3:
            return "angle"

        elif len(self.atoms) == 4:
            assert len(self.weights) == 2
            return "transfer"

        elif len(self.atoms) == 8:
            return "plane-plane"

        else:
            mout.errorOut("Could not determine restraint type", fatal=True)

    def determine_name(self):
        if self.type == "distance":
            return "dist(" + self.atoms[0].name + "-" + self.atoms[1].name + ")"
        elif self.type == "transfer":
            return (
                "RC("
                + self.atoms[0].name
                + "-"
                + self.atoms[1].name
                + ","
                + self.atoms[2].name
                + "-"
                + self.atoms[3].name
                + ")"
            )
        elif self.type == "angle":
            return (
                "ang("
                + self.atoms[0].name
                + ","
                + self.atoms[1].name
                + ","
                + self.atoms[2].name
                + ")"
            )
        elif self.type == "plane-plane":
            return "plane-plane(...)"
        else:
            return "restraint"

    def value(self, system, scale_angle=None):
        import numpy as np
        import mout

        this_atoms = []
        for atom in self.atoms:
            this_atom = system.atoms[atom.index]
            assert str(this_atom) == str(atom)
            this_atoms.append(this_atom)

        if self.type == "distance":
            return np.linalg.norm(this_atoms[0].np_pos - this_atoms[1].np_pos)

        elif self.type == "transfer":
            dist1 = np.linalg.norm(this_atoms[0].np_pos - this_atoms[1].np_pos)
            dist2 = np.linalg.norm(this_atoms[2].np_pos - this_atoms[3].np_pos)
            return self.weights[0] * dist1 + self.weights[1] * dist2

        elif self.type == "angle":
            vec1 = this_atoms[0].np_pos - this_atoms[1].np_pos
            vec2 = this_atoms[2].np_pos - this_atoms[1].np_pos
            dot = np.dot(vec1, vec2)
            norms = np.linalg.norm(vec1) * np.linalg.norm(vec2)
            if scale_angle == None:
                return np.arccos(dot / norms) / np.pi * 180
            else:
                return np.arccos(dot / norms) / np.pi * scale_angle

        elif self.type == "plane-plane":
            # (r1 - r2) × (r3 - r4) and (r5 - r6) × (r7 - r8)
            vec1 = this_atoms[0].np_pos - this_atoms[1].np_pos
            vec2 = this_atoms[2].np_pos - this_atoms[3].np_pos
            vec3 = this_atoms[4].np_pos - this_atoms[5].np_pos
            vec4 = this_atoms[6].np_pos - this_atoms[7].np_pos
            cross1 = np.cross(vec1, vec2)
            cross2 = np.cross(vec3, vec4)
            dot = np.dot(cross1, cross2)
            norms = np.linalg.norm(cross1) * np.linalg.norm(cross2)
            if scale_angle == None:
                return np.arccos(dot / norms) / np.pi * 180
            else:
                return np.arccos(dot / norms) / np.pi * scale_angle

        else:
            mout.errorOut("Unsupported restraint type", fatal=True)

    def values(self, system_list=None, scale_angle=None):
        if system_list is None:
            return self._values
        else:
            from .system import System

            assert isinstance(system_list, list)
            values = []
            for system in system_list:
                assert isinstance(system, System)
                values.append(self.value(system, scale_angle=scale_angle))
            self._values = values
            return values

    def equi_values(self, n):
        assert self._values is not None

        from scipy.interpolate import interp1d

        x = [i for i in range(len(self._values))]

        inter_func = interp1d(x, self._values)

        equi = []

        for i in range(n):
            this_x = x[0] + i * (x[-1] - x[0]) / (n - 1)

            equi.append(inter_func(this_x))

        return equi

    def set_values(self, values):
        self._values = values

    def fetch_atoms(self, new_system, new_residues, res_map=None, verbosity=1):

        import mout

        count = 0

        new_atoms = []

        for atom in self.atoms:

            found = False

            # print(atom.residue,new_residues)

            if res_map is not None:
                for pair in res_map:
                    # if atom.residue in pair[0]:
                    # 	print(atom.residue,pair[0])
                    # 	atom.residue = pair[1]
                    # 	print(atom.residue,pair[0])
                    # 	found_map = True
                    # 	break
                    if atom.residue == pair[1]:
                        # print(atom.residue)
                        atom.residue = pair[0]
                        # print(atom.residue)
                        break

            for res in new_residues:
                if res.name == atom.residue:
                    new_atoms.append(res.get_atom(atom.name))
                    found = True
                    count += 1

            if not found:
                mout.varOut("atoms", self.atoms)
                mout.varOut("new_residues", new_residues)
                mout.errorOut(
                    "Could not find match for "
                    + atom.name
                    + " ("
                    + atom.residue
                    + ")"
                    + " in residues "
                    + str(new_residues),
                    fatal=True,
                    code="Restraint.fetch_atoms",
                )

        if verbosity > 0:
            mout.out(
                "Retrieved "
                + str(count)
                + " atoms from "
                + new_system.name
                + "'s "
                + str(new_residues)
            )

        self.atoms = new_atoms

    def amber_block(self, index, harmonic_width=10.0, debug=False):

        value = self.values()[index]

        end = "\n"

        rst_buffer = "&rst" + end

        rst_buffer += "iat="

        for atom in self.atoms:
            rst_buffer += str(atom.pdb_index) + ","

        if debug:
            rst_buffer += " #"

            for atom in self.atoms:
                rst_buffer += str(atom.name) + ", "

        rst_buffer += end

        if self.type == "transfer":
            rst_buffer += "rstwt="
            rst_buffer += str(self.weights[0]) + ","
            rst_buffer += str(self.weights[1]) + "," + end

        r1 = value - harmonic_width
        r2 = value
        r3 = value
        r4 = value + harmonic_width

        rk2 = self.force_constant
        rk3 = self.force_constant

        rst_buffer += "r1=" + str(r1) + "," + end
        rst_buffer += "r2=" + str(r2) + "," + end
        rst_buffer += "r3=" + str(r3) + "," + end
        rst_buffer += "r4=" + str(r4) + "," + end
        rst_buffer += "rk2=" + str(rk2) + "," + end
        rst_buffer += "rk3=" + str(rk3) + "," + end

        rst_buffer += "/" + end

        return rst_buffer

    def copy(self):
        import copy

        return copy.deepcopy(self)

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name
