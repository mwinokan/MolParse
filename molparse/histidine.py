from .amino import AminoAcid


class Histidine(AminoAcid):
    """Class for Histidine Amino Acid Residue

    These objects should not be created by the user,
    but constructed automatically when parsing a
    coordinate file via amp.parsePDB or otherwise"""

    def __init__(self, name, index, number, chain):
        super(Histidine, self).__init__(name, index, number, chain)

        assert len(name) == 3
        assert name in ["HIS", "HID", "HIE", "HSE", "HSD", "HSP"]

    @property
    def ND(self):
        """Get Nitrogen-Delta Atom"""
        return self.get_atom("ND1", verbosity=False)

    @property
    def CG(self):
        """Get Carbon-Gamma Atom"""
        return self.get_atom("CG", verbosity=False)

    @property
    def CD(self):
        """Get Carbon-Delta Atom"""
        return self.get_atom("CD2", verbosity=False)

    @property
    def NE(self):
        """Get Nitrogen-Epsilon Atom"""
        return self.get_atom("NE2", verbosity=False)

    @property
    def HD(self):
        """Get Hydrogen-Delta Atom"""
        return self.get_atom("HD1", verbosity=False)

    @property
    def HE(self):
        """Get Hydrogen-Epsilon Atom"""
        return self.get_atom("HE2", verbosity=False)

    @property
    def has_HD(self):
        """Get Hydrogen-Delta Atom"""
        return bool(self.HD)

    @property
    def has_HE(self):
        """Get Hydrogen-Epsilon Atom"""
        return bool(self.HE)

    def protonate(self, delta=True, epsilon=True, verbosity=1):
        """Mutate this histidine residue to another protonation state"""

        import mout
        import mcol
        import numpy as np
        from .atom import Atom

        if verbosity > 0:
            if delta and epsilon:
                mout.headerOut(
                    f"Protonating {mcol.arg}{self.longname}{mcol.clear + mcol.bold} --> {mcol.arg} Histidine (Protonated)"
                )
            elif delta and not epsilon:
                mout.headerOut(
                    f"Protonating {mcol.arg}{self.longname}{mcol.clear + mcol.bold} --> {mcol.arg} Histidine (Delta)"
                )
            elif epsilon and not delta:
                mout.headerOut(
                    f"Protonating {mcol.arg}{self.longname}{mcol.clear + mcol.bold} --> {mcol.arg} Histidine (Epsilon)"
                )
            else:
                mout.errorOut(f"Unsupported protonation {delta=} {epsilon=}")
                return

        # warning
        if not delta and not epsilon:
            mout.warningOut("Neither site protonated!")

        # add delta-hydrogen
        if delta and not self.has_HD:

            if verbosity > 1:
                mout.warningOut("Protonating delta site")

            length = 1.00  # Ang

            site = self.ND
            opposite = (self.NE + self.CD) / 2
            vec = site.np_pos - opposite
            vec = length * vec / np.linalg.norm(vec)
            pos = site + vec

            if not epsilon and self.has_HE:
                atom = self.HE
                atom.name = "HD1"
            else:
                atom = Atom("HD1")
                self.add_atom(atom)
            atom.position = pos

        # add epsilon-hydrogen
        if epsilon and not self.has_HE:

            if verbosity > 1:
                mout.warningOut("Protonating epsilon site")

            length = 1.00  # Ang

            site = self.NE
            opposite = (self.ND + self.CG) / 2
            vec = site.np_pos - opposite
            vec = length * vec / np.linalg.norm(vec)
            pos = site + vec

            if not delta and self.has_HD:
                atom = self.HD
                atom.name = "HE2"
            else:
                atom = Atom("HE2")
                self.add_atom(atom)
            atom.position = pos

        # remove delta-hydrogen
        if not delta and self.has_HD:
            self.delete_atom("HD1")

        # remove epsilon-hydrogen
        if not epsilon and self.has_HE:
            self.delete_atom("HE2")

        # rename the residue
        if self.has_HD and not self.has_HE:
            self.name = "HSD"
        elif self.has_HE and not self.has_HD:
            self.name = "HSE"
        elif self.has_HD and self.has_HE:
            self.name = "HSP"
