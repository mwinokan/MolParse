from .residue import Residue, res_type
from .group import AtomGroup


def alphabet(name):
    if name == "DA":
        return "A"
    elif name == "DT":
        return "T"
    elif name == "DG":
        return "G"
    elif name == "DC":
        return "C"
    elif name == "DA3":
        return "A"
    elif name == "DT3":
        return "T"
    elif name == "DG3":
        return "G"
    elif name == "DC3":
        return "C"
    elif name == "DA5":
        return "A"
    elif name == "DT5":
        return "T"
    elif name == "DG5":
        return "G"
    elif name == "DC5":
        return "C"
    else:
        import mout

        mout.errorOut(f"Unsupported nucleic acid residue name: {name}")
        return None


def longname(name):
    if name == "DA":
        return "Adenine"
    elif name == "DT":
        return "Thymine"
    elif name == "DG":
        return "Guanine"
    elif name == "DC":
        return "Cytosine"
    elif name == "DA3":
        return "Adenine 3-Terminus"
    elif name == "DT3":
        return "Thymine 3-Terminus"
    elif name == "DG3":
        return "Guanine 3-Terminus"
    elif name == "DC3":
        return "Cytosine 3-Terminus"
    elif name == "DA5":
        return "Adenine 5-Terminus"
    elif name == "DT5":
        return "Thymine 5-Terminus"
    elif name == "DG5":
        return "Guanine 5-Terminus"
    elif name == "DC5":
        return "Cytosine 5-Terminus"
    else:
        import mout

        mout.errorOut(f"Unsupported nucleic acid residue name: {name}")
        return None


class NucleicAcid(Residue):
    """Class for Nucleic Acid Residue

    These objects should not be created by the user,
    but constructed automatically when parsing a
    coordinate file via amp.parsePDB or otherwise"""

    def __init__(
        self,
        name: str,
        index: int = None,
        number: int = None,
        chain: str = None,
        atoms=None,
    ):
        assert res_type(name) == "DNA"
        super(NucleicAcid, self).__init__(name, index, number, chain, atoms)

        self._backbone = None
        self._nucleobase = None

    @property
    def is_purine(self):
        return self.name in ["DG", "DA"]

    @property
    def is_pyrimidine(self):
        return self.name in ["DC", "DT"]

    @property
    def letter(self):
        """Nucleic acid alphabet"""
        return alphabet(self.name)

    @property
    def longname(self):
        """Nucleic acid chemical name"""
        return longname(self.name)

    @property
    def type(self):
        """Fixed typing"""
        return "DNA"

    @property
    def backbone(self):
        """Get backbone atoms"""
        if not self._backbone:
            self._backbone = self.get_atom(self.bb_names)
        return self._backbone

    def remove_backbone(self, add_link=True):
        """Remove backbone atoms"""

        if add_link:
            self.get_atom("C1'").set_name("HLNK")

        for atom in self.backbone:
            if not atom:
                continue
            self.delete_atom(atom.name)

    @property
    def nucleobase(self):
        """Get nucleobase atoms"""
        if not self._nucleobase:
            self._nucleobase = self.get_atom(self.nonbb_names)
        return self._nucleobase

    @property
    def bb_names(self):
        return [
            "H5T",
            "O5'",
            "C5'",
            "H5'",
            "H5''",
            "C4'",
            "H4'",
            "O4'",
            "C1'",
            "H1'",
            "C2'",
            "H2'",
            "H2''",
            "C3'",
            "H3'",
            "O3'",
            "P",
            "O1P",
            "O2P",
            "H3T",
        ]

    @property
    def nonbb_names(self):
        """Get nucleobase atoms"""
        return [n for n in self.atom_names() if n not in self.bb_names]

    def make_5ter(residue):
        """Create a 5-Terminus"""

        # Append 5 to residue name
        if not residue.name.endswith("5"):
            residue.name += "5"

        # Remove HTER/H5T
        residue.delete_atom("HTER")
        residue.delete_atom("H5T")
        residue.delete_atom("HO5'")
        residue.delete_atom("OXT")
        residue.delete_atom("O5T")
        residue.delete_atom("O1P")
        residue.delete_atom("O2P")
        residue.delete_atom("OP1")
        residue.delete_atom("OP2")

        # Rename P->H5T
        atom = residue.get_atom("P")
        if atom is not None:
            atom.name = "H5T"

    def make_3ter(residue):
        """Create a 3-Terminus"""

        # Append 3 to residue name
        if not residue.name.endswith("3"):
            residue.name += "3"

        # Remove HTER/H3T
        residue.delete_atom("O1P3")
        residue.delete_atom("O2P3")
        residue.delete_atom("O3T")
        residue.delete_atom("H3T")
        residue.delete_atom("HO3'")
        residue.delete_atom("HCAP")

        # Rename P3->H3T
        atom = residue.get_atom("P3")
        if atom is not None:
            atom.name = "H3T"

    def flip(self):

        import mout
        import mcol

        mout.header(
            f"Flipping {mcol.arg}{self.longname} {mcol.clear}{mcol.bold}({mcol.arg}{self.name_number_str}{mcol.clear}{mcol.bold})"
        )

        nucleobase_group = AtomGroup.from_any(
            f"{self.name}.nucleobase", self.nucleobase
        )

        if self.is_purine:
            source = nucleobase_group.atoms["N9"][0]
            target = nucleobase_group.atoms["O6"][0]

        vector = target - source

        nucleobase_group.rotate(180, vector, source.np_pos)

        for a1, a2 in zip(self.nucleobase, nucleobase_group.atoms):
            a1.position = a2.position

    def mutate(self, newname, show=False):
        """Mutate this nucleic acid to another"""

        import mout
        import mcol

        mout.headerOut(
            f"Mutating {mcol.arg}{self.longname}{mcol.clear + mcol.bold} --> {mcol.arg}{longname(newname)}{mcol.clear + mcol.bold} ({mcol.result}{self.letter}{self.number}{alphabet(newname)}{mcol.clear + mcol.bold})"
        )

        if self.name == newname:
            mout.warningOut("Skipping mutation with self == target!")
            return

        import os

        amp_path = os.path.dirname(__file__)

        from .io import parse

        if newname in ["DA", "DT"]:
            ref_sys = parse(f"{amp_path}/ref/AT_FLAT.pdb", verbosity=0)
        elif newname in ["DG", "DC"]:
            ref_sys = parse(f"{amp_path}/ref/GC_FLAT.pdb", verbosity=0)
        else:
            mout.errorOut("Unsupported target mutation")

        rem = [c for c in ref_sys.chain_names if c not in newname[-1]][0]
        ref_sys.remove_chain(rem, verbosity=0)
        ref_sys.fix_indices()
        ref_res = ref_sys.residues[0]

        if self.is_purine and ref_res.is_purine:

            # get atoms to align by
            names = ["N9", "C6", "C2"]
            start = ref_res.get_atom(names)
            target = self.get_atom(names)

            # align the system
            extra = ref_sys.align_by_posmap([start, target])

            # delete intersecting atoms
            ref_res.delete_atom("H9", verbosity=0)

        elif self.is_pyrimidine and ref_res.is_pyrimidine:

            # get atoms to align by
            names = ["C6", "C4", "C2"]
            start = ref_res.get_atom(names)
            target = self.get_atom(names)

            # align the system
            extra = ref_sys.align_by_posmap([start, target])

            # shift to connect to the backbone
            ref_res.translate(self.get_atom("N1") - ref_res.get_atom("N1"))

            # delete intersecting atoms
            ref_res.delete_atom("H1", verbosity=0)

        elif self.is_purine and ref_res.is_pyrimidine:
            mout.warningOut("Purine --> Pyrimidine")

            # get atoms to align by
            start = ref_res.get_atom(["C6", "C4", "C2"])
            target = self.get_atom(["C4", "C6", "C2"])

            # align the system
            extra = ref_sys.align_by_posmap([start, target])

            # shift to connect to the backbone
            ref_res.translate(self.get_atom("N9") - ref_res.get_atom("N1"))

            # delete intersecting atoms
            ref_res.delete_atom("H1", verbosity=0)

        elif self.is_pyrimidine and ref_res.is_purine:
            mout.warningOut("Pyrimidine --> Purine")

            # get atoms to align by
            start = ref_res.get_atom(["C4", "C6", "C2"])
            target = self.get_atom(["C6", "C4", "C2"])

            # align the system
            extra = ref_sys.align_by_posmap([start, target])

            # shift to connect to the backbone
            ref_res.translate(self.get_atom("N1") - ref_res.get_atom("N9"))

            # delete intersecting atoms
            ref_res.delete_atom("H9", verbosity=0)

        self.name = newname

        if show:
            view_sys = ref_sys.copy()
            for atom in self.atoms:
                view_sys.add_atom(atom)
            view_sys.plot3d(extra, 1.0)

        # delete intersecting atoms
        self.delete_atom(self.nonbb_names, verbosity=0)

        # add the reference nucleobase
        for atom in ref_res.atoms:
            atom._NUMBER = None
            self.add_atom(atom)

        if self.letter == "G":
            from .guanine import Guanine

            self = Guanine.from_nucleic(self)

        return self

    def oxidise(self):

        assert self.name == "DG"

        import numpy as np
        from .atom import Atom

        import mout
        import mcol

        mout.header(
            f"Oxidising {mcol.arg}{self.longname}{mcol.header} --> {mcol.result}{mcol.bold}8-Oxo-7,8-dihydroguanine (DOG)"
        )

        # H8 --> O8
        self.get_atom("H8").name = "O8"

        # +H7
        N7 = self.get_atom("N7")
        N9 = self.get_atom("N9")
        C4 = self.get_atom("C4")
        vec = N7 - (N9 + C4) / 2
        vec /= np.linalg.norm(vec)
        atom = Atom("H7")
        atom.position = N7 + vec
        self.add_atom(atom)

        # +H2
        N2 = self.get_atom("N2")
        C2 = self.get_atom("C2")
        N3 = self.get_atom("N3")
        vec = np.cross(N2 - C2, N3 - C2)
        vec /= np.linalg.norm(vec)
        atom = Atom("H2")
        atom.position = C2 + vec
        self.add_atom(atom)

        # +H3
        C6 = self.get_atom("C6")
        vec = N3 - C6
        vec /= np.linalg.norm(vec)
        atom = Atom("H3")
        atom.position = N3 + vec
        self.add_atom(atom)

        self.name = "DOG"
