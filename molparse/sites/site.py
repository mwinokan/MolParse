# Site = namedtuple('Site',['types','atoms','position','sidechain','res_name','res_number'])

import numpy as np
import mout
import mcol
from ..residue import Residue
from ..monte import CompoundVolume, CappedCone
import mgo


class Site:
    """Protein Interaction Site"""

    def __init__(
        self,
        interaction_types: list,
        atoms: list,
        position: np.ndarray,
        sidechain: bool,
        res_name: str,
        res_number: int,
        res_chain: str,
    ):

        self.types = interaction_types
        self.atoms = atoms
        self.position = position
        self.sidechain = sidechain
        self.res_name = res_name
        self.res_number = res_number
        self.res_chain = res_chain

        self.atom_numbers = [a.number for a in self.atoms]

        self.accessible = None
        self.in_screen = None

    @property
    def type_str(self):
        return ", ".join(self.types)

    @property
    def atom_names(self):
        return [a.name for a in self.atoms]

    @property
    def atom_str(self):
        return ", ".join(self.atom_names)

    @property
    def name(self):
        return f"{self.res_name} {self.res_number} {self.res_chain} {self.type_str} {self.atom_str}"

    def __repr__(self):
        return self.name

    def summary(self):
        mout.out(f"{mcol.varName}{self.name}{mcol.clear}: {self.position}")

    def accessibility(
        self,
        residue: Residue,
        candidates: list,
        cutoff: float = 5.0,
        proximity: float = 1.5,
        samples: int = 4000,
        skip_hydrogen: bool = True,
        fig=None,
    ):
        """

        residue: parent mp.Residue object of the Site
        candidates: list of mp.Residue objects to consider for occlusion
        cutoff: the radius within which to look for valid binding pocket points
        proximity: minimum radius of binding pocket point
        samples: number of monte carlo samples for volume calculation
        skip_hydrogen: do not calculate occlusion by hydrogen atoms

        """

        self.nearby = residue.get_nearby(candidates, cutoff=cutoff)

        self.forbidden_volume = CompoundVolume(self.position)

        for res in [residue] + self.nearby:
            for atom in res.atoms:
                if skip_hydrogen and atom.species == "H":
                    continue

                d = np.linalg.norm(atom.np_pos - self.position)

                if d < 0.5:
                    continue

                if d > atom.vdw_radius + cutoff:
                    continue

                cone = CappedCone(
                    self.position,
                    atom.np_pos,
                    atom.vdw_radius,
                    name=f"{res.name_number_str} {atom.name_number_str}",
                )

                self.forbidden_volume.add_volume(cone)

        from ..monte import mc_spherical

        if fig is not None:

            fraction, valid_volume, points_in, points_out = mc_spherical(
                self.forbidden_volume,
                cutoff,
                proximity,
                self.position,
                samples,
                boolean=False,
                inverse=True,
                points=True,
            )

            import plotly.graph_objects as go

            trace = mgo.point_trace(points_in, "valid")
            fig.add_trace(trace)
            trace = mgo.point_trace(points_out, "invalid")
            fig.add_trace(trace)

            self.accessible = bool(points_in)

        else:

            self.accessible = mc_spherical(
                self.forbidden_volume,
                cutoff,
                proximity,
                self.position,
                samples,
                boolean=True,
                inverse=True,
            )

        return self.accessible


# add __eq__ to compare to PLIP interaction

# from collections import UserDict

# class SiteList(UserDict):
