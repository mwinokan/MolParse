import os
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem
import numpy as np
from .mol import mol_to_pdb_block
from ..group import AtomGroup

FDEF = AllChem.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)
FEATURE_FAMILIES = FDEF.GetFeatureFamilies()

COMPLEMENTARY_FEATURES = {
    "Donor": ["Acceptor"],  # hydrogen bond
    "Acceptor": ["Donor"],  # hydrogen bond
    "NegIonizable": ["PosIonizable"],  # electrostatic
    "PosIonizable": ["NegIonizable", "Aromatic"],  # electrostatic, pi-cation
    "Aromatic": ["Aromatic", "PosIonizable"],  # pi-stacking, pi-cation
    "Hydrophobe": ["Hydrophobe", "LumpedHydrophobe"],  # hydrophobic
    "LumpedHydrophobe": ["Hydrophobe", "LumpedHydrophobe"],  # hydrophobic
}

INTERACTION_TYPES = {
    ("Hydrophobe", "Hydrophobe"): "Hydrophobic",
    ("LumpedHydrophobe", "Hydrophobe"): "Hydrophobic",
    ("Hydrophobe", "LumpedHydrophobe"): "Hydrophobic",
    ("LumpedHydrophobe", "LumpedHydrophobe"): "Hydrophobic",
    ("Donor", "Acceptor"): "Hydrogen Bond",
    ("Acceptor", "Donor"): "Hydrogen Bond",
    ("NegIonizable", "PosIonizable"): "Electrostatic",
    ("PosIonizable", "NegIonizable"): "Electrostatic",
    ("Aromatic", "Aromatic"): "π-stacking",
    ("Aromatic", "PosIonizable"): "π-cation",
    ("PosIonizable", "Aromatic"): "π-cation",
}


def features_from_mol(mol, protonate=True, group=None):
    raw_features = raw_features_from_mol(mol, protonate=protonate)

    if not group:
        group = AtomGroup.from_pdb_block(mol_to_pdb_block(mol))

    feature_list = []
    for feat in raw_features:

        indices = feat.GetAtomIds()
        family = feat.GetFamily()

        atoms = [group.atoms[i] for i in indices]

        # position from indices
        if len(indices) == 1:
            position = atoms[0].np_pos
        else:
            position = np.mean([a.np_pos for a in atoms], axis=0)

        # f_dict = dict(family=family,position=position,indices=indices,x=position[0],y=position[1],z=position[2])

        feat_obj = Feature(
            family=family,
            atoms=atoms,
            position=position,
            res_name=None,
            res_number=None,
            res_chain=None,
        )

        feature_list.append(feat_obj)

    return feature_list


def raw_features_from_mol(mol, protonate=True):
    # protonate
    if protonate:
        m_pdb_prot = Chem.AddHs(mol)
    else:
        m_pdb_prot = mol

    # solve for a pose
    ps = AllChem.ETKDGv3()
    AllChem.EmbedMolecule(m_pdb_prot, ps)

    # feature factory
    fdef = AllChem.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    )

    # get the features
    raw_features = fdef.GetFeaturesForMol(m_pdb_prot)

    return raw_features


def features_from_group(group, protonate=True):
    m_pdb = group.rdkit_mol
    raw_features = raw_features_from_mol(m_pdb, protonate=protonate)

    feature_list = []

    for feat in raw_features:

        indices = feat.GetAtomIds()
        family = feat.GetFamily()

        # position from indices
        if len(indices) == 1:
            atoms = [group.atoms[indices[0]]]
            position = atoms[0].np_pos
        else:
            atoms = [group.atoms[i] for i in indices]
            position = np.mean([group.atoms[i].np_pos for i in indices], axis=0)

        feat_obj = Feature(
            family=family,
            atoms=atoms,
            position=position,
            res_name=None,
            res_number=None,
            res_chain=None,
        )

        feature_list.append(feat_obj)

    return feature_list


import numpy as np


class Feature(object):
    """Object representing an interaction feature (usually on a protein)"""

    def __init__(
        self,
        family: list,
        atoms: list,
        position: np.ndarray,
        res_name: str,
        res_number: int,
        res_chain: str,
    ):

        self.family = family
        self.atoms = atoms
        self.position = position
        self.res_name = res_name
        self.res_number = res_number
        self.res_chain = res_chain

        self.atom_numbers = [a.number for a in self.atoms]

    @property
    def x(self):
        """x coordinate"""
        return self.position[0]

    @property
    def y(self):
        """y coordinate"""
        return self.position[1]

    @property
    def z(self):
        """z coordinate"""
        return self.position[2]

    @property
    def dict(self):
        """Dictionary containing family and coordinates"""
        return dict(
            family=self.family,
            x=self.x,
            y=self.y,
            z=self.z,
        )

    def __repr__(self):
        return f"{self.family} @ {self.x:.2f} {self.y:.2f} {self.z:.2f}"

    @property
    def name_number_chain_str(self):
        """String containing the residue name, residue number, and chain"""
        return f"{self.res_name} {self.res_number} {self.res_chain}"

    @property
    def family_name_number_chain_str(self):
        """String containing the feature family, residue name, residue number, and chain"""
        return f"{self.family} {self.res_name} {self.res_number} {self.res_chain}"

    @property
    def family_name_number_chain_atoms_str(self):
        """String containing the feature family, residue name, residue number, chain, and atoms"""
        return f"{self.family} {self.res_name} {self.res_number} {self.res_chain} {self.atoms}"

    def __sub__(self, other):
        if isinstance(other, np.ndarray):
            return self.position - other
        else:
            return self.position - other.position
