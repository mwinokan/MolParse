
import os
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem
import numpy as np


FDEF = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'))
FEATURE_FAMILIES = FDEF.GetFeatureFamilies()

COMPLEMENTARY_FEATURES = {
    "Donor": "Acceptor",
    "Acceptor": "Donor",
    "NegIonizable": "PosIonizable",
    "PosIonizable": "NegIonizable",
    "Aromatic": "Aromatic",
    "Aromatic": "PosIonizable",
    "PosIonizable": "Aromatic",
    "Hydrophobe": "Hydrophobe", 
}

def features_from_group(group):
    m_pdb = group.rdkit_mol

    # protonate
    m_pdb_prot = Chem.AddHs(m_pdb)

    # solve for a pose
    ps = AllChem.ETKDGv3()
    AllChem.EmbedMolecule(m_pdb_prot,ps)

    # feature factory
    fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'))

    # get the features
    raw_features = fdef.GetFeaturesForMol(m_pdb_prot)

    feature_list = []

    for feat in raw_features:

        indices = feat.GetAtomIds()
        family = feat.GetFamily()

        # position from indices
        if len(indices) == 1:
            position = group.atoms[indices[0]].np_pos
        else:
            position = np.mean([group.atoms[i].np_pos for i in indices],axis=0)

        f_dict = dict(family=family,position=position,indices=indices,x=position[0],y=position[1],z=position[2])

        feature_list.append(f_dict)

    return feature_list

import numpy as np

class Feature(object):

    def __init__(self, family: list, atoms: list, position: np.ndarray, sidechain: bool, res_name: str, res_number: int, res_chain: str):
            
        self.family = family
        self.atoms = atoms
        self.position = position
        self.sidechain = sidechain
        self.res_name = res_name
        self.res_number = res_number
        self.res_chain = res_chain

        self.atom_numbers = [a.number for a in self.atoms]

    @property
    def x(self):
        return self.position[0]

    @property
    def y(self):
        return self.position[1]

    @property
    def z(self):
        return self.position[2]
    
    @property
    def dict(self):
        return dict(
            family=self.family,
            x=self.x,
            y=self.y,
            z=self.z,
        )
    
    def __repr__(self):
        return str(list(self.dict.values()))

    @property
    def name_number_chain_str(self):
        return f'{self.res_name} {self.res_number} {self.res_chain}'

    @property
    def family_name_number_chain_str(self):
        return f'{self.family} {self.res_name} {self.res_number} {self.res_chain}'

    def __sub__(self,other):
        assert isinstance(other,Feature)
        return self.position - other.position
