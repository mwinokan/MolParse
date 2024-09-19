from rdkit import Chem


def strip_quotes(val):
    if val[0] == '"':
        val = val[1:]
    if val[-1] == '"':
        val = val[:-1]
    return val


BOND_TYPES = {
    "single": Chem.rdchem.BondType.SINGLE,
    "double": Chem.rdchem.BondType.DOUBLE,
    "triple": Chem.rdchem.BondType.TRIPLE,
    "SINGLE": Chem.rdchem.BondType.SINGLE,
    "DOUBLE": Chem.rdchem.BondType.DOUBLE,
    "TRIPLE": Chem.rdchem.BondType.TRIPLE,
    "aromatic": Chem.rdchem.BondType.AROMATIC,
    "deloc": Chem.rdchem.BondType.SINGLE,
    "SING": Chem.rdchem.BondType.SINGLE,
    "DOUB": Chem.rdchem.BondType.DOUBLE,
    "TRIP": Chem.rdchem.BondType.TRIPLE,
}
