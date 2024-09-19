import mout
from rdkit import Chem, Geometry


# class Molecule(Chem.rdchem.Mol):
#   def __init__(self,*args,**kwargs):
#       super(ClassName, self).__init__(*args,**kwargs)


def mol_from_pdb_block(pdb_block):
    return Chem.rdmolfiles.MolFromPDBBlock(pdb_block)


def mol_to_smiles(mol):
    return Chem.MolToSmiles(mol)


def mol_from_smiles(smiles):
    return Chem.MolFromSmiles(smiles)


def mol_to_pdb_block(mol, hetatm_only: bool = False):
    if hetatm_only:
        buffer = Chem.MolToPDBBlock(mol)

        new = []
        for line in buffer.split("\n"):
            if line.startswith("HETATM"):
                new.append(line)
        return "\n".join(new)

    else:
        return Chem.MolToPDBBlock(mol)


def mol_to_AtomGroup(mol):
    from ..group import AtomGroup

    group = AtomGroup.from_pdb_block(mol_to_pdb_block(mol))
    if hasattr(mol, "_Name"):
        group.name = mol._Name
    return group


def mol_from_AtomGroup(group):
    return group.mol


def protonate(mol, embed=True, align=True, verbosity=1):
    mol_prot = Chem.AddHs(mol)
    if embed and align:
        if verbosity:
            mout.warning("May lose exact pose", code="mp.rdkit.protonate.1")
        mol = Chem.AllChem.ConstrainedEmbed(mol_prot, mol)
    elif embed:
        if verbosity:
            mout.warning(
                "Disregarding original coordinates", code="mp.rdkit.protonate.2"
            )
        ps = Chem.AllChem.ETKDGv3()
        Chem.AllChem.EmbedMolecule(mol_prot, ps)
    return mol_prot


def copy_mol(mol):
    return Chem.Mol(mol)


def mol_from_cif(cif_file):

    from gemmi import cif
    from .xca_utils import strip_quotes, BOND_TYPES

    mol = Chem.RWMol()
    conf = Chem.Conformer()

    doc = cif.read(str(cif_file))

    # Diamond CIFs have two blocks, but the one we want will be named data_comp_LIG
    block = doc.find_block("comp_LIG")

    # Other CIFs have unpredictable block names, so let's hope there is only one
    if not block:
        block = doc.sole_block()
    if not block:
        print("sole block not found")
        return None

    comp_ids = block.find_loop("_chem_comp_atom.comp_id")
    atom_ids = block.find_loop("_chem_comp_atom.atom_id")
    atom_symbols = block.find_loop("_chem_comp_atom.type_symbol")
    # coordinates are sometimes called "x" and sometimes "model_Cartn_x" etc.
    x = block.find_loop("_chem_comp_atom.x")
    if not x:
        x = block.find_loop("_chem_comp_atom.model_Cartn_x")
    y = block.find_loop("_chem_comp_atom.y")
    if not y:
        y = block.find_loop("_chem_comp_atom.model_Cartn_y")
    z = block.find_loop("_chem_comp_atom.z")
    if not z:
        z = block.find_loop("_chem_comp_atom.model_Cartn_z")
    charges = [0] * len(atom_ids)
    if block.find_loop("_chem_comp_atom.charge"):
        charges = list(block.find_loop("_chem_comp_atom.charge"))
    elif block.find_loop("_chem_comp_atom.partial_charge"):
        charges = list(block.find_loop("_chem_comp_atom.partial_charge"))

    atoms = {}
    ligand_name = None
    for name, s, id, px, py, pz, charge in zip(
        comp_ids, atom_symbols, atom_ids, x, y, z, charges
    ):
        # sometimes that atom ids are wrapped in double quotes
        if ligand_name is None:
            ligand_name = name
        elif name != ligand_name:
            print(
                "WARNING: ligand name has changed from {} to {}. Old name will be used.".format(
                    ligand_name, name
                )
            )

        id = strip_quotes(id)

        if len(s) == 2:
            s = s[0] + s[1].lower()

        atom = Chem.Atom(s)
        atom.SetFormalCharge(round(float(charge)))
        atom.SetProp("atom_id", id)
        idx = mol.AddAtom(atom)
        atom.SetIntProp("idx", idx)
        atoms[id] = atom

        point = Geometry.Point3D(float(px), float(py), float(pz))
        conf.SetAtomPosition(idx, point)

    atom1 = block.find_loop("_chem_comp_bond.atom_id_1")
    atom2 = block.find_loop("_chem_comp_bond.atom_id_2")
    bond_type = block.find_loop("_chem_comp_bond.type")
    if not bond_type:
        bond_type = block.find_loop("_chem_comp_bond.value_order")

    for a1, a2, bt in zip(atom1, atom2, bond_type):
        mol.AddBond(
            atoms[strip_quotes(a1)].GetIntProp("idx"),
            atoms[strip_quotes(a2)].GetIntProp("idx"),
            BOND_TYPES[bt],
        )

    Chem.SanitizeMol(mol)
    mol.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(mol)
    mol = Chem.RemoveAllHs(mol)

    mol.SetProp("_Name", ligand_name)

    return mol
