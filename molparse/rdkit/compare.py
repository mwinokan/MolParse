def compare_mols(m1, m2, names=["m1", "m2"], ignore_hydrogen=False):

    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    f1 = CalcMolFormula(m1)
    f2 = CalcMolFormula(m2)

    from molparse.atomtypes import formula_to_atomtype_dict, subtract_atomtype_dict

    d1 = formula_to_atomtype_dict(f1)
    d2 = formula_to_atomtype_dict(f2)

    diff1 = subtract_atomtype_dict(d1, d2, ignore_hydrogen=ignore_hydrogen)
    diff2 = subtract_atomtype_dict(d2, d1, ignore_hydrogen=ignore_hydrogen)

    return diff1, diff2
