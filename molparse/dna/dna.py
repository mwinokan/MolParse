def get_hbond_pairs(residue1, residue2):
    import mout
    import mcol

    legal = check_legality(residue1, residue2)

    if not legal:
        mout.errorOut(
            "Incompatible DNA residues " + str(residue1) + " & " + str(residue2),
            fatal=True,
        )

    if residue1.name.startswith("DA") or residue1.name == "ADE":
        DA_N1 = residue1.get_atom(name="N1")
        DA_N6 = residue1.get_atom(name="N6")
        DA_H61 = residue1.get_atom(name="H61")
        DA_H62 = residue1.get_atom(name="H62")
        DT_O4 = residue2.get_atom(name="O4")
        DT_N3 = residue2.get_atom(name="N3")
        DT_H3 = residue2.get_atom(name="H3")
        return [[DA_N1, DT_N3, DT_H3], [DT_O4, DA_N6, DA_H61, DA_H62]]
    elif residue1.name.startswith("DT") or residue1.name == "THY":
        DA_N1 = residue2.get_atom(name="N1")
        DA_N6 = residue2.get_atom(name="N6")
        DA_H61 = residue2.get_atom(name="H61")
        DA_H62 = residue2.get_atom(name="H62")
        DT_O4 = residue1.get_atom(name="O4")
        DT_N3 = residue1.get_atom(name="N3")
        DT_H3 = residue1.get_atom(name="H3")
        return [[DA_N1, DT_N3, DT_H3], [DT_O4, DA_N6, DA_H61, DA_H62]]
    elif residue1.name.startswith("DC") or residue1.name == "CYT":
        DC_O2 = residue1.get_atom(name="O2")
        DC_N3 = residue1.get_atom(name="N3")
        DC_N4 = residue1.get_atom(name="N4")
        DC_H41 = residue1.get_atom(name="H41")
        DC_H42 = residue1.get_atom(name="H42")
        DG_N1 = residue2.get_atom(name="N1")
        DG_H1 = residue2.get_atom(name="H1")
        DG_N2 = residue2.get_atom(name="N2")
        DG_H21 = residue2.get_atom(name="H21")
        DG_H22 = residue2.get_atom(name="H22")
        DG_C6 = residue2.get_atom(name="C6")
        DG_O6 = residue2.get_atom(name="O6")
        return [
            [DC_O2, DG_N2, DG_H21, DG_H22],
            [DC_N3, DG_N1, DG_H1],
            [DG_O6, DC_N4, DC_H41, DC_H42],
        ]
    elif residue1.name.startswith("DG") or residue1.name == "GUA":
        DC_O2 = residue2.get_atom(name="O2")
        DC_N3 = residue2.get_atom(name="N3")
        DC_N4 = residue2.get_atom(name="N4")
        DC_H41 = residue2.get_atom(name="H41")
        DC_H42 = residue2.get_atom(name="H42")
        DG_N1 = residue1.get_atom(name="N1")
        DG_H1 = residue1.get_atom(name="H1")
        DG_N2 = residue1.get_atom(name="N2")
        DG_H21 = residue1.get_atom(name="H21")
        DG_H22 = residue1.get_atom(name="H22")
        DG_C6 = residue1.get_atom(name="C6")
        DG_O6 = residue1.get_atom(name="O6")
        return [
            [DC_O2, DG_N2, DG_H21, DG_H22],
            [DC_N3, DG_N1, DG_H1],
            [DG_O6, DC_N4, DC_H41, DC_H42],
        ]

    return False


def check_legality(residue1, residue2):
    if residue1.name.startswith("DA") or residue1.name == "ADE":
        if residue2.name.startswith("DT") or residue2.name == "THY":
            return True
    elif residue1.name.startswith("DT") or residue1.name == "THY":
        if residue2.name.startswith("DA") or residue2.name == "ADE":
            return True
    elif residue1.name.startswith("DG") or residue1.name == "GUA":
        if residue2.name.startswith("DC") or residue2.name == "CYT":
            return True
    elif residue1.name.startswith("DC") or residue1.name == "CYT":
        if residue2.name.startswith("DG") or residue2.name == "GUA":
            return True

    return False


def fix_termini(chain):
    chain.residues[0].make_5ter()
    chain.residues[-1].make_3ter()


def make_5ter(residue):
    import mout

    mout.warningOut(
        "amp.dna.make_5ter is deprecated. Use NucleicAcid.make_5ter instead"
    )

    from ..residue import res_type

    print(residue.name, res_type(residue.name))
    assert res_type(residue.name) == "DNA"

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
    try:
        residue.get_atom("P").name = "H5T"
    except AttributeError:
        pass


def make_3ter(residue):
    import mout

    mout.warningOut(
        "amp.dna.make_3ter is deprecated. Use NucleicAcid.make_5ter instead"
    )

    from ..residue import res_type

    assert res_type(residue.name) == "DNA"

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


def prep4gmx(system, verbosity=1):
    system.rename_residues("ADE", "DA", verbosity=verbosity)
    system.rename_residues("THY", "DT", verbosity=verbosity)
    system.rename_residues("CYT", "DC", verbosity=verbosity)
    system.rename_residues("GUA", "DG", verbosity=verbosity)

    # Rename DNA Backbone atoms
    system.rename_atoms("OP1", "O1P", res_filter="D", verbosity=verbosity)
    system.rename_atoms("OP2", "O2P", res_filter="D", verbosity=verbosity)
    system.rename_atoms("C7", "C5M", res_filter="D", verbosity=verbosity)

    system.rename_atoms("H5'1", "H5'", res_filter="D", verbosity=verbosity)
    system.rename_atoms("H5'2", "H5''", res_filter="D", verbosity=verbosity)
    system.rename_atoms("H2'1", "H2'", res_filter="D", verbosity=verbosity)
    system.rename_atoms("H2'2", "H2''", res_filter="D", verbosity=verbosity)

    # Deal with DNA termini
    for chain in system.chains:
        if chain.type == "DNA":
            fix_termini(chain)


def get_dna_basepairs(chain1, chain2):
    # checks
    assert chain1.type == "DNA"
    assert chain2.type == "DNA"
    assert len(chain1) == len(chain2)

    return [[a, b] for a, b in zip(chain1.residues, reversed(chain2.residues))]
