def printEnergy(atoms, perAtom=True, precision=8, printScript=False):
    import mcol  # https://github.com/mwinokan/MPyTools
    import mout  # https://github.com/mwinokan/MPyTools

    if perAtom:
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
    else:
        epot = atoms.get_potential_energy()
        ekin = atoms.get_kinetic_energy()

    mout.varOut(
        "E_pot",
        epot,
        unit="eV/atom",
        valCol=mcol.result,
        precision=precision,
        printScript=printScript,
        end=", ",
    )
    mout.varOut(
        "E_kin",
        ekin,
        unit="eV/atom",
        valCol=mcol.result,
        precision=precision,
        printScript=False,
    )
