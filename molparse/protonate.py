import mrich as logger


def protonate(
    sys: "System",
    minimise: bool = False,
    pH: float = 7.8,
    remove_residues: list[str] | None = ["DMS", "TRS", "LIG", "CL"],
    trim_terminal_residues: int = 0,
    return_file: bool = False,
) -> "System | str":
    """Protonate a system using pdbfixer and openmm

    :param sys: Input :class:`.System`
    :param minimise: Perform an energy minimisation?
    :param pH: System pH
    :param remove_residues: list of residue names to remove
    :param trim_terminal_residues: number of residues to trim from each end of all protein chains
    :returns: Output :class:`.System`
    """

    from .amber import prep4amber
    from .io import parsePDB
    from openmm.app import PME, ForceField, Modeller, PDBFile, Simulation
    from tempfile import NamedTemporaryFile
    from pdbfixer import PDBFixer

    orig_sys = sys.copy()

    if remove_residues:
        orig_sys.remove_residues(names=remove_residues, no_summary=True)

    if trim_terminal_residues:
        for chain in orig_sys.chains:
            if chain.type != "PRO":
                continue
            residues = [r for r in chain.residues[:trim_terminal_residues]] + [
                r for r in chain.residues[-trim_terminal_residues:]
            ]
            indices = [r.index for r in residues]
            chain.remove_residues(indices=indices)

    # prepare IO files
    pdb_orig = NamedTemporaryFile(mode="w+t", suffix=".pdb")
    pdb_prot = NamedTemporaryFile(mode="w+t", suffix=".pdb")

    # write the original
    orig_sys.write(pdb_orig.name, verbosity=0)

    # fix the PDB termini and hydrogens
    fixer = PDBFixer(pdb_orig.name)
    fixer.findMissingResidues()
    if fixer.missingResidues:
        logger.warning(f"{fixer.missingResidues=}")
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        logger.warning(f"{fixer.nonstandardResidues=}")
    fixer.findMissingAtoms()
    if fixer.missingAtoms:
        logger.warning(f"{fixer.missingAtoms=}")
    if fixer.missingTerminals:
        logger.print(f"{fixer.missingTerminals=}")
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(pH)
    PDBFile.writeFile(fixer.topology, fixer.positions, pdb_prot)
    # PDBFile.writeFile(fixer.topology, fixer.positions, "pdb_prot.pdb")

    if minimise:

        pdb_mini = NamedTemporaryFile(mode="w+t", suffix=".pdb")

        # openmm objects
        pdb = PDBFile(pdb_prot.name)
        forcefield = ForceField("amber99sb.xml", "tip3p.xml")
        modeller = Modeller(pdb.topology, pdb.positions)

        # prepare the simulation
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        integrator = VerletIntegrator(0.001 * picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        # minimise proton positions
        simulation.minimizeEnergy(maxIterations=100)

        # update positions
        positions = simulation.context.getState(getPositions=True).getPositions()

        # write files
        PDBFile.writeFile(simulation.topology, positions, pdb_mini)

        # read in mp.System
        sys = parsePDB(pdb_mini.name, verbosity=0)
        pdb_mini.close()

    else:
        # read in mp.System
        sys = parsePDB(pdb_prot.name, verbosity=0)
        pdb_prot.close()

    # fix residue numbering and chain naming
    for orig_chain, new_chain in zip(orig_sys.chains, sys.chains):

        if orig_chain.sequence != new_chain.sequence:

            assert len(orig_chain.sequence) > len(
                new_chain.sequence
            ), f"Sequences don't match: {orig_chain.sequence} {new_chain.sequence}"

            offset = orig_chain.sequence.find(new_chain.sequence)
            assert offset >= 0

        else:
            offset = 0

        new_chain.name = orig_chain.name

        for i, new_residue in enumerate(new_chain.residues):
            orig_residue = orig_chain.residues[i + offset]
            assert orig_residue.name == new_residue.name
            new_residue.number = orig_residue.number

    # close files
    pdb_orig.close()

    if return_file:
        pdb_out = NamedTemporaryFile(mode="w+t", suffix=".pdb")
        sys.write(pdb_out.name, verbosity=0)
        return sys, pdb_out

    return sys
