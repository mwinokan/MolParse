from mlog import setup_logger

logger = setup_logger("MolParse")


def protonate(
    sys: "System",
    minimise: bool = False,
    pH: float = 7.8,
    remove_residues: list[str] | None = ["DMS", "TRS", "LIG", "CL"],
    trim_terminal_residues: int = 0,
) -> "System":
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

    sys = sys.copy()

    if remove_residues:
        sys.remove_residues(names=remove_residues)

    if trim_terminal_residues:
        for chain in sys.chains:
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
    pdb_mini = NamedTemporaryFile(mode="w+t", suffix=".pdb")

    # write the original
    sys.write(pdb_orig.name, verbosity=0)

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
        logger.info(f"{fixer.missingTerminals=}")
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(pH)
    PDBFile.writeFile(fixer.topology, fixer.positions, pdb_prot)
    PDBFile.writeFile(fixer.topology, fixer.positions, "pdb_prot.pdb")

    if minimise:

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
        sys = parsePDB(pdb_mini.name)

    else:
        # read in mp.System
        sys = parsePDB(pdb_prot.name)

    # close files
    pdb_orig.close()
    pdb_prot.close()
    pdb_mini.close()

    return sys
