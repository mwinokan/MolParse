def view(atoms, **kwargs):
    from .system import System

    if isinstance(atoms, System):
        names = [a.name for a in atoms.atoms]
        atoms = atoms.ase_atoms
        for n, a in zip(names, atoms):
            if len(n) > 1:
                try:
                    a.tag = int(n[1:])
                except ValueError:
                    continue
    from ase import visualize

    return visualize.view(atoms, **kwargs)
