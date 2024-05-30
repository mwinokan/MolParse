
Aligning molecules
==================

It is quite common to need to align two molecular systems to each other. The :class:`.System` class uses `ase.build.rotate` to minimise rotation and translation between groups of atoms.

In general if you have two :class:`.System` objects, you can align one to the other using:

::

	system1.align_to(system2)

N.B. This will fail or at the very least give unexpected results when the number and type of atoms in the two systems are different.

Aligning by protein sub-System
------------------------------

In cases where you want to align the systems by considering only their proteins use the `protein_only` argument. The `backbone_only` and `heavy_atom_only` options can further be used to restrict the aligning calculation.

::

	system1.align_to(system2, protein_only=True, backbone_only=True)

This code was intended to align single similar protein chains onto one another, i.e. mapping similar sequences onto each other. Differences between the protein sequences will be resolved by keeping only common residues.

Aligning proteins with multiple chains has not been tested, it is recommended to separately align each chain to a target.
 