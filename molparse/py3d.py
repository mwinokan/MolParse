import py3Dmol


def render(sys, protein="cartoon", ligand="stick", protein_color="spectrum"):

    from tempfile import NamedTemporaryFile

    # protein
    view = py3Dmol.view(data=sys.protein_system.pdb_block)
    view.setStyle({protein: {"color": protein_color}})

    # ligands
    if ligand:
        for res in sys.ligand_residues:
            view.addModel(res.pdb_block, "pdb")
            view.setStyle({"model": -1}, {ligand: {}})

    return view
