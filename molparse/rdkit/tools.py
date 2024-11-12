from mlog import setup_logger

logger = setup_logger("MolParse")


def sdf_align(
    sdf_path: str,
    source_protein_path: str,
    target_protein_path: str,
    mol_col: str = "ROMol",
    name_col: str = "ID",
    backbone_only: bool = True,
    heavy_atoms_only: bool = True,
):
    """Align molecules in an SDF

    :param sdf_path: input file
    :param source_protein_path: source protein conformation
    :param target_protein_path: target protein conformation
    :param mol_col: name of molecule object column
    :param name_col: name of molecule ID column
    :param backbone_only: align using protein backbone atoms only
    :param heavy_atoms_only: align using protein heavy atoms only

    """

    from ..io import parse
    from numpy import array, double
    from rdkit.Chem import PandasTools, rdMolTransforms, Mol
    from .mol import mol_from_smiles, mol_to_smiles, mol_to_pdb_block

    # parse the reference files
    source = parse(source_protein_path)
    target = parse(target_protein_path)

    # align the proteins
    source_center, target_center, rotation = source.align_to(
        target,
        protein_only=True,
        return_transformations=True,
        backbone_only=backbone_only,
        heavy_atoms_only=heavy_atoms_only,
    )

    # create the transformation matrices

    translation_1 = array(
        [
            [1, 0, 0, -source_center[0]],
            [0, 1, 0, -source_center[1]],
            [0, 0, 1, -source_center[2]],
            [0, 0, 0, 1],
        ],
        dtype=double,
    )

    translation_2 = array(
        [
            [1, 0, 0, target_center[0]],
            [0, 1, 0, target_center[1]],
            [0, 0, 1, target_center[2]],
            [0, 0, 0, 1],
        ],
        dtype=double,
    )

    rotation = array(
        [
            [*rotation[0], 0],
            [*rotation[1], 0],
            [*rotation[2], 0],
            [0, 0, 0, 1],
        ],
        dtype=double,
    )

    df = PandasTools.LoadSDF(sdf_path)

    for i, row in df.iterrows():
        if row[name_col] == "ver_1.2":
            logger.warning("Skipping Fragalysis header molecule")
            continue

        mol = row[mol_col]

        # translate to origin
        rdMolTransforms.TransformConformer(mol.GetConformer(0), translation_1)

        # apply rotation
        rdMolTransforms.TransformConformer(mol.GetConformer(0), rotation)

        # translate to target
        rdMolTransforms.TransformConformer(mol.GetConformer(0), translation_2)

    out_path = sdf_path.replace(".sdf", "_aligned.sdf")

    logger.writing(out_path)
    PandasTools.WriteSDF(df, out_path, mol_col, name_col, list(df.columns))


def sdf_combine(files: list, output: str, mol_col: str = "ROMol", name_col: str = "ID"):
    """Concatenate SDF files

    :param files: generator or iterable of file names/paths
    :param output: file name/path of output SDF
    :param mol_col: name of molecule object column
    :param name_col: name of molecule ID column

    """

    from rdkit.Chem import PandasTools
    from pathlib import Path
    from pandas import concat
    import mrich

    out_path = Path(output)

    dfs = []
    for file in mrich.track(files, prefix="parsing SDFs"):
        in_path = Path(file)

        df = PandasTools.LoadSDF(in_path)

        dfs.append(df)

    df = concat(dfs, ignore_index=True)

    mrich.writing(out_path)
    PandasTools.WriteSDF(
        df, out_path, molColName=mol_col, idName=name_col, properties=df.columns
    )
