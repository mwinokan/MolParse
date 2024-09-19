from rdkit.Chem import Fragments, MolFromSmarts, Draw

import os
from rdkit import RDConfig

import logging

logger = logging.getLogger("MolParse")


def classify_mol(mol, draw=True):
    """Find RDKit Fragments within the molecule and draw them (or just return a list of (descriptor, count) tuples)"""

    mols = []
    highlights = []
    counts = []

    for name, descriptor in FRAGMENT_DESCRIPTORS.items():

        f = getattr(Fragments, name)

        n = f(mol)

        if not n:
            continue

        smarts = FRAGMENT_SMARTS[name]

        pattern = MolFromSmarts(smarts)

        matches = mol.GetSubstructMatches(pattern, uniquify=True)

        mols.append(mol)
        counts.append((descriptor, n))

        highlight = []

        for match in matches:
            for index in match:
                highlight.append(index)

        highlights.append(highlight)

    if draw:
        legends = [f"{n} x {descriptor}" for descriptor, n in counts]
        drawing = Draw.MolsToGridImage(
            mols, highlightAtomLists=highlights, legends=legends
        )
        display(drawing)
    else:
        return counts


defaultPatternFileName = os.path.join(RDConfig.RDDataDir, "FragmentDescriptors.csv")

FRAGMENT_SMARTS = {}

for line in open(defaultPatternFileName).readlines():

    if not len(line) or line[0] == "#":
        continue

    split = line.split("\t")

    if len(split) < 3:
        continue

    key = split[0].replace("=", "_").replace("-", "_")
    smarts = split[2]

    FRAGMENT_SMARTS[key] = smarts

FRAGMENT_DESCRIPTORS = {
    "fr_Al_COO": "aliphatic carboxylic acid",
    "fr_Al_OH": "aliphatic hydroxyl",
    "fr_Al_OH_noTert": "aliphatic hydroxyl (not tertiary)",
    "fr_ArN": "nitrogen group on aromatic",
    "fr_Ar_COO": "aromatic carboxylic acid",
    "fr_Ar_N": "aromatic nitrogen",
    "fr_Ar_NH": "aromatic amine",
    "fr_Ar_OH": "aromatic hydroxyl",
    "fr_COO": "carboxylic acid",
    "fr_COO2": "carboxylic acid",
    "fr_C_O": "carbonyl",
    "fr_C_O_noCOO": "carbonyl (not carboxylic acid)",
    "fr_C_S": "thiocarbonyl",
    "fr_HOCCN": "C(OH)CCN-Ctert-alkyl or C(OH)CCNcyclic",
    "fr_Imine": "imine",
    "fr_NH0": "tertiary amine",
    "fr_NH1": "secondary amine",
    "fr_NH2": "primary amine",
    "fr_N_O": "hydroxylamine",
    "fr_Ndealkylation1": "XCCNR",
    "fr_Ndealkylation2": "tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)",
    "fr_Nhpyrrole": "H-pyrrole",
    "fr_SH": "thiol",
    "fr_aldehyde": "aldehyde",
    "fr_alkyl_carbamate": "alkyl carbamate (subject to hydrolysis)",
    "fr_alkyl_halide": "alkyl halide",
    "fr_allylic_oxid": "allylic oxidation sites excluding steroid dienone",
    "fr_amide": "amide",
    "fr_amidine": "amidine",
    "fr_aniline": "aniline",
    "fr_aryl_methyl": "aryl methyl sites for hydroxylation",
    "fr_azide": "azide",
    "fr_azo": "azo",
    "fr_barbitur": "barbiturate",
    "fr_benzene": "benzene",
    "fr_benzodiazepine": "benzodiazepines with no additional fused rings",
    "fr_bicyclic": "bicyclic",
    "fr_diazo": "diazo",
    "fr_dihydropyridine": "dihydropyridine",
    "fr_epoxide": "epoxide",
    "fr_ester": "ester",
    "fr_ether": "ether oxygen (including phenoxy)",
    "fr_furan": "furan",
    "fr_guanido": "guanidine",
    "fr_halogen": "halogen",
    "fr_hdrzine": "hydrazine",
    "fr_hdrzone": "hydrazone",
    "fr_imidazole": "imidazole",
    "fr_imide": "imide",
    "fr_isocyan": "isocyanate",
    "fr_isothiocyan": "isothiocyanate",
    "fr_ketone": "ketone",
    "fr_ketone_Topliss": "ketone excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha",
    "fr_lactam": "beta-lactam",
    "fr_lactone": "cyclic ester (lactone)",
    "fr_methoxy": "methoxy",
    "fr_morpholine": "morpholine",
    "fr_nitrile": "nitrile",
    "fr_nitro": "nitro",
    "fr_nitro_arom": "nitro benzene ring substituent",
    "fr_nitro_arom_nonortho": "non-ortho nitro benzene ring substituent",
    "fr_nitroso": "nitroso group, excluding NO2",
    "fr_oxazole": "oxazole",
    "fr_oxime": "oxime",
    "fr_para_hydroxylation": "para-hydroxylation site",
    "fr_phenol": "phenol",
    "fr_phenol_noOrthoHbond": "phenolic OH excluding ortho intramolecular Hbond substituent",
    "fr_phos_acid": "phosphoric acid",
    "fr_phos_ester": "phosphoric ester",
    "fr_piperdine": "piperdine ring",
    "fr_piperzine": "piperzine ring",
    "fr_priamide": "primary amide",
    "fr_prisulfonamd": "primary sulfonamide",
    "fr_pyridine": "pyridine ring",
    "fr_quatN": "quaternary nitrogen",
    "fr_sulfide": "thioether",
    "fr_sulfonamd": "sulfonamide",
    "fr_sulfone": "sulfone",
    "fr_term_acetylene": "terminal acetylene",
    "fr_tetrazole": "tetrazole ring",
    "fr_thiazole": "thiazole ring",
    "fr_thiocyan": "thiocyanate",
    "fr_thiophene": "thiophene ring",
    "fr_unbrch_alkane": "unbranched alkane of at least 4 members (excludes halogenated alkanes)",
    "fr_urea": "urea",
}
