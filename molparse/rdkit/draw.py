
import py3Dmol
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFMCS, AllChem, Draw
IPythonConsole.ipython_3d = True

def draw_mol(m, feats=None, p=None, confId=-1, hydrogen=True):
    if p is None:
        p = py3Dmol.view(width=400, height=400)
    p.removeAllModels()
    if not hydrogen:
        m = Chem.RemoveHs(m)
    IPythonConsole.addMolToView(m,p,confId=confId)
    feats = feats or []
    for feat in feats:
        pos = feat.GetPos()
        clr = featColors.get(feat.GetFamily(),(.5,.5,.5))
        p.addSphere({'center':{'x':pos.x,'y':pos.y,'z':pos.z},'radius':.5,'color':colorToHex(clr)});
    p.zoomTo()
    return p.show()

def draw_mols(ms, p=None, confId=-1, hydrogen=True,colors=('cyanCarbon','redCarbon','blueCarbon')):
    if p is None:
        p = py3Dmol.view(width=400, height=400)
    p.removeAllModels()
    for i,m in enumerate(ms):
        if not hydrogen:
            m = Chem.RemoveHs(m)
        IPythonConsole.addMolToView(m,p,confId=confId)
    for i,m in enumerate(ms):
        p.setStyle({'model':i,},
                        {'stick':{'colorscheme':colors[i%len(colors)]}})
    p.zoomTo()
    return p.show()

def draw_grid(mols, labels=None, find_mcs=False, align_substructure=True):

    mols = [Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol for mol in mols]

    if not labels:
        labels = [Chem.MolToSmiles(mol) for mol in mols]

    if not find_mcs:
        return Chem.Draw.MolsToGridImage(mols, legends=labels)

    else:

        res = rdFMCS.FindMCS(mols)
        # mcs_smarts = res.smartsString
        mcs_mol = Chem.MolFromSmarts(res.smartsString)
        smarts = res.smartsString
        smart_mol = Chem.MolFromSmarts(smarts)
        smarts_and_mols = [smart_mol] + mols

        smarts_legend = "Max. substructure match"

        legends =  [smarts_legend] + labels

        matches = [""] + [mol.GetSubstructMatch(mcs_mol) for mol in mols]

        subms = [x for x in smarts_and_mols if x.HasSubstructMatch(mcs_mol)]

        AllChem.Compute2DCoords(mcs_mol)

        if align_substructure:
            for m in subms:
                _ = AllChem.GenerateDepictionMatching2DStructure(m, mcs_mol)

        drawing = Chem.Draw.MolsToGridImage(smarts_and_mols, highlightAtomLists=matches, legends=legends)

        return drawing

def draw_mcs(smiles: list[str] or dict[str, str], align_substructure: bool = True, verbose: bool = False, **kwargs):
     """
     Convert a list (or dictionary) of SMILES strings to an RDKit grid image of the maximum common substructure (MCS) match between them

     :returns: RDKit grid image, and (if verbose=True) MCS SMARTS string and molecule, and list of molecules for input SMILES strings
     :rtype: RDKit grid image, and (if verbose=True) string, molecule, and list of molecules
     :param molecules: The SMARTS molecules to be compared and drawn
     :type molecules: List of (SMARTS) strings, or dictionary of (SMARTS) string: (legend) string pairs
     :param align_substructure: Whether to align the MCS substructures when plotting the molecules; default is True
     :type align_substructure: boolean
     :param verbose: Whether to return verbose output (MCS SMARTS string and molecule, and list of molecules for input SMILES strings); default is False so calling this function will present a grid image automatically
     :type verbose: boolean
     """
     mols = [Chem.MolFromSmiles(smile) for smile in smiles]
     res = rdFMCS.FindMCS(mols, **kwargs)
     mcs_smarts = res.smartsString
     mcs_mol = Chem.MolFromSmarts(res.smartsString)
     smarts = res.smartsString
     smart_mol = Chem.MolFromSmarts(smarts)
     smarts_and_mols = [smart_mol] + mols

     smarts_legend = "Max. substructure match"

     # If user supplies a dictionary, use the values as legend entries for molecules
     if isinstance(smiles, dict):
          mol_legends = [smiles[molecule] for molecule in smiles]
     else:
          mol_legends = ["" for mol in mols]

     legends =  [smarts_legend] + mol_legends
    
     matches = [""] + [mol.GetSubstructMatch(mcs_mol) for mol in mols]

     subms = [x for x in smarts_and_mols if x.HasSubstructMatch(mcs_mol)]

     AllChem.Compute2DCoords(mcs_mol)

     if align_substructure:
          for m in subms:
               _ = AllChem.GenerateDepictionMatching2DStructure(m, mcs_mol)

     drawing = Draw.MolsToGridImage(smarts_and_mols, highlightAtomLists=matches, legends=legends)

     if verbose:
          return drawing, mcs_smarts, mcs_mol, mols
     else:
          return drawing
