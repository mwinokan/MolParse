
import py3Dmol
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
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
