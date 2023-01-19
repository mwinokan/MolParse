from ase import io

import mcol              # https://github.com/mwinokan/MPyTools
import mout              # https://github.com/mwinokan/MPyTools
import molparse as amp # https://github.com/mwinokan/AseMolPlot

# # load some molecule
# from ase.data.pubchem import pubchem_atoms_search, pubchem_atoms_conformer_search
# atoms = pubchem_atoms_search(smiles='N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O')

# custom_style = amp.styles.standard
# custom_style['rotation'] = ''
# custom_style['canvas_width'] = 500

# # make some single images:
# amp.makeImage('amp',atoms,**amp.styles.standard)
# amp.makePovImage('pov',atoms,**amp.styles.standard)

# system = amp.parsePDB("AT_PBE_GOLD.pdb")
system = amp.parsePDB("DNA9_AT_TER.pdb")

import cs4qb

# ff,system = cs4qb.buildGulpFF("AT_PBE_GOLD.pdb",output_ff="custom.ff",output_cjson="system.cjson",output_psf="custom.psf")

# system.check_indices()

amp.write("test.pdb",system)
