
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import copy

from .chain import Chain
# from .bondlist import Connectivity

class System:

  def __init__(self,name):

    self.name = name
    self.chains = []

    self.bondlist=None

  def check_indices(self):

    for index,atom in enumerate(self.atoms):
      if index != atom.index:
        print(index,atom.index,atom.name,atom.residue)

  def fix_indices(self):
    for index,atom in enumerate(self.atoms):
      atom.index = index

  def addChain(self,chain):
    self.chains.append(chain)

  @property
  def num_atoms(self):
    num_atoms = 0
    for chain in self.chains:
      num_atoms += chain.num_atoms
    return num_atoms

  @property
  def num_chains(self):
    return len(self.chains)
    
  def atom_names(self,wRes=False,noPrime=False):
    names_list = []
    for chain in self.chains:
      names_list.append(chain.atom_names(wRes=wRes,noPrime=noPrime))
    return names_list

  @property
  def atomic_numbers(self):
    number_list = []
    for chain in self.chains:
      number_list.append(chain.atomic_numbers)
    return number_list

  @property
  def positions(self):
    positions_list = []
    for chain in self.chains:
      positions_list += chain.positions
    return positions_list

  @property
  def charges(self):
    charges = []
    for chain in self.chains:
      charges += chain.charges
    return charges

  @property
  def atoms(self):
    atoms = []
    for chain in self.chains:
      atoms += chain.atoms
    return atoms

  @property
  def QM_indices(self):
    index_list = []
    for index,atom in enumerate(self.atoms):
      if atom.QM:
        index_list.append(index)
    return index_list

  @property
  def residues(self):
    residues = []
    for chain in self.chains:
      residues += chain.residues
    return residues

  @property
  def res_names(self):
    return [res.name for res in self.residues]

  @property
  def num_residues(self):
    return len(self.residues)

  @property
  def FF_atomtypes(self):
    atomtype_list = []
    for chain in self.chains:
      atomtype_list += chain.FF_atomtypes
    return atomtype_list

  def write_CJSON(self,filename,use_atom_types=False,gulp_names=False):
    from .io import writeCJSON
    writeCJSON(filename,self,use_atom_types=use_atom_types,gulp_names=gulp_names)

  def set_coordinates(self,reference):

    if type(reference) is str:
      from ase.io import read
      atoms = read(reference)
    else:
      atoms = reference

    for index,atom in enumerate(self.atoms):
      atom.position = atoms[index].position

  def copy(self):
    return copy.deepcopy(self)
