
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import copy

from .chain import Chain
# from .bondlist import Connectivity

import numpy as np

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

    for index,residue in enumerate(self.residues):
      residue.number = index

  def add_chain(self,chain):
    self.chains.append(chain)

  def add_system(self,system):

    for chain in system.chains:
      self.add_chain(chain)

    self.fix_indices()

  @property
  def num_atoms(self):
    num_atoms = 0
    for chain in self.chains:
      num_atoms += chain.num_atoms
    return num_atoms

  @property
  def charge(self):
    charge = 0
    for atom in self.atoms:
      charge += atom.charge
    return charge

  def summary(self,res_limit=10):
    mout.headerOut("\nSystem "+mcol.arg+self.name+
                   mcol.clear+mcol.bold+" contains "+
                   mcol.result+str(self.num_chains)+
                   mcol.clear+mcol.bold+" chains:")
    for chain in self.chains:
      mout.headerOut("Chain "+mcol.result+chain.name+
                     mcol.clear+mcol.bold+" contains "+
                     mcol.result+str(chain.num_residues)+
                     mcol.clear+mcol.bold+" residues:")
      names = ""
      for name in chain.res_names[:res_limit]:
        names += name+" "
      if chain.num_residues > res_limit:
        names += "..."
      mout.out(names)
    # mout.varOut("Total Charge",self.charge)

  @property
  def num_chains(self):
    return len(self.chains)

  @property
  def chain_names(self):
    names = []
    for chain in self.chains:
      names.append(chain.name)
    return names

  def rename_atoms(self,old,new,res_filter=None,verbosity=2):
    count=0
    for residue in self.residues:
      if res_filter is not None and res_filter not in residue.name:
        continue
      for atom in residue.atoms:
        if atom.name == old:
          count += 1
          atom.set_name(new,verbosity=verbosity-1)
    if verbosity > 0:
      if res_filter is None:
        mout.warningOut("Renamed "+mcol.result+str(count)+mcol.warning+" atoms from "+mcol.arg+old+
                        mcol.warning+" to "+mcol.arg+new)
      else:
        mout.warningOut("Renamed "+mcol.result+str(count)+mcol.warning+" atoms from "+mcol.arg+old+
                        mcol.warning+" to "+mcol.arg+new+mcol.warning+" with res_filter "+mcol.arg+res_filter)
    return count

  def rename_residues(self,old,new,verbosity=2):
    count=0
    for residue in self.residues:
      if residue.name == old:
        # residue.name = new
        residue.rename(new,verbosity=verbosity-1)
        count += 1
    if verbosity > 0:
      mout.warningOut("Renamed "+mcol.result+str(count)+mcol.warning+" residues from "+mcol.arg+old+
                      mcol.warning+" to "+mcol.arg+new)
    return count

  def get_chain(self,name):
    for chain in self.chains:
      if chain.name == name:
        return chain
    mout.errorOut("Chain with name "+mcol.arg+name+mcol.error+" not found.",fatal=True)

  def remove_chain(self,name,verbosity=1):
    for index,chain in enumerate(self.chains):
      if chain.name == name:
        if verbosity > 0:
          mout.warningOut("Removing chain "+mcol.arg+name+str([index]))
        del self.chains[index]
        return chain
    mout.errorOut("Chain with name "+mcol.arg+name+mcol.error+" not found.",fatal=True)
    
  def remove_heterogens(self,verbosity=1):
    del_list = []
    atoms = self.atoms
    for index,atom in enumerate(atoms):
      if atom.heterogen:
        del_list.append(index)
    number_deleted = self.remove_atoms_by_index(del_list,verbosity=verbosity-1)
    if verbosity > 0:
      mout.warningOut("Removed "+mcol.result+str(number_deleted)+mcol.warning+" heterogens")

  def remove_atoms_by_index(self,del_list,verbosity=1):
    number_deleted=0
    self.fix_indices()
    for chain in self.chains:
      for residue in chain.residues:
        for index,atom in reversed(list(enumerate(residue.atoms))):
          if atom.index in del_list:
            del residue.atoms[index]
            number_deleted += 1
            if verbosity > 1:
              mout.warningOut("Removed atom "+
                              mcol.arg+atom.name+str([atom.index])+
                              mcol.warning+" in residue "+
                              mcol.arg+residue.name+str([residue.number]))
    return number_deleted

  def remove_atoms(self,name,res_filter=None,verbosity=2):
    number_deleted=0
    self.fix_indices()
    for chain in self.chains:
      for residue in chain.residues:
        if res_filter is not None and res_filter not in residue.name:
          continue
        for index,atom in reversed(list(enumerate(residue.atoms))):
          if atom.name == name:
            del residue.atoms[index]
            number_deleted += 1
            if verbosity > 1:
              mout.warningOut("Removed atom "+
                              mcol.arg+atom.name+str([atom.index])+
                              mcol.warning+" in residue "+
                              mcol.arg+residue.name+str([residue.number]))
    if verbosity > 0:
      if res_filter is None:
        mout.warningOut("Removed "+mcol.result+str(number_deleted)+mcol.warning+" atoms named "+mcol.arg+name)
      else:
        mout.warningOut("Removed "+mcol.result+str(number_deleted)+mcol.warning+" atoms named "+mcol.arg+name+mcol.warning+" with res_filter "+mcol.arg+res_filter)
    return number_deleted

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

  def centre_of_mass(self,set=None,shift=None):
    return CoM(set=set,shift=shift)

  def CoM(self,set=None,shift=None):

    position_list = self.positions

    centre_of_mass = np.array([sum([pos[0] for pos in position_list])/len(position_list),
                               sum([pos[1] for pos in position_list])/len(position_list),
                               sum([pos[2] for pos in position_list])/len(position_list)])

    mout.varOut("CoM of "+self.name,
                centre_of_mass,
                unit="Angstroms",precision=4)

    if set is not None:

      assert np.ndim(set) == 1

      try:
        new_com = np.array([set[0],set[1],set[2]])
      except:
        mout.errorOut("Incompatible array input",code="amp.system.CoM.1")

      shift_com = new_com - centre_of_mass

    else:

      shift_com = np.array([0,0,0])

    if shift is not None:

      assert np.ndim(shift) == 1

      try:
        new_shift = np.array([shift[0],shift[1],shift[2]])
      except:
        mout.errorOut("Incompatible array input",code="amp.system.CoM.2")

      shift_com = shift_com + new_shift
      
    if set is not None or shift is not None:

      for atom in self.atoms:
        atom.position = atom.position + shift_com

    return centre_of_mass

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

  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name
