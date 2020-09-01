
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import copy

from .atom import Atom

class Residue:
  def __init__(self,name,number=None,chain=None):

    self.name = name
    self.chain = chain
    self.number = number
    self._atoms = []

  @property
  def num_atoms(self):
    return len(self._atoms)

  def atom_names(self,wRes=False,noPrime=False):
    names_list = []
    for atom in self._atoms:
      names_list.append(atom.get_name(wRes=wRes,noPrime=noPrime))
    return names_list

  @property
  def positions(self):
    positions_list = []
    for atom in self._atoms:
      positions_list.append(atom.position)
    return positions_list

  @property
  def atomic_numbers(self):
    number_list = []
    for atom in self._atoms:
      number_list.append(atom.atomic_number)
    return number_list

  @property
  def FF_atomtypes(self):
    atomtype_list = []
    for atom in self._atoms:
      atomtype_list.append(atom.FF_atomtype)
    return atomtype_list

  def addAtom(self,atom):
    self._atoms.append(atom)

  def print(self):

    mout.varOut("Residue Name",self.name)
    mout.varOut("Residue Chain",self.chain)
    mout.varOut("Residue Number",self.number)
    mout.varOut("Number of Atoms",self.num_atoms)

  @property
  def species(self):
    species_list = []
    for atom in self._atoms:
      species_list.append(atom.species)
    return ''.join(species_list)
  
  @property
  def charges(self):
    charges = []
    for atom in self._atoms:
      charges.append(atom.charge)
    return charges

  @property
  def atoms(self):
    return self._atoms

  def get_atom(self,name=None):
    for atom in self._atoms:
      if atom.name == name: return atom
    mout.errorOut("No atom "+
                  mcol.arg+name+
                  mcol.error+" in residue "+
                  mcol.arg+self.name,fatal=True,code="Residue.1")

  def index_from_name(self,name):
    for index,atom in enumerate(self._atoms):
      if atom.name == name: return index
    mout.errorOut("No atom "+
                  mcol.arg+name+
                  mcol.error+" in residue "+
                  mcol.arg+self.name,fatal=True,code="Residue.2")

  def copy(self):
    return copy.deepcopy(self)

  def get_distance(self,i,j):
    if type(i) is str:
      i = self.get_atom(name=i)
    elif type(i) is int:
      i = self._atoms[i]
    if type(j) is str:
      j = self.get_atom(name=j)
    elif type(j) is int:
      j = self._atoms[j]

    return displacement(i.position,j.position)

  def get_angle(self,i,j,k):
    if type(i) is str:
      i = self.get_atom(name=i)
    elif type(i) is int:
      i = self._atoms[i]
    if type(j) is str:
      j = self.get_atom(name=j)
    elif type(j) is int:
      j = self._atoms[j]
    if type(k) is str:
      k = self.get_atom(name=k)
    elif type(k) is int:
      k = self._atoms[k]

    return angle(i.position,j.position,k.position)