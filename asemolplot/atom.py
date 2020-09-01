
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import copy

from ase.data import atomic_numbers as ase_atomic_numbers # Atom only

class Atom:

  def __init__(self,
               name,
               index,
               pdb_index,
               position,
               residue,
               chain=None,
               res_number=None,
               charge=0.0,
               FF_atomtype=None,
               mass=None,
               LJ_sigma=None,
               LJ_epsilon=None,
               QM=False,
               occupancy=None,
               temp_factor=None,
               heterogen=None,
               charge_str=None):

    self._name = name
    self._atomic_number = None
    self.species = name[0]
    self.index = index
    self.pdb_index = pdb_index
    self._position = position
    self.residue = residue
    self.chain=chain
    self.res_number = res_number
    self.charge = charge
    self.FF_atomtype = FF_atomtype
    self.mass = mass
    self.LJ_sigma = LJ_sigma
    self.LJ_epsilon = LJ_epsilon
    self.QM = QM
    self.occupancy = occupancy
    self.temp_factor = temp_factor
    self.heterogen = heterogen
    self.charge_str = charge_str

  def print(self):

    mout.varOut("Atom Name",self._name)
    mout.varOut("Atom Species",self.species)
    mout.varOut("Atom (Residue) Index",self.index)
    mout.varOut("Atom (PDB) Index",self.pdb_index)
    mout.varOut("Atom Position",self._position)
    mout.varOut("Atom Residue",self.residue)
    mout.varOut("Atom Chain",self.chain)
    mout.varOut("Atom Residue Number",self.res_number)

  def get_name(self,wRes=False,noPrime=False):
    if wRes:
      namestring = self.residue+"_"+self._name
    else:
      namestring = self._name
    if noPrime:
      return namestring.replace("'", "P")
    else:
      return namestring

  @property
  def name(self):
    return self.get_name()

  @property
  def atomic_number(self):
    return ase_atomic_numbers[self.species]
  
  @property
  def position(self):
    return self._position

  @property
  def x(self):
    return self._position[0]

  @property
  def y(self):
    return self._position[1]

  @property
  def z(self):
    return self._position[2]

  @position.setter
  def position(self,pos):
    self._position = pos

  def copy(self):
    return copy.deepcopy(self)