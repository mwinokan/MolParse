
class Atom:
  """Fundamental unit for molecular systems.

  These objects should not be created by the user, 
  but constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""

  def __init__(self,name,index=None,pdb_index=None,position=None,residue=None,chain=None,res_number=None,charge=0.0,FF_atomtype=None,mass=None,LJ_sigma=None,LJ_epsilon=None,QM=False,occupancy=None,temp_factor=None,heterogen=None,charge_str=None,velocity=None):

    # necessary upon init
    self._name = name
    self.index = index
    self.pdb_index = pdb_index
    self._position = position
    self.residue = residue

    name_symbol_dict = {'MG': 'Mg'}

    if name in name_symbol_dict.keys():
      self.species = name_symbol_dict[name]
      self.symbol = name_symbol_dict[name]
    else:
      self.species = name[0]
      self.symbol = name[0]

    self.chain=chain
    self.chain_number = None
    self._atomic_number = None
    if res_number:
      self.res_number = int(res_number)
    else:
      self.res_number = None
    
    # optional
    self.FF_atomtype = FF_atomtype
    
    self.charge = charge
    self._mass = mass
    self.LJ_sigma = LJ_sigma
    self.LJ_epsilon = LJ_epsilon
    self.QM = QM
    self.occupancy = occupancy
    self.temp_factor = temp_factor
    self.heterogen = heterogen
    self.charge_str = charge_str
    self.ter_line = None
    self.terminal = None
    self._velocity = velocity

  def __deepcopy__(self, memodict={}):
    copy_object = Atom(self.name, self.index, self.pdb_index, self.position, self.residue)
    copy_object.chain = self.chain
    copy_object.res_number = self.res_number
    copy_object.QM = self.QM
    copy_object.occupancy = self.occupancy
    copy_object.temp_factor = self.temp_factor
    copy_object.heterogen = self.heterogen
    copy_object.charge_str = self.charge_str
    copy_object.charge = self.charge
    copy_object.velocity = self.velocity
    copy_object.chain_number = self.chain_number

    return copy_object

  def print(self):
    """Verbose output of the atom's properties"""
    import mout
    mout.varOut("Atom Name",self._name)
    mout.varOut("Atom Species",self.species)
    mout.varOut("Atom (ASE) Index",self.index)
    mout.varOut("Atom (PDB) Index",self.pdb_index)
    mout.varOut("Atom Position",self._position)
    mout.varOut("Atom Residue",self.residue)
    mout.varOut("Atom Chain",self.chain)
    mout.varOut("Atom Residue Number",self.res_number)

  def summary(self):
    """Summarised output of the atom's properties"""
    print(f'Atom {self.name}, index={self.index}, pdb_index={self.pdb_index}, res={self.residue}, res_number={self.res_number}, chain={self.chain}, chain_number={self.chain_number}')

  def get_name(self,wRes=False,noPrime=False):
    """Returns a string of resname_atomname"""
    if wRes:
      namestring = self.residue+"_"+self._name
    else:
      namestring = self._name
    if noPrime:
      return namestring.replace("'", "P")
    else:
      return namestring

  def set_name(self,name,verbosity=1):
    """Rename the atom"""
    if verbosity > 0:
      import mout
      import mcol
      mout.out("Renaming atom "+
               mcol.arg+self.name+str([self.index])+mcol.clear+" of residue "+
               mcol.varName+self.residue+str([self.res_number])+
               mcol.clear+" to "+mcol.arg+name)
    self._name = name
    self.species = name[0]

  @property
  def name(self):
    """Name (str) property"""
    return self.get_name()

  @name.setter
  def name(self,name):
    self.set_name(name=name,verbosity=0)

  @property
  def atomic_number(self):
    """Atomic number (int) property"""
    from ase.data import atomic_numbers as ase_atomic_numbers # Atom only
    return ase_atomic_numbers[self.species]

  @property
  def mass(self):
    """Atomic mass (float) property"""
    if self._mass is None:
      from ase.data import atomic_masses as ase_atomic_masses # Atom only
      return ase_atomic_masses[self.atomic_number]
    else:
      return self._mass
  
  @property
  def position(self):
    """Cartesian position (list)"""
    return self._position

  @property
  def velocity(self):
    """Cartesian velocity (list)"""
    return self._velocity

  @property
  def np_pos(self):
    """Cartesian position (numpy array)"""
    import numpy as np
    return np.array((self._position[0] ,self._position[1], self._position[2]))

  @property
  def np_vel(self):
    """Cartesian velocity (numpy array)"""
    import numpy as np
    return np.array((self._velocity[0] ,self._velocity[1], self._velocity[2]))

  @property
  def x(self):
    """X-coordinate"""
    return self._position[0]

  @property
  def y(self):
    """Y-coordinate"""
    return self._position[1]

  @property
  def z(self):
    """Z-coordinate"""
    return self._position[2]

  @position.setter
  def position(self,pos):
    self._position = pos

  @velocity.setter
  def velocity(self,vel):
    self._velocity = vel

  def copy(self):
    """Returns a deepcopy of the Atom object"""
    import copy
    return copy.deepcopy(self)

  def is_in_residue(self,residue):
    from .residue import Residue
    assert isinstance(residue, Residue)

    # print(f"self: {self.residue}, {self.res_number}, {self.chain}")
    # print(f"other: {residue.name}, {residue.number}, {residue.chain}")

    if self.residue != residue.name:
      return False
    if self.res_number != residue.number:
      return False
    if self.chain != residue.chain:
      return False
    return True

  @property
  def type(self):
    from .residue import res_type
    return res_type(self.residue)

  @property
  def ase_atom(self):
    from ase import Atom
    return Atom(self.species,self.position)

  def __float__(self):
    return self.np_pos

  def __sub__(self,other):
    assert isinstance(other,Atom)
    return self.np_pos - other.np_pos

  def __add__(self,other):
    import numpy as np
    if isinstance(other,np.ndarray):
      return self.np_pos + other
    if isinstance(other,Atom):
      return self.np_pos + other.np_pos
    raise TypeError

  @property
  def children(self):
    return None

  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name