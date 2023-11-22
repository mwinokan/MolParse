
import mout

class Atom:
  """Fundamental unit for molecular systems.

  These objects should not be created by the user, 
  but constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""
  
  name_symbol_dict = {'MG': 'Mg','Mg': 'Mg', 'LIT': 'Li', 'Li': 'Li', 'SOD': 'Na', 'POT':'K'}

  def __init__(self,name,index=None,pdb_index=None,position=None,residue=None,chain=None,res_number=None,charge=None,FF_atomtype=None,mass=None,LJ_sigma=None,LJ_epsilon=None,occupancy=None,temp_factor=None,heterogen=None,charge_str=None,velocity=None,alternative_site=None,res_index=None,element=None):

    # necessary upon init
    self._name = name
    self.index = index
    if pdb_index is not None:
      self._NUMBER = int(pdb_index)
    else:
      self._NUMBER = None
    self._position = position
    self.residue = residue

    if element:
      self._element = element
    else:
      mout.warning('Guessing element from first character of atom name!')
      self._element = name[0]
    assert self.symbol is not None

    self.chain=chain
    self.chain_number = None
    self._atomic_number = None
    if res_number is not None:
      self._res_number = int(res_number)
    else:
      self._res_number = None
    if res_index:
      self.res_index = int(res_index)
    else:
      self.res_index = None
    
    # optional
    self.FF_atomtype = FF_atomtype
    
    if charge:
      self.charge = charge
    elif charge_str:
      if charge_str.endswith('-'):
        self.charge = -1.0*float(charge_str[:-1])
      else:
        self.charge = float(charge_str[:-1])
    else:
      self.charge = 0.0

    self._mass = mass
    self.LJ_sigma = LJ_sigma
    self.LJ_epsilon = LJ_epsilon
    self.occupancy = occupancy
    self.temp_factor = temp_factor
    self.heterogen = heterogen
    self.charge_str = charge_str
    self.ter_line = None
    self.terminal = None
    self._alternative_site = alternative_site
    self._velocity = velocity

    self._parent = None

  @property
  def element(self):
    return self._element

  @property
  def symbol(self):
    return self._element

  @property
  def species(self):
    return self._element

  @property
  def alternative_site(self):
    return self._alternative_site
  
  @property
  def parent(self):
    return self._parent

  @parent.setter
  def parent(self,obj):
    # import weakref
    # self._parent = weakref.ref(obj)
    self._parent = obj

  def __deepcopy__(self, memodict={}):
    copy_object = Atom(self.name, self.index, self.pdb_index, self.position, self.residue, res_number=self.res_number, res_index=self.res_index, element=self.element)
    copy_object.chain = self.chain
    copy_object.occupancy = self.occupancy
    copy_object.temp_factor = self.temp_factor
    copy_object.heterogen = self.heterogen
    copy_object.charge_str = self.charge_str
    copy_object.charge = self.charge
    copy_object.velocity = self.velocity
    copy_object.chain_number = self.chain_number
    copy_object._alternative_site = self.alternative_site

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
    print(f'Atom {self.name} ({self.element}), index={self.index}, pdb_index={self.pdb_index}, res={self.residue}, res_number={self.res_number}, chain={self.chain}, chain_number={self.chain_number}')

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

  def set_name(self,name,element,verbosity=1):
    """Rename the atom"""
    if verbosity > 0:
      import mout
      import mcol
      mout.out("Renaming atom "+
               mcol.arg+self.name+str([self.index])+mcol.clear+" of residue "+
               mcol.varName+self.residue+str([self.res_number])+
               mcol.clear+" to "+mcol.arg+name)
    self._name = name
    
    assert element
    self._element = element
    assert self.symbol is not None

  @property
  def res_number(self):
    """number or pdb_index are the fixed residue index from the original structure file"""
    return self._res_number
  
  @property
  def number(self):
    """number or pdb_index are the fixed atom index from the original structure file"""
    return self._NUMBER
  
  @property
  def pdb_index(self):
    """number or pdb_index are the fixed atom index from the original structure file"""
    return self._NUMBER

  @property
  def name_number_str(self):
    return f'{self.name} {self.number}'

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
  def covalent_radius(self):
    """Covalent radius (float) property"""
    from ase.data import covalent_radii as ase_covalent_radii # Atom only
    return ase_covalent_radii[self.atomic_number]

  @property
  def vdw_radius(self):
    """vdW radius (float) property"""
    from ase.data import vdw_radii as ase_vdw_radii # Atom only
    return ase_vdw_radii[self.atomic_number]

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
    return np.array((self.velocity[0] ,self.velocity[1], self.velocity[2]))

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
    if isinstance(other,Atom):
      return self.np_pos - other.np_pos
    else:
      return self.np_pos - other

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

  # def __del__(self):
  #   mout.debug(f'Atom({self},id({id(self)})) was deleted')
  