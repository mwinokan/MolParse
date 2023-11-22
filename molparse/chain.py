
from .group import AtomGroup

class Chain(AtomGroup):
  """Class containing covalently bonded Residue objects

  These objects should not be created by the user, 
  but constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""

  def __init__(self,name):

    super(Chain, self).__init__(name)

    from .list import NamedList
    self.residues = NamedList()

    self.index = None
    self.fix_names()

    self._parent = None

  @property
  def parent(self):
    return self._parent

  @parent.setter
  def parent(self,obj):
    # import weakref
    # self._parent = weakref.ref(obj)
    self._parent = obj

  def fix_names(self):
    """Ensure child-classes have correct parent name"""
    for residue in self.residues:
      residue.chain = self._name
    for atom in self.atoms:
      atom.chain = self._name

  @property
  def num_residues(self):
    """Number of child residues (int)"""
    return len(self.residues)

  def add_residue(self,residue):
    """Add a Residue to the chain"""
    import mcol
    import mout
    from .residue import Residue
    assert isinstance(residue,Residue)
    residue.set_chain_number(self.index)
    residue.set_chain_char(self.name)
    residue.parent = self

    # remove any termini
    if self.residues:
      for atom in self.residues[-1].atoms:
        atom.terminal = None
        atom.ter_line = None

    self.residues.append(residue)
    if self.num_residues > 1 and self.residues[-1].type != self.type:
      mout.errorOut(f'Differing residue types in same chain {self.residues[-1]} ({self.residues[-1].type}), {self.residues[0]} ({self.type})')

  @property
  def type(self):
    """Classification of Residues in this Chain"""
    return self.residues[0].type

  @property
  def res_names(self):
    """All child Residue names (list)"""
    return [res.name for res in self.residues]

  def atom_names(self,wRes=False,noPrime=False):
    """All child Atom names (list)"""
    names_list = []
    for residue in self.residues:
      names_list += residue.atom_names(wRes=wRes,noPrime=noPrime)
    return names_list

  @property
  def atoms(self):
    """All child Atom objects (list)"""
    atoms = []
    for residue in self.residues:
      atoms += residue.atoms
    return atoms

  @property
  def FF_atomtypes(self):
    """All child Atom atomtypes (list)"""
    atomtype_list = []
    for residue in self.residues:
      atomtype_list += residue.FF_atomtypes
    return atomtype_list

  def index_from_name(self,namestring):
    """Find the index of an atom by its name.
    Assumes formatting resname_atomname"""
    import mcol
    import mout

    search_residue, search_atom = namestring.split("_")

    for index,atom in enumerate(self.atoms):
      if atom.residue == search_residue and atom.name == search_atom:
        return index

    mout.errorOut("Atom "+
              mcol.arg+search_atom+
              mcol.error+" could not be found in residue"+
              mcol.arg+search_residue+" of chain "+mcol.arg+self.name,fatal=True,code="Chain.1")

  def copy(self,fast=False):
    """return a deepcopy of the Chain"""
    if fast:
      new_chain = Chain(self.name)
      for residue in self.residues:
        new_chain.add_residue(residue.copy(fast=fast))
      return new_chain
    else:
      import copy
      return copy.deepcopy(self)

  def add_atom(self,atom):
    """add an atom to the chain"""
    from .residue import Residue

    if self.residues and atom.is_in_residue(self.residues[-1]): 
      self.residues[-1].add_atom(atom)
    else:

      residue = Residue(atom.residue,number=atom.res_number)
      residue.index = atom.res_index
      residue.chain = atom.chain
      residue.add_atom(atom)
      self.add_residue(residue)

  @property
  def children(self):
    return self.residues
    