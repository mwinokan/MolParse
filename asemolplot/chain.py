
class Chain:
  """Class containing covalently bonded Residue objects

  These objects should not be created by the user, 
  but constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""

  def __init__(self,name):
    self._name = name
    self.residues = []
    self.fix_names()

  @property
  def name(self):
    """Name (str) property"""
    return self._name

  @name.setter
  def name(self,name):
    self._name = name
    self.fix_names()

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
    self.residues.append(residue)
    if self.num_residues > 1 and self.residues[-1].type != self.type:
      mout.errorOut(f'Differing residue types in same chain {self.residues[-1]} ({self.residues[-1].type}), {self.residues[0]} ({self.type})')

  @property
  def num_atoms(self):
    """Number of child Atoms (int)"""
    num_atoms = 0
    for residue in self.residues:
      num_atoms += residue.num_atoms
    return num_atoms

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
      names_list.append(residue.atom_names(wRes=wRes,noPrime=noPrime))
    return names_list

  @property
  def atomic_numbers(self):
    """All child atomic numbers (list)"""
    number_list = []
    for residue in self.residues:
      number_list.append(residue.atomic_numbers)
    return number_list

  @property
  def positions(self):
    """All child Atom positions (list)"""
    positions_list = []
    for residue in self.residues:
      positions_list += residue.positions
    return positions_list

  @property
  def charges(self):
    """All child Atom charges (list)"""
    charges = []
    for residue in self.residues:
      charges += residue.charges
    return charges

  @property
  def masses(self):
    """All child Atom masses (list)"""
    masses = []
    for residue in self.residues:
      masses += residue.masses
    return masses

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
  @property
  def __repr__(self):
    return self.name
  @property
  def __str__(self):
    return self.name
  @property
  def __len__(self):
    return len(self.residues)
    