
class Chain:

  def __init__(self,name):

    self._name = name
    self.residues = []
    self.fix_names()

  @property
  def name(self):
    return self._name

  @name.setter
  def name(self,name):
    self._name = name
    self.fix_names()

  def fix_names(self):
    for residue in self.residues:
      residue.chain = self._name
    for atom in self.atoms:
      atom.chain = self._name

  @property
  def num_residues(self):
    return len(self.residues)

  def add_residue(self,residue):
    import mcol
    import mout
    from .residue import Residue
    assert isinstance(residue,Residue)
    self.residues.append(residue)
    if self.num_residues > 1 and self.residues[-1].type != self.type:
      mout.errorOut(f'Differing residue types in same chain {self.residues[-1]} ({self.residues[-1].type}), {self.residues[0]} ({self.type})')

  @property
  def num_atoms(self):
    num_atoms = 0
    for residue in self.residues:
      num_atoms += residue.num_atoms
    return num_atoms

  @property
  def type(self):
    return self.residues[0].type

  @property
  def res_names(self):
    return [res.name for res in self.residues]

  def atom_names(self,wRes=False,noPrime=False):
    names_list = []
    for residue in self.residues:
      names_list.append(residue.atom_names(wRes=wRes,noPrime=noPrime))
    return names_list

  @property
  def atomic_numbers(self):
    number_list = []
    for residue in self.residues:
      number_list.append(residue.atomic_numbers)
    return number_list

  @property
  def positions(self):
    positions_list = []
    for residue in self.residues:
      positions_list += residue.positions
    return positions_list

  @property
  def charges(self):
    charges = []
    for residue in self.residues:
      charges += residue.charges
    return charges

  @property
  def masses(self):
    masses = []
    for residue in self.residues:
      masses += residue.masses
    return masses

  @property
  def atoms(self):
    atoms = []
    for residue in self.residues:
      atoms += residue.atoms
    return atoms

  @property
  def FF_atomtypes(self):
    atomtype_list = []
    for residue in self.residues:
      atomtype_list += residue.FF_atomtypes
    return atomtype_list

  def index_from_name(self,namestring):
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
    if fast:
      new_chain = Chain(self.name)
      for residue in self.residues:
        new_chain.add_residue(residue.copy(fast=fast))
      return new_chain
    else:
      import copy
      return copy.deepcopy(self)

  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name
  def __len__(self):
    return len(self.residues)