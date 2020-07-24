
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

from ase.data import atomic_numbers as ase_atomic_numbers

class System:

  def __init__(self,name):

    self.name = name
    self.chains = []

  def addChain(self,chain):
    self.chains.append(chain)

  @property
  def num_atoms(self):
    num_atoms = 0
    for chain in self.chains:
      num_atoms += chain.num_atoms
    return num_atoms

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
  def residues(self):
    residues = []
    for chain in self.chains:
      residues += chain.residues
    return residues

  @property
  def FF_atomtypes(self):
    atomtype_list = []
    for chain in self.chains:
      atomtype_list += chain.FF_atomtypes
    return atomtype_list

  def write_CJSON(self,filename,use_atom_types=False,gulp_names=False):
    from .io import writeCJSON
    writeCJSON(filename,self,use_atom_types=use_atom_types,gulp_names=gulp_names)

class Chain:

  def __init__(self,name):

    self.name = name
    self.residues = []

  def addResidue(self,residue):
    self.residues.append(residue)

  @property
  def num_atoms(self):
    num_atoms = 0
    for residue in self.residues:
      num_atoms += residue.num_atoms
    return num_atoms

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
                  mcol.arg+self.name,fatal=True)

class Atom:

  def __init__(self,
               name,
               index,
               pdb_index,
               position,
               residue,
               chain=None,
               res_number=None,
               charge=0.0,FF_atomtype=None):

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
  