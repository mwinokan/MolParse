
# from importlib import reload

class Residue:
  def __init__(self,name,number=None,chain=None):

    self._name = name
    self.chain = chain
    self.number = number
    self._atoms = []

  def rename(self,new,verbosity=1):
    if verbosity > 0:
      import mcol
      import mout
      mout.out("Renaming residue "+
               mcol.arg+self.name+str([self.number])+mcol.clear+" of chain "+
               mcol.varName+self.chain+
               mcol.clear+" to "+mcol.arg+new)
    self._name = new
    self.fix_names()

  @property
  def num_atoms(self):
    return len(self._atoms)

  @property
  def name(self):
    return self._name

  @name.setter
  def name(self,name):
    self._name = name
    self.fix_names()

  def fix_names(self):
    for atom in self.atoms:
      atom.residue = self._name

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

  @property
  def type(self):
    if self.name.startswith(('DA','DT','DC','DG')):
      this_type = "DNA"
    elif self.name.startswith(('SOL','WAT','TIP','T3P')):
      this_type = "SOL"
    elif self.name.startswith(('ION','MG','CL','NA','SOD','POT','CAL','LIT')):
      this_type = "ION"
    elif self.name.startswith(('DPPC','POPC')):
      this_type = "LIP"
    elif self.name.startswith(("ALA","ARG","ASN","ASP","ATP",
                               "CYS","GLN","GLU","GLY","HSD",
                               "HSE","HIS","ILE","LEU","LYS",
                               "MET","PHE","PRO","SER","THR",
                               "TRP","TYR","VAL","HID")):
      this_type = "PRO"
    else:
      import mcol
      import mout
      mout.warningOut("Unknown residue type for "+mcol.arg+self.name)
      this_type = None
    return this_type

  def addAtom(self,atom):
    from .atom import Atom
    assert isinstance(atom,Atom)
    self._atoms.append(atom)

  def translate(self,vector):
    for atom in self._atoms:
      atom.position = atom.np_pos + vector

  def print(self):
    import mout
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
  def masses(self):
    masses = []
    for atom in self._atoms:
      masses.append(atom.mass)
    return masses

  @property
  def atoms(self):
    return self._atoms

  def get_atom(self,name):
    for atom in self._atoms:
      if atom.name == name: return atom

    import mout
    import mcol
    mout.errorOut("No atom "+
                  mcol.arg+name+
                  mcol.error+" in residue "+
                  mcol.arg+self.name+str([self.number]),fatal=False,code="Residue.1")
    return None

  def delete_atom(self,name,verbosity=1):
    import mcol
    import mout
    for index,atom in enumerate(self._atoms):
      if atom.name == name:
        del self._atoms[index]
        if verbosity > 0:
          mout.warningOut("Deleted "+mcol.result+"1"+
                          mcol.warning+" atom of name "+
                          mcol.arg+name+
                          mcol.warning+" in residue "+
                          mcol.arg+self.name)
        return
    if verbosity > 1:
      mout.warningOut("Deleted "+mcol.result+"0"+
                      mcol.warning+" atom of name "+
                      mcol.arg+name+
                      mcol.warning+" in residue "+
                      mcol.arg+self.name)

  def copy(self):
    import copy
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

  def centre_of_mass(self,verbosity=1):
    return self.CoM(verbosity=verbosity)

  def summary(self):
    import mout
    mout.headerOut(f'\nResidue {self.name}, index={self.number}, #atoms={self.num_atoms}')

    for atom in self.atoms:
      print(atom.name,atom.index,atom.pdb_index)
    
  def CoM(self,verbosity=1):
    import numpy as np

    position_list = self.positions

    centre_of_mass = np.array([sum([pos[0] for pos in position_list])/len(position_list),
                               sum([pos[1] for pos in position_list])/len(position_list),
                               sum([pos[2] for pos in position_list])/len(position_list)])

    if verbosity > 0: 
      mout.varOut("CoM of "+self.name,
                  centre_of_mass,
                  unit="Angstroms",precision=4)

    return centre_of_mass

  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name
  def __len__(self):
    return len(self.atoms)
