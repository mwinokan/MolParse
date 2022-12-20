
class Residue:
  """Class for chemical Residue

  These objects should not be created by the user, 
  but constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""

  def __init__(self,name: str,number: int=None,chain: str=None,atoms=None):
    self._name = name
    self.chain = chain
    self.number = number
    self._atoms = []
    
    if atoms:

      from .atom import Atom

      # import from ASE atoms object:

      for atom in atoms:
        self.add_atom(Atom(name=atom.symbol,position=atom.position,residue=self.name))

        # print(atom)

  def rename(self,new: str,verbosity: int=1):
    """Rename the residue"""
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
    """Number of child atoms (int)"""
    return len(self._atoms)

  @property
  def name(self):
    return self._name

  @name.setter
  def name(self,name: str):
    self._name = name
    self.fix_names()

  def fix_names(self):
    """Ensure child Atoms have correct parent name"""
    for atom in self.atoms:
      atom.residue = self.name

  def fix_indices(self):
    """Ensure child Atoms have correct res_number"""
    for atom in self.atoms:
      atom.res_number = self.number

  def atom_names(self,wRes: bool=False,noPrime: bool=False):
    """Returns names of all child Atoms (list)"""
    names_list = []
    for atom in self._atoms:
      names_list.append(atom.get_name(wRes=wRes,noPrime=noPrime))
    return names_list

  @property
  def positions(self):
    """Returns positions of all child Atoms (list)"""
    positions_list = []
    for atom in self._atoms:
      positions_list.append(atom.position)
    return positions_list

  @property
  def atomic_numbers(self):
    """Returns atomic numbers of all child Atoms (list)"""
    number_list = []
    for atom in self._atoms:
      number_list.append(atom.atomic_number)
    return number_list

  @property
  def FF_atomtypes(self):
    """Returns atomtypes of all child Atoms (list)"""
    atomtype_list = []
    for atom in self._atoms:
      atomtype_list.append(atom.FF_atomtype)
    return atomtype_list

  @property
  def type(self):
    """Guess type from residue name"""
    return res_type(self.name)

  def add_atom(self,atom):
    """add an Atom"""
    self.addAtom(atom)

  def addAtom(self,atom):
    """add an Atom"""
    from .atom import Atom
    assert isinstance(atom,Atom)
    self._atoms.append(atom)

  def translate(self,vector):
    """Adjust the position of all atoms in the residue"""
    for atom in self._atoms:
      atom.position = atom.np_pos + vector

  def print(self):
    """Verbose output of the Residue's properties"""
    import mout
    mout.varOut("Residue Name",self.name)
    mout.varOut("Residue Chain",self.chain)
    mout.varOut("Residue Number",self.number)
    mout.varOut("Number of Atoms",self.num_atoms)

  @property
  def species(self):
    """Returns species of all child Atoms (list)"""
    species_list = []
    for atom in self._atoms:
      species_list.append(atom.species)
    return ''.join(species_list)
  
  @property
  def charges(self):
    """Returns charges of all child Atoms (list)"""
    charges = []
    for atom in self._atoms:
      charges.append(atom.charge)
    return charges

  @property
  def masses(self):
    """Returns masses of all child Atoms (list)"""
    masses = []
    for atom in self._atoms:
      masses.append(atom.mass)
    return masses

  @property
  def indices(self):
    """Returns indices of all child Atoms (list)"""
    indices = []
    for atom in self._atoms:
      indices.append(atom.index)
    return indices

  @property
  def atoms(self):
    """Returns all child Atoms (list)"""
    return self._atoms

  def get_atom(self,name):
    """Get a child Atom by name"""
    if isinstance(name,list):
      return [self.get_atom(n) for n in name]

    for atom in self._atoms:
      if atom.name == name: return atom

    import mout
    import mcol
    mout.errorOut("No atom "+
                  mcol.arg+name+
                  mcol.error+" in residue "+
                  mcol.arg+self.name+str([self.number]),fatal=False,code="Residue.1")
    return None

  def delete_atom(self,name: str,verbosity: int=1):
    """Delete a child Atom by name"""
    import mcol
    import mout
    if isinstance(name,list):
      for n in name:
        self.delete_atom(n,verbosity=verbosity-1)
      return
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

  def copy(self,fast: bool=False):
    """Return a deepcopy of the Residue"""
    if fast:
      new_residue = Residue(self.name,self.number,self.chain)
      for atom in self.atoms:
        new_residue.addAtom(atom.copy(fast=fast))
      return new_residue
    else:
      import copy
      return copy.deepcopy(self)

  def get_distance(self,i: str,j: str):
    """Get distance between two named atoms"""
    if type(i) is str:
      i = self.get_atom(name=i)
    elif type(i) is int:
      i = self._atoms[i]
    if type(j) is str:
      j = self.get_atom(name=j)
    elif type(j) is int:
      j = self._atoms[j]

    return displacement(i.position,j.position)

  def get_angle(self,i: str,j: str,k: str):
    """Get angle between three named atoms"""
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

  def centre_of_mass(self,verbosity: int=1):
    """Calculate the centre of mass of the Residue"""
    return self.CoM(verbosity=verbosity)

  def summary(self):
    """Summarised output of the Residue"""
    import mout
    mout.headerOut(f'\nResidue {self.name}, index={self.number}, #atoms={self.num_atoms}')

    for atom in self.atoms:
      print(atom.name,atom.index,atom.pdb_index)
    
  def CoM(self,verbosity: int=1):
    """Centre of mass of the Residue"""
    import numpy as np
    import mout

    position_list = self.positions

    centre_of_mass = np.array([sum([pos[0] for pos in position_list])/len(position_list),
                               sum([pos[1] for pos in position_list])/len(position_list),
                               sum([pos[2] for pos in position_list])/len(position_list)])

    if verbosity > 0: 
      mout.varOut("CoM of "+self.name,
                  centre_of_mass,
                  unit="Angstroms",precision=4)

    return centre_of_mass

  def is_same_as(self,residue):
    assert isinstance(residue, Residue)

    if self.name != residue.name:
      return False
    if self.index != residue.index:
      return False
    if self.chain != residue.chain:
      return False
    return True

  def set_chain_number(self,index):
    for atom in self.atoms:
      atom.chain_number = index

  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name
  def __len__(self):
    return len(self.atoms)

  @property
  def bbox(self):
    """Bounding box of the residue"""
    import numpy as np
    x = [min([a.np_pos[0] for a in self.atoms]),max([a.np_pos[0] for a in self.atoms])]
    y = [min([a.np_pos[1] for a in self.atoms]),max([a.np_pos[1] for a in self.atoms])]
    z = [min([a.np_pos[2] for a in self.atoms]),max([a.np_pos[2] for a in self.atoms])]
    return [x,y,z]

  @property
  def bbox_norm(self):
    """Length of bounding box diagonal"""
    import numpy as np
    return np.linalg.norm([x[1]-x[0] for x in self.bbox])

def res_type(resname):
  """Guess type from residue name"""
  if resname.startswith(('DA','DT','DC','DG')):
    this_type = "DNA"
  elif resname.startswith(('SOL','WAT','TIP','T3P')):
    this_type = "SOL"
  elif resname.startswith(('ION','MG','CL','NA','SOD','POT','CAL','LIT','Na+','Cl-')):
    this_type = "ION"
  elif resname.startswith(('DPPC','POPC','DAG','TAG')):
    this_type = "LIP"
  elif resname.startswith(('ATP','GTP')):
    this_type = "LIG"
  elif resname.startswith(("ALA","ARG","ASN","ASP",
                             "CYS","GLN","GLU","GLY","HSD",
                             "HSE","HIS","ILE","LEU","LYS",
                             "MET","PHE","PRO","SER","THR",
                             "TRP","TYR","VAL","HID","HIE","HIP","HSP")):
    this_type = "PRO"
  else:
    import mcol
    import mout
    mout.warningOut("Unknown residue type for "+mcol.arg+resname)
    this_type = None
  return this_type
