
from .group import AtomGroup

# from enum import Enum
# class UnknownSiteHandlingMethod(Enum):
#   IGNORE = 0 # atoms with no alternative sites are ignored when separating
#   UNION = 1  # all atoms with no alternative sites are added to each separated residue     

class Residue(AtomGroup):
  """Class for chemical Residue

  These objects should not be created by the user, 
  but constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""

  def __init__(self,name: str,index: int=None, number: int=None,chain: str=None,atoms=None):
    
    super(Residue, self).__init__(name)

    self.chain = chain
    self.index = index
    self._number = number

    from .list import NamedList
    self._atoms = NamedList()
    
    if atoms:

      from .atom import Atom

      for atom in atoms:
        self.add_atom(
          Atom(
            name=atom.name,
            index=atom.index,
            pdb_index=atom.number,
            position=atom.position,
            residue=self.name,
            res_number=self.number,
            chain=self.chain,
        ))

    self._parent = None

  @property
  def alternative_sites(self):
    if self.contains_alternative_sites:

      if self.name == 'LIG':
        import mout
        atoms_w_alt_site = [a.alternative_site for a in self.atoms if a.alternative_site is not None]
        # mout.debug(f'{atoms_w_alt_site=}')
        # mout.debug(f'{list(set(atoms_w_alt_site))=}')

      return list(set([a.alternative_site for a in self.atoms if a.alternative_site is not None]))

    return []

  def split_by_site(self,unknown_site_handling='ignore'):
    """return a list of Residue objects, one for each site"""
    if self.contains_alternative_sites:
      ignore_count = 0
      separated_residues = []
      import mout
      # mout.debug(f'{self.alternative_sites=}')
      for site in self.alternative_sites:
        
        res = Residue(self.name, index=self.index, number=self.number, chain=self.chain)
        for atom in self.atoms:
          if atom.alternative_site == site:
            res.add_atom(atom)
          elif atom.alternative_site is None:
            if unknown_site_handling == 'union':
              res.add_atom(atom)  
            else:
              ignore_count += 1

        separated_residues.append(res)

      if ignore_count:
        import mout
        mout.warning(f'Ignored {ignore_count} atoms with no alternative site [{self.name_number_str}]')

      return separated_residues
    else:
      return [self]

  @property
  def name_number_str(self):
    return f'{self.name} {self.number}'

  @property
  def name_number_chain_str(self):
    return f'{self.name} {self.number} {self.chain}'

  @property
  def parent(self):
    return self._parent

  @parent.setter
  def parent(self,obj):
    # import weakref
    # self._parent = weakref.ref(obj)
    self._parent = obj

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
  def resid(self):
    """resid and number are the residue index from the original file"""
    return self._number
  
  @property
  def number(self):
    """resid and number are the residue index from the original file"""
    return self._number

  @property
  def ase_atoms(self):
    from ase import Atoms
    atoms = self.atoms
    symbols = [a.symbol for a in atoms]
    positions = [a.position for a in atoms]
    return Atoms(symbols=symbols,cell=None, pbc=None,positions=positions)

  def view(self):
    """View the system with ASE"""
    from .gui import view
    view(self.ase_atoms)

  def fix_names(self):
    """Ensure child Atoms have correct parent name"""
    for atom in self.atoms:
      atom.residue = self.name

  def fix_indices(self):
    """Ensure child Atoms have correct res_index"""
    for atom in self.atoms:
      atom.res_index = self.index

  def atom_names(self,wRes: bool=False,noPrime: bool=False):
    """Returns names of all child Atoms (list)"""
    names_list = []
    for atom in self._atoms:
      names_list.append(atom.get_name(wRes=wRes,noPrime=noPrime))
    return names_list

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

  def addAtom(self,atom,copy=True):
    """add an Atom"""
    from .atom import Atom
    assert isinstance(atom,Atom)

    if copy:
      atom_new = atom.copy()
      # import mout
      # mout.debug(f'copied atom: {atom} {id(atom)} --> {id(atom_new)}')
      atom = atom_new

    atom.chain = self.chain
    atom.residue = self.name
    atom.res_index = self.index
    atom._res_number = self.number

    atom.parent = self

    if self.atoms:
      # remove any termini
      self.atoms[-1].terminal = None
      self.atoms[-1].ter_line = None

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
  def atoms(self):
    """Returns all child Atoms (list)"""
    return self._atoms

  def get_atom(self,name,verbosity=1):
    """Get a child Atom by name"""
    if isinstance(name,list):
      return [self.get_atom(n,verbosity=verbosity-1) for n in name]

    for atom in self._atoms:
      if atom.name == name: return atom

    import mout
    import mcol
    if verbosity > 0:
      mout.errorOut("No atom "+
                    mcol.arg+name+
                    mcol.error+" in residue "+
                    mcol.arg+self.name+str([self.number]),fatal=False,code="Residue.1")
    return None

  def set_positions(self,positions):
    """Set positions for all child atoms"""
    for atom,pos in zip(self.atoms,positions):
      atom.position = pos

  def align_to(self,target,names):
    """Minimise rotation w.r.t. a target"""

    import numpy as np
    from ase.build.rotate import rotation_matrix_from_points

    assert isinstance(target,Residue)

    # get subset of the residue to align
    self_copy = self.copy()
    self_copy.delete_atom([n for n in self_copy.atom_names() if n not in names])
    
    # get subset of the target
    target_copy = target.copy()
    target_copy.delete_atom([n for n in target_copy.atom_names() if n not in names])

    # position lists
    p = self_copy.ase_atoms.get_positions()
    p0 = target_copy.ase_atoms.get_positions()

    # centeroids
    c = np.mean(p, axis=0)
    c0 = np.mean(p0, axis=0)
        
    # subtract centroids
    p -= c
    p0 -= c0

    # Compute rotation matrix
    R = rotation_matrix_from_points(p.T, p0.T)

    # compute displacement
    p_new = np.dot(self_copy.ase_atoms.get_positions(), R.T) + c0
    d = p0[0] + c0 - p_new[0]

    # set the new positions
    self.set_positions(np.dot(self.ase_atoms.get_positions(), R.T) + c0 + d)

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

  def summary(self):
    """Summarised output of the Residue"""
    import mout
    import mcol
    mout.headerOut(f'\nResidue {self.name}, index={self.index}, number={self.number}, chain={self.chain}, #atoms={self.num_atoms}')

    mout.out(f"{mcol.underline}{'Sy':2} {'NAME':4} {'INDEX':>6} {'NUMBER':>6} {'X':>7} {'Y':>7} {'Z':>7} {'Alt.':>4}{mcol.clear} ")
    for atom in self.atoms:
      mout.out(f'{atom.symbol:2} {atom.name:4} {atom.index if atom.index is not None else " ":>6} {atom.number if atom.number is not None else " ":>6} {atom.x:>7.2f} {atom.y:>7.2f} {atom.z:>7.2f} {atom.alternative_site or " ":>4}')
    
  def is_same_as(self,residue):
    assert isinstance(residue, Residue)

    if self.name != residue.name:
      return False
    if self.number != residue.number:
      return False
    if self.chain != residue.chain:
      return False
    return True

  def set_chain_number(self,index):
    for atom in self.atoms:
      atom.chain_number = index

  def set_chain_char(self,char):
    for atom in self.atoms:
      atom.chain = char

  @property
  def children(self):
    return self.atoms

RES_TYPES = {
  'DNA': ('DA','DT','DC','DG','ADE9','THMN','GUA9','CTSN','DOG'),
  'SOL': ('SOL','WAT','TIP','T3P','HOH','PEG','SO4','DMS','H2S'),
  'ION': ('ION','MG','CL','NA','SOD','POT','CAL','LIT','Na+','Cl-','CA','ZN'),
  'LIP': ('DPPC','POPC','DAG','TAG'),
  'LIG': ('ATP','GTP','LIG','UNL','CTP','TTP','OGTP'),
  'N/A': ('QM','MM'),
  'PRO': ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HSD","HSE","HIS",
          "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
          "HID","HIE","HIP","HSP"),
}

def res_type(resname):
  """Guess type from residue name"""

  for k,v in RES_TYPES.items():
    if resname in v:
      return k

  # Amber style amino acid termini
  if resname[0] in ("N","C") and len(resname) == 4 and resname[1:] in RES_TYPES['PRO']:
    return 'PRO'
  
  # Chromophore
  if resname == 'CRO':
    return 'PRO'

  import mcol
  import mout
  mout.warningOut("Unknown residue type for "+mcol.arg+resname)
  return None
