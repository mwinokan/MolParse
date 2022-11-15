
class System:

  def __init__(self,name):

    self.name = name
    self.description = None
    self.chains = []

    self.bondlist=None
    self.box = None

  def autoname_chains(self):
    import mout

    pro_count = 0
    lig_count = 0
    dna_count = 0
    sol_count = 0

    for i,chain in enumerate(self.chains):
      if chain.type == 'PRO':
        if lig_count > 4:
          mout.warningOut("Too many protein chains!")
          chain.name = 'P'
        else:
          chain.name = 'ABCDE'[pro_count]
          pro_count += 1
      elif chain.type == 'DNA':
        if lig_count > 4:
          mout.warningOut("Too many DNA chains!")
          chain.name = 'X'
        else:
          chain.name = 'XYZ'[dna_count]
          dna_count += 1
      elif chain.type == 'LIG':
        if lig_count > 1:
          mout.warningOut("Multiple ligand chains!")
        chain.name = 'L'
        lig_count += 1
      elif chain.type == 'SOL':
        if sol_count > 1:
          mout.warningOut("Multiple solvent chains!")
        chain.name = 'W'
        sol_count += 1
      elif chain.type == 'ION':
        chain.name = chain.residues[0].name[0]

      if chain.name in [c.name for c in self.chains[0:i]]:
        mout.warningOut(f"Ambiguous naming! Multiple chains named {chain.name}!")

  def check_indices(self):

    for index,atom in enumerate(self.atoms):
      if index != atom.index:
        print(index,atom.index,atom.name,atom.residue)

  def fix_indices(self):
    for index,atom in enumerate(self.atoms):
      atom.index = index

    for index,residue in enumerate(self.residues):
      residue.number = index

  def fix_atomnames(self,verbosity=1):
    import mout
    count=0
    for index,atom in enumerate(self.atoms):
      if atom.name[0].isnumeric():
        old = atom.name
        new = old[1:]+old[0]
        atom.set_name(new,verbosity=verbosity-1)
        if verbosity > 1:
          mout.out(old+" -> "+new)
        count+=1
    if verbosity > 0 and count != 0:
      mout.warningOut("Fixed "+str(count)+" atom names which appeared to have cycled.")

  def add_chain(self,chain):
    from .chain import Chain
    assert isinstance(chain,Chain)
    self.chains.append(chain)

  def add_system(self,system):

    for chain in system.chains:
      self.add_chain(chain)

    self.fix_indices()

  @property
  def num_atoms(self):
    num_atoms = 0
    for chain in self.chains:
      num_atoms += chain.num_atoms
    return num_atoms

  @property
  def charge(self):
    charge = 0
    for atom in self.atoms:
      charge += atom.charge
    return charge

  def summary(self,res_limit=10):
    import mout
    import mcol
    reset = mcol.clear+mcol.bold
    if self.description is not None:
      mout.headerOut(f'\n"{mcol.underline}{self.description}{reset}"')
    mout.headerOut("\nSystem "+mcol.arg+self.name+
                   mcol.clear+mcol.bold+" contains "+
                   mcol.result+str(self.num_chains)+
                   mcol.clear+mcol.bold+" chains:")
    for i,chain in enumerate(self.chains):
      mout.headerOut(f'Chain[{mcol.arg}{i}{reset}] {mcol.result}{chain.name}{reset} ({mcol.varType}{chain.type}{reset}) [#={mcol.result}{chain.num_residues}{reset}] =',end=' ')

      names = ""
      for name in chain.res_names[:res_limit]:
        names += name+" "
      if chain.num_residues > res_limit:
        names += "..."
      mout.out(names)
    # mout.varOut("Total Charge",self.charge)

  @property
  def num_chains(self):
    return len(self.chains)

  @property
  def chain_names(self):
    names = []
    for chain in self.chains:
      names.append(chain.name)
    return names

  def rename_atoms(self,old,new,res_filter=None,verbosity=2):
    import mcol
    import mout
    count=0
    for residue in self.residues:
      if res_filter is not None and res_filter not in residue.name:
        continue
      for atom in residue.atoms:
        if atom.name == old:
          count += 1
          atom.set_name(new,verbosity=verbosity-1)
    if verbosity > 0:
      if res_filter is None:
        mout.warningOut("Renamed "+mcol.result+str(count)+mcol.warning+" atoms from "+mcol.arg+old+
                        mcol.warning+" to "+mcol.arg+new)
      else:
        mout.warningOut("Renamed "+mcol.result+str(count)+mcol.warning+" atoms from "+mcol.arg+old+
                        mcol.warning+" to "+mcol.arg+new+mcol.warning+" with res_filter "+mcol.arg+res_filter)
    return count

  def rename_residues(self,old,new,verbosity=2):
    import mcol
    import mout
    count=0
    for residue in self.residues:
      if residue.name == old:
        # residue.name = new
        residue.rename(new,verbosity=verbosity-1)
        count += 1
    if verbosity > 0:
      mout.warningOut("Renamed "+mcol.result+str(count)+mcol.warning+" residues from "+mcol.arg+old+
                      mcol.warning+" to "+mcol.arg+new)
    return count

  def get_chain(self,name):
    import mout
    import mcol
    for chain in self.chains:
      if chain.name == name:
        return chain
    mout.errorOut("Chain with name "+mcol.arg+name+mcol.error+" not found.",fatal=True)

  def remove_chain(self,name,verbosity=1):
    import mcol
    import mout
    for index,chain in enumerate(self.chains):
      if chain.name == name:
        if verbosity > 0:
          mout.warningOut("Removing chain "+mcol.arg+name+str([index]))
        del self.chains[index]
        return chain
    mout.errorOut("Chain with name "+mcol.arg+name+mcol.error+" not found.",fatal=True)
    
  def remove_heterogens(self,verbosity=1):
    import mcol
    import mout
    del_list = []
    atoms = self.atoms
    for index,atom in enumerate(atoms):
      if atom.heterogen:
        del_list.append(index)
    number_deleted = self.remove_atoms_by_index(del_list,verbosity=verbosity-1)
    if verbosity > 0:
      mout.warningOut("Removed "+mcol.result+str(number_deleted)+mcol.warning+" heterogens")

  def remove_atoms_by_index(self,del_list,verbosity=1):
    import mcol
    import mout

    number_deleted=0
    self.fix_indices()
    for chain in self.chains:
      for residue in chain.residues:
        for index,atom in reversed(list(enumerate(residue.atoms))):
          if atom.index in del_list:
            del residue.atoms[index]
            number_deleted += 1
            if verbosity > 1:
              mout.warningOut("Removed atom "+
                              mcol.arg+atom.name+str([atom.index])+
                              mcol.warning+" in residue "+
                              mcol.arg+residue.name+str([residue.number]))
    return number_deleted

  def get_atom_by_index(self,index,pdb=True):
    for atom in self.atoms:
      if pdb:
        if atom.pdb_index == index:
          return atom
      else:
        if atom.index == index:
          return atom

  def remove_atoms(self,name,res_filter=None,verbosity=2):
    import mcol
    import mout

    number_deleted=0
    self.fix_indices()
    for chain in self.chains:
      for residue in chain.residues:
        if res_filter is not None and res_filter not in residue.name:
          continue
        for index,atom in reversed(list(enumerate(residue.atoms))):
          if atom.name == name:
            del residue.atoms[index]
            number_deleted += 1
            if verbosity > 1:
              mout.warningOut("Removed atom "+
                              mcol.arg+atom.name+str([atom.index])+
                              mcol.warning+" in residue "+
                              mcol.arg+residue.name+str([residue.number]))
    if verbosity > 0:
      if res_filter is None:
        mout.warningOut("Removed "+mcol.result+str(number_deleted)+mcol.warning+" atoms named "+mcol.arg+name)
      else:
        mout.warningOut("Removed "+mcol.result+str(number_deleted)+mcol.warning+" atoms named "+mcol.arg+name+mcol.warning+" with res_filter "+mcol.arg+res_filter)
    return number_deleted

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
  def masses(self):
    masses = []
    for chain in self.chains:
      masses += chain.masses
    return masses

  @property
  def atoms(self):
    atoms = []
    for chain in self.chains:
      atoms += chain.atoms
    return atoms

  def centre_of_mass(self,set=None,shift=None):
    return CoM(set=set,shift=shift)

  def CoM(self,set=None,shift=None):
    import mout
    import numpy as np

    position_list = self.positions

    centre_of_mass = np.array([sum([pos[0] for pos in position_list])/len(position_list),
                               sum([pos[1] for pos in position_list])/len(position_list),
                               sum([pos[2] for pos in position_list])/len(position_list)])

    mout.varOut("CoM of "+self.name,
                centre_of_mass,
                unit="Angstroms",precision=4)

    if set is not None:

      assert np.ndim(set) == 1

      try:
        new_com = np.array([set[0],set[1],set[2]])
      except:
        mout.errorOut("Incompatible array input",code="amp.system.CoM.1")

      shift_com = new_com - centre_of_mass

    else:

      shift_com = np.array([0,0,0])

    if shift is not None:

      assert np.ndim(shift) == 1

      try:
        new_shift = np.array([shift[0],shift[1],shift[2]])
      except:
        mout.errorOut("Incompatible array input",code="amp.system.CoM.2")

      shift_com = shift_com + new_shift
      
    if set is not None or shift is not None:

      for atom in self.atoms:
        atom.position = atom.position + shift_com

    return centre_of_mass

  @property
  def QM_indices(self):
    index_list = []
    for index,atom in enumerate(self.atoms):
      if atom.QM:
        index_list.append(index)
    return index_list

  @property
  def residues(self):
    residues = []
    for chain in self.chains:
      residues += chain.residues
    return residues

  @property
  def res_names(self):
    return [res.name for res in self.residues]

  @property
  def num_residues(self):
    return len(self.residues)

  @property
  def FF_atomtypes(self):
    atomtype_list = []
    for chain in self.chains:
      atomtype_list += chain.FF_atomtypes
    return atomtype_list

  @property
  def ase_atoms(self):
    from .io import write, read
    write("__temp__.pdb",self,verbosity=0)
    return read("__temp__.pdb",verbosity=0)

  def write_CJSON(self,filename,use_atom_types=False,gulp_names=False):
    from .io import writeCJSON
    writeCJSON(filename,self,use_atom_types=use_atom_types,gulp_names=gulp_names)

  def set_coordinates(self,reference):
    if type(reference) is str:
      from ase.io import read
      atoms = read(reference)
    elif isinstance(reference,list):
      for index,atom in enumerate(self.atoms):
        atom.position = reference[index]
      return
    else:
      atoms = reference

    for index,atom in enumerate(self.atoms):
      atom.position = atoms[index].position

  def rotate(self,angle,vector,center=(0,0,0)):
    ase_atoms = self.ase_atoms
    ase_atoms.rotate(angle,vector,center=center)
    self.set_coordinates(ase_atoms)

  def view(self):
    from .gui import view
    view(self)

  def auto_rotate(self):
    from .manipulate import auto_rotate
    ase_atoms = self.ase_atoms
    ase_atoms = auto_rotate(ase_atoms)
    self.set_coordinates(ase_atoms)

  def align_to(self,target):
    from ase.build import minimize_rotation_and_translation
    if isinstance(target,System):
      target = target.ase_atoms
    atoms = self.ase_atoms
    minimize_rotation_and_translation(target,atoms)
    self.set_coordinates(atoms)
    return atoms

  def copy(self,fast=False):
    if fast:
      new_sys = System(self.name)
      for chain in self.chains:
        new_sys.add_chain(chain.copy(fast=fast))
      new_sys.box = self.box
      return new_sys
    else:
      import copy
      return copy.deepcopy(self)

  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name
