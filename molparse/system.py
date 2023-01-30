
from .group import AtomGroup

class System(AtomGroup):

  """Top-level object for molecular systems

  These objects are rarely created by the user, 
  instead constructed automatically when parsing a 
  coordinate file via amp.parsePDB or otherwise"""

  def __init__(self,name: str):

    super(System, self).__init__(name)

    self.description = None
    self.chains = []

    self.bondlist=None
    self.box = None

  def autoname_chains(self,verbosity=1):
    """Automatically name chains"""
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
        if lig_count > 3:
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

      if verbosity > 0 and chain.name in [c.name for c in self.chains[0:i]]:
        mout.warningOut(f"Ambiguous naming! Multiple chains named {chain.name}!")

  def check_indices(self):
    """Print all child Atoms who's indices are incorrect"""

    for index,atom in enumerate(self.atoms):
      if index != atom.index:
        print(index,atom.index,atom.name,atom.residue)

  def fix_indices(self,verbosity=0):
    """Fix all child Atoms' indices"""
    if verbosity:
      import mout
      exclude = ['SOL','ION']
    for index,atom in enumerate(self.atoms):
      if verbosity == 2 and atom.index != index:
        mout.warningOut(f"Re-indexing atom {atom} (#{atom.index} --> #{index})")
      elif verbosity == 1 and atom.type not in exclude and atom.index != index:
        mout.warningOut(f"Re-indexing atom {atom} (#{atom.index} --> #{index})")
      atom.index = index

    for index,residue in enumerate(self.residues):
      if verbosity == 2 and residue.number != index:
        mout.warningOut(f"Re-indexing residue {residue} (#{residue.number} --> #{index})")
      elif verbosity == 2 and residue.type not in exclude and residue.number != index:
        mout.warningOut(f"Re-indexing residue {residue} (#{residue.number} --> #{index})")
      residue.number = index

      residue.fix_names()
      residue.fix_indices()

  def clear_pdbindices(self):
    """Clear all pdb indices from the system"""
    for atom in self.atoms:
      atom.pdb_index = None

  def fix_atomnames(self,verbosity=1):
    """Attempt to fix all child Atom names"""
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
    """Add a child Chain"""
    from .chain import Chain
    assert isinstance(chain,Chain)
    chain.index = len(self.chains)
    self.chains.append(chain)

  def add_system(self,system,same_chain=False):
    """Merge another system to this one"""

    if same_chain:
      for residue in system.residues:
        self.chains[-1].add_residue(residue.copy())

    else:
      for chain in system.chains:
        self.add_chain(chain.copy())

    self.fix_indices()

  @property
  def bbox(self):
    """Bounding box of the system"""
    import numpy as np
    x = [min([a.np_pos[0] for a in self.atoms]),max([a.np_pos[0] for a in self.atoms])]
    y = [min([a.np_pos[1] for a in self.atoms]),max([a.np_pos[1] for a in self.atoms])]
    z = [min([a.np_pos[2] for a in self.atoms]),max([a.np_pos[2] for a in self.atoms])]
    return [x,y,z]

  @property
  def bbox_sides(self):
    """Bounding box of the system"""
    import numpy as np
    x = max([a.np_pos[0] for a in self.atoms])-min([a.np_pos[0] for a in self.atoms])
    y = max([a.np_pos[1] for a in self.atoms])-min([a.np_pos[1] for a in self.atoms])
    z = max([a.np_pos[2] for a in self.atoms])-min([a.np_pos[2] for a in self.atoms])
    return [x,y,z]

  @property
  def bbox_norm(self):
    """Length of bounding box diagonal"""
    import numpy as np
    return np.linalg.norm([x[1]-x[0] for x in self.bbox])

  def check_intersection(self,system,radius=1,by_residue=True,boolean=False,chain=None):
    """Return a list of indices that are intersecting within a given radius between this system and another."""

    import mout
    # mout.debugOut(f"amp.System.check_intersection({system},radius={radius},by_residue={by_residue},boolean={boolean},chain={chain})")
    import numpy as np

    indices = []

    if by_residue:

      if chain:
        residues = self.get_chain(chain).residues
      else:
        residues = self.residues

      system_CoM = system.CoM(verbosity=0)
      r_system = system.bbox_norm
      num_residues = len(residues)

      # for each residue
      for i,res in enumerate(residues):
        if num_residues > 500 and i%100==0:
          mout.progress(i,num_residues,prepend="Calculating intersection",append=" of residues checked")

        residue_CoM = res.CoM(verbosity=0)
        d = np.linalg.norm(residue_CoM-system_CoM)
        r_residue = res.bbox_norm

        if d > r_system + r_residue:
          continue

        for atom2 in system.atoms:
          
          d = np.linalg.norm(residue_CoM-atom2.np_pos)

          if d > r_system:
            continue

          for atom1 in res.atoms:

            d = np.linalg.norm(atom1.np_pos - atom2.np_pos)

            if d <= radius:

              if boolean:
                return True
              else:
                indices.append(res.number)
                break

      if num_residues > 500:
        mout.progress(num_residues,num_residues,prepend="Calculating intersection",append=" of residues checked. Done.")

    else:

      if chain:
        atoms = self.get_chain(chain).atoms
      else:
        atoms = self.atoms

      for atom1 in self.atoms:
        for atom2 in system.atoms:
          if np.linalg.norm(atom1.np_pos - atom2.np_pos) <= radius:
            if boolean:
              return True
            else:
              indices.append(atom1.index)

    if boolean:
      return False
    else:
      return indices

  @property
  def num_atoms(self):
    """Number of child Atoms (int)"""
    num_atoms = 0
    for chain in self.chains:
      num_atoms += chain.num_atoms
    return num_atoms

  @property
  def charge(self):
    """Total charge (float)"""
    charge = 0
    for atom in self.atoms:
      charge += atom.charge
    return charge

  def summary(self,res_limit=10):
    """Print a summary of the System"""
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
      mout.headerOut(f'Chain[{mcol.arg}{i}{reset}] {mcol.result}{chain.name}{reset} ({mcol.varType}{chain.type}{reset}) [#={mcol.result}{chain.num_atoms}{reset}] =',end=' ')

      if chain.type == "PRO":
        names = ""
        for r in chain.residues[:res_limit]:
            names += r.letter + " "
        if chain.num_residues > res_limit:
          names += "..."
        mout.out(names)
      else:
        names = ""
        for name in chain.res_names[:res_limit]:
          names += name+" "
        if chain.num_residues > res_limit:
          names += "..."
        mout.out(names)
    # mout.varOut("Total Charge",self.charge)

  @property
  def num_chains(self):
    """Number of child Chains (int)"""
    return len(self.chains)

  @property
  def chain_names(self):
    """Get all Chain names (list)"""
    names = []
    for chain in self.chains:
      names.append(chain.name)
    return names

  def rename_atoms(self,old:str,new:str,res_filter:str=None,verbosity:int=2):
    """Rename all matching atoms"""
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

  def rename_residues(self,old:str,new:str,verbosity=2):
    """Rename all matching residues"""
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

  def get_chain(self,name:str):
    """Get Chain by name"""
    import mout
    import mcol
    for chain in self.chains:
      if chain.name == name:
        return chain
    mout.errorOut("Chain with name "+mcol.arg+name+mcol.error+" not found.",fatal=True)

  def remove_chain(self,name,verbosity=1):
    """Delete Chain by name"""
    import mcol
    import mout
    for index,chain in enumerate(self.chains):
      if chain.name == name:
        if verbosity > 0:
          mout.warningOut("Removing chain "+mcol.arg+name+str([index]))
        del self.chains[index]
        return chain
    mout.errorOut("Chain with name "+mcol.arg+name+mcol.error+" not found.",fatal=True)
    
  def remove_heterogens(self,verbosity:int=1):
    """Remove HETATM entries"""
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

  def remove_atoms_by_index(self,del_list:list,verbosity:int=1):
    """Remove Atoms by their index"""
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

  def remove_residues_by_index(self,del_list:list,verbosity:int=1):
    """Remove Atoms by their index"""
    import mcol
    import mout

    number_deleted=0
    self.fix_indices()
    for chain in self.chains:
      for index,residue in reversed(list(enumerate(chain.residues))):
          if residue.number in del_list:
            del chain.residues[index]
            number_deleted += 1
            if verbosity > 1:
              mout.warningOut("Removed residue "+
                              mcol.arg+residue.name+str([residue.number])+
                              mcol.warning+" in chain "+
                              mcol.arg+chain.name)
    return number_deleted

  def get_atom_by_index(self,index:int,use_pdb_index:bool=True):
    """Get Atom by its index"""
    for atom in self.atoms:
      if use_pdb_index:
        if atom.pdb_index == index:
          return atom
      else:
        if atom.index == index:
          return atom

  def remove_atoms(self,name:str,res_filter:str=None,verbosity:int=2):
    """Remove Atoms by their name"""
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
    """Get all child Atom names (list)"""
    names_list = []
    for chain in self.chains:
      names_list += chain.atom_names(wRes=wRes,noPrime=noPrime)
    return names_list

  @property
  def atomic_numbers(self):
    """Get all child Atom atomic numbers (list)"""
    number_list = []
    for chain in self.chains:
      number_list.append(chain.atomic_numbers)
    return number_list

  @property
  def positions(self):
    """Get all child Atom positions (list)"""
    positions_list = []
    for chain in self.chains:
      positions_list += chain.positions
    return positions_list

  @property
  def charges(self):
    """Get all child Atom charges (list)"""
    charges = []
    for chain in self.chains:
      charges += chain.charges
    return charges

  @property
  def symbols(self):
    """Get all child Atom symbols (list)"""
    symbols = []
    for atom in self.atoms:
      symbols.append(atom.symbol)  
    return symbols

  @property
  def masses(self):
    """Get all child Atom masses (list)"""
    masses = []
    for chain in self.chains:
      masses += chain.masses
    return masses

  @property
  def atoms(self):
    """Get all child Atoms (list)"""
    atoms = []
    for chain in self.chains:
      atoms += chain.atoms
    return atoms

  def centre_of_mass(self,set=None,shift=None):
    """Calculate centre of mass"""
    return self.CoM(set=set,shift=shift)

  def center(self):
    """Move the system's CoM to the origin"""
    return self.CoM(set=[0,0,0])

  def CoM(self,set=None,shift=None,verbosity=1):
    """Calculate or manipulate the system's centre of mass. 

    if not set and not shift: return CoM
    if set: move the System to the CoM
    if shift: move the System by the specified vector"""

    import mout
    import numpy as np

    position_list = self.positions

    centre_of_mass = np.array([sum([pos[0] for pos in position_list])/len(position_list),
                               sum([pos[1] for pos in position_list])/len(position_list),
                               sum([pos[2] for pos in position_list])/len(position_list)])

    if verbosity > 0: 
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
    """Return list of Atom indices with the QM flag"""
    index_list = []
    for index,atom in enumerate(self.atoms):
      if atom.QM:
        index_list.append(index)
    return index_list

  @property
  def residues(self):
    """Get all child Residues (list)"""
    residues = []
    for chain in self.chains:
      residues += chain.residues
    return residues

  @property
  def res_names(self):
    """Get all child Residue names (list)"""
    return [res.name for res in self.residues]

  @property
  def num_residues(self):
    """Get number of child Residue (int)"""
    return len(self.residues)

  @property
  def FF_atomtypes(self):
    """Get all child Atom atomtypes (list)"""
    atomtype_list = []
    for chain in self.chains:
      atomtype_list += chain.FF_atomtypes
    return atomtype_list

  @property
  def ase_atoms(self):
    """Construct an equivalent ase.Atoms object"""
    
    # from .io import write, read
    # write("__temp__.pdb",self,verbosity=0)
    # return read("__temp__.pdb",verbosity=0)
    from ase import Atoms
    return Atoms(symbols=self.symbols,cell=None, pbc=None,positions=self.positions)

  def write_CJSON(self,filename,use_atom_types=False,gulp_names=False):
    """Export a CJSON"""
    from .io import writeCJSON
    writeCJSON(filename,self,use_atom_types=use_atom_types,gulp_names=gulp_names)

  def set_coordinates(self,reference):
    """Set all coordinates according to a reference ase.Atoms object"""
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
    """Rotate the system (see ase.Atoms.rotate)"""
    ase_atoms = self.ase_atoms
    ase_atoms.rotate(angle,vector,center=center)
    self.set_coordinates(ase_atoms)

  def view(self,**kwargs):
    """View the system with ASE"""
    from .gui import view
    view(self,**kwargs)

  def plot(self,ax=None,color=None,center_index=0,show=False,offset=None,padding=1,zeroaxis=None,frame=False,labels=True,textdict={"horizontalalignment":"center","verticalalignment":"center"}):
    """Render the system with matplotlib"""
    
    import numpy as np
    from ase.visualize.plot import plot_atoms

    # copy the system so as not to alter it
    copy = self.copy()

    # create new axes if none provided
    if ax is None:
      import matplotlib.pyplot as plt
      fig,ax = plt.subplots()

    # create colours list
    if color is not None:
      color = [color for a in copy.atoms]

    # create the canvas unit size from bounding box
    canvas = np.array(self.bbox_sides)*1.5

    # center the render if needed
    if offset is None:
      offset=[0,0]
      if center_index is not None:
        vec = - copy.atoms[center_index].np_pos + canvas
        copy.CoM(shift=vec,verbosity=0)

    # use ASE to render the atoms
    plot_atoms(copy.ase_atoms,ax,colors=color,offset=offset,bbox=[0,0,canvas[0]*2,canvas[1]*2])

    # do the labels
    if labels:
      for atom in copy.atoms:
        if atom.species != 'H':
          ax.text(atom.position[0],atom.position[1],atom.name,**textdict)

    # crop the plot
    xmax = copy.bbox[0][1] + padding
    xmin = copy.bbox[0][0] - padding
    ymax = copy.bbox[1][1] + padding
    ymin = copy.bbox[1][0] - padding
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # render the centering axes
    if zeroaxis is not None:
      ax.axvline(canvas[0],color=zeroaxis)
      ax.axhline(canvas[1],color=zeroaxis)

    if not frame:
      ax.axis('off')

    if show:
      plt.show()

    return ax, copy

  def plot3d(self,extra=[],alpha=1.0):
    """Render the system with plotly graph objects. 
    extra can contain pairs of coordinates to be shown as vectors."""
    from .go import plot3d
    return plot3d(self.atoms,extra,alpha)

  def auto_rotate(self):
    """Rotate the system into the XY plane"""
    from .manipulate import auto_rotate
    ase_atoms = self.ase_atoms
    ase_atoms = auto_rotate(ase_atoms)
    self.set_coordinates(ase_atoms)

  def align_to(self,target):
    """Minimise rotation w.r.t. a target"""
    from ase.build import minimize_rotation_and_translation
    if isinstance(target,System):
      target = target.ase_atoms
    atoms = self.ase_atoms
    minimize_rotation_and_translation(target,atoms)
    self.set_coordinates(atoms)
    return atoms

  def translate(self,displacement):
    """Translate the system"""
    self.CoM(shift=displacement,verbosity=0)

  def align_by_posmap(self,map):

    """Align the system to a target by superimposing three shared atoms:

      a --> A
      b --> B
      c --> C

      (a,b,c) are atoms from this system
      (A,B,C) are atoms from the target system

      map should contain Atoms: [[a,b,c],[A,B,C]]

    """

    import numpy as np

    a = map[0][0]
    b = map[0][1]
    c = map[0][2]

    A = map[1][0]
    B = map[1][1]
    C = map[1][2]

    # TRANSLATION

    displacement = A - a
    self.translate(displacement)

    # ROTATION 1

    d = (b+c)/2
    D = (B+C)/2
    self.rotate(d-A.np_pos,D-A.np_pos,center=A.np_pos)

    # ROTATION 2

    d = (b+c)/2
    D = (B+C)/2

    v_bc = c - b
    v_BC = C - B

    def unit_vector(a):
      return a/np.linalg.norm(a)

    def remove_parallel_component(vector,reference):
      unit_reference = unit_vector(reference)
      projection = np.dot(vector,unit_reference)
      return vector - projection

    v_ad = d - a.np_pos
    v_bc_normal_to_ad = remove_parallel_component(v_bc, v_ad)
    v_BC_normal_to_ad = remove_parallel_component(v_BC, v_ad)
    self.rotate(v_bc_normal_to_ad, v_BC_normal_to_ad,center=A.np_pos)

    # extra stuff for plotly
    extra = [[A.np_pos,A.np_pos+v_ad]]
    return extra

  def align_by_pairs(self,target,index_pairs,alt=False):
    """Align the system (i) to the target (j) by consider the vectors:

        a --> b
        a --> c

        where index pairs contains the indices for the three atoms: 
        a,b,c in the respective systems:

        index_pairs = [[i_a,j_a],[i_b,j_b],[i_c,j_c]]

        Alternatively you can pass the positions j_a, j_b, j_c as target,
        and index_pairs can be i_a, i_b, i_c.
        """

    assert len(index_pairs) == 3
    import numpy as np

    if isinstance(target,System):
      for pair in index_pairs:
        assert self.atoms[pair[0]].name == target.atoms[pair[1]].name

      pos_0_0 = self.atoms[index_pairs[0][0]].np_pos
      pos_1_0 = self.atoms[index_pairs[1][0]].np_pos
      pos_2_0 = self.atoms[index_pairs[2][0]].np_pos

      pos_0_1 = target.atoms[index_pairs[0][1]].np_pos
      pos_1_1 = target.atoms[index_pairs[1][1]].np_pos
      pos_2_1 = target.atoms[index_pairs[2][1]].np_pos

    else:

      pos_0_0 = self.atoms[index_pairs[0]].np_pos
      pos_1_0 = self.atoms[index_pairs[1]].np_pos
      pos_2_0 = self.atoms[index_pairs[2]].np_pos

      pos_0_1 = target[0]
      pos_1_1 = target[1]
      pos_2_1 = target[2]

    vec = pos_0_1 - pos_0_0
    self.CoM(shift=vec,verbosity=0)

    a = pos_1_0 - pos_0_0
    b = pos_1_1 - pos_0_1
    self.rotate(a,b,pos_0_1)

    a = pos_1_0 - pos_0_0
    a_hat = a/np.linalg.norm(a)

    b = pos_2_0 - pos_0_0
    c = pos_2_1 - pos_0_1
    d = b - np.dot(a,b)*a_hat
    e = c - np.dot(a,c)*a_hat

    d_hat = d/np.linalg.norm(d)
    e_hat = e/np.linalg.norm(e)
    ang = np.arccos(np.clip(np.dot(d_hat, e_hat), -1.0, 1.0))

    if alt:
      self.rotate((ang/np.pi*180),a,pos_0_1)
    else:
      self.rotate(-(90-ang/np.pi*180),a,pos_0_1)

  def guess_names(self,target):
    """Try and set the atom names of the system by looking for 
    the closest atom of the same species in the target system"""

    import numpy as np

    positions = [b.np_pos for b in target.atoms]
    species = [b.species for b in target.atoms]

    for a in self.atoms:
      a_pos = a.np_pos
      distances = [np.linalg.norm(a_pos - p) if s == a.species else 999 for s,p in zip(species,positions)]
      index = np.argmin(distances)
      b = target.atoms[index]
      a.set_name(b.name,verbosity=0)

  def rmsd(self,reference):
    """Calculate the RMS displacement between this system and a reference"""
    
    import numpy as np
    assert len(self.atoms) == len(reference.atoms)
    displacements = np.array([np.linalg.norm(a.np_pos-b.np_pos) for a,b in zip(self.atoms,reference.atoms)])
    return np.sqrt(np.mean(displacements**2))

  def reorder_atoms(self,reference,map:dict=None):
    new_sys = self.copy()
    for atom in reference.atoms:
      res = self.get_residue(atom.residue,map=map)
      new_sys.add_atom(res.get_atom(atom.name).copy())
    new_sys.fix_indices()
    self = new_sys

  def get_residue(self,name:str,map:dict=None):
    """ return residues with matching name"""
    import mout
    if map is not None:
      if name in map.keys():
        name = map[name]
    matches = [r for r in self.residues if r.name == name]
    if len(matches) == 0:
      mout.errorOut(f"No residue found with name {name}")
      return []
    elif len(matches) == 1:
      return matches[0]
    else:
      mout.warning(f"Multiple residues found with name {name}")
      return matches

  def copy(self,disk=False,alt=False):
    """Return a deepcopy of the System"""
    if disk:
      from .io import write, parseGRO
      write(f"__temp__.gro",self,verbosity=0)
      system = parseGRO(f"__temp__.gro",verbosity=0)
      system.name = self.name
      return system
    elif alt:
      if self.num_chains != len(set([str(c) for c in self.chains])):
        import mout
        mout.errorOut("System has duplicate chain names so this copying method will have incorrect chains!")
      copy_system = System(self.name + " (copy)")
      copy_system.box = self.box
      for atom in self.atoms:
        copy_system.add_atom(atom)
      return copy_system
    else:
      import copy
      return copy.deepcopy(self)

  def subsystem(self,indices,use_pdb_index=True):
    """Extract a subset of the system by the atom indices"""

    new_system = System(self.name+" (filtered)")
    new_system.box = self.box

    for index in indices:

      atom = self.get_atom_by_index(index,use_pdb_index=use_pdb_index)

      # print(index,atom)
      if atom is None:
        continue

      new_system.add_atom(atom)

    return new_system

  def add_atom(self,atom):
    """add an atom to the system"""
    from .chain import Chain
    from .residue import Residue, res_type

    if self.chains and atom.chain == self.chains[-1].name and res_type(atom.residue) == self.chains[-1].residues[-1].type:
      self.chains[-1].add_atom(atom)
    else:
      chain = Chain(atom.chain)
      chain.add_atom(atom)
      self.add_chain(chain)

  @property
  def children(self):
    return self.chains