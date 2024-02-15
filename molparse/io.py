
def write(filename,image,verbosity=1,**parameters):
  from ase import io
  from ase import atoms as aseatoms
  import mcol
  import mout
  import json
  from .system import System
  from .group import AtomGroup
  import plotly.graph_objects as go
  import pickle

  filename = str(filename)

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             end='') # user output

  # Different behaviour depending on format:
  try:
    # pickle:
    if filename.endswith(".pickle"):
      with open(filename,'wb') as f:
        pickle.dump(image,f)

    # JSON:
    elif filename.endswith(".json"):
      with open(filename,'wt') as f:
        json.dump(image,f,indent=2)

    # CJSON:
    elif filename.endswith(".cjson"):
      writeCJSON(filename,image)
    # PBD and amp.System:
    elif (filename.endswith(".pdb") or filename.endswith(".pdb2")) and isinstance(image,System):
      writePDB(filename,image,verbosity=verbosity-1)
    elif (filename.endswith(".pdb") or filename.endswith(".pdb2")) and isinstance(image,AtomGroup):
      writePDB(filename,image,verbosity=verbosity-1)
    # GRO and amp.System:
    elif filename.endswith(".gro") and isinstance(image,System):
      writeGRO(filename,image,verbosity=verbosity-1)
    # PDB and list of amp.System's:
    elif filename.endswith(".pdb") and isinstance(image,list):

      ### List of Systems
      if all([isinstance(i,System) for i in image ]):

        for index,frame in enumerate(image):
          if index == 0:
            writePDB(filename,frame,verbosity=verbosity-1)
          else:
            writePDB(filename,frame,verbosity=verbosity-1,append=True,model=index+1)

      elif all([isinstance(i,aseatoms.Atoms) for i in image ]):
        io.write(filename,image,**parameters)
      else:
        mout.errorOut("Unsupported",fatal=True)

    elif isinstance(image,System):
      mout.errorOut(f"Filetype {mcol.file}{filename.split('.')[-1]}{mcol.error} not supported for mp.System object",code="amp.io.write[1]")
      return None

    elif isinstance(image,go.Figure):
      if filename.endswith('.html'):
        image.write_html(filename)
      else:
        if filename.endswith('.pdf'):
          import plotly.io as pio   
          if pio.kaleido.scope is not None:
            pio.kaleido.scope.mathjax = None
        image.write_image(filename)

    # Others:
    else:
      io.write(filename,image,**parameters)
  except TypeError as e:
    mout.out("Fail.")
    mout.error(e)
    mout.error(f"{type(image)} could not be written to "+mcol.file+filename,code="amp.io.write[1]")
    return None

  if (verbosity > 0):
    mout.out("Done.") # user output

def read(filename,index=None,verbosity=1,tagging=True,tagByResidue=False,**parameters):
  from ase import io
  import mcol
  import mout

  if (verbosity > 0):
    mout.out("reading "+mcol.file+
             filename+
             mcol.clear+" ... ",
             end='') # user output

  try:
    atoms = io.read(filename,index,**parameters)
  except FileNotFoundError:
    mout.out("Fail.")
    mout.errorOut("Could not read missing file "+mcol.file+filename,code="amp.io.read[1]",end="")
    return None

  if tagging and filename.endswith(".pdb"):
    getAndSetTags(filename,atoms,byResidue=tagByResidue)

  if (verbosity > 0):
    mout.out("Done.") # user output

  return atoms

def getAndSetTags(pdb,atoms,byResidue=False):
  from ase import atoms as aseatoms
  import mout

  # Parse the first model in the PDB to get the tags
  taglist=[]
  with open(pdb,"r") as input_pdb:
    searching = True
    for line in input_pdb:
      if searching:
        if line.startswith("MODEL"):
          searching = False
        elif line.startswith("ATOM"):
          searching = False
          taglist.append(tagFromLine(line,byResidue=byResidue))
      else:
        if line.startswith("ENDMDL"):
          break
        if line.startswith("END"):
          break
        if line.startswith("TER"):
          continue
        else:
          # get the tags:
          taglist.append(tagFromLine(line,byResidue=byResidue))

  # set the tags
  if isinstance(atoms,aseatoms.Atoms):
    for index,tag in enumerate(taglist):
      atoms[index].tag = tag
  else:
    for atms in atoms:
      for index,tag in enumerate(taglist):
        atms[index].tag = tag

def tagFromLine(line,byResidue):
  try:
    if byResidue:
      return int(line[24:27])
    else:
      return int(''.join(filter(lambda i: i.isdigit(), line.strip().split()[2])))
  except:
    return 0

def parse(file,verbosity=1):
  from .ndx import parseNDX
  file = str(file)
  extension = file.split(".")[-1]
  if extension == "pdb":
    return parsePDB(file,verbosity=verbosity)
  elif extension == "gro":
    return parseGRO(file,verbosity=verbosity)
  elif extension == "xyz":
    return parseXYZ(file,verbosity=verbosity)
  elif extension == "mol":
    return parseMol(file,verbosity=verbosity)
  elif extension == "ndx":
    return parseNDX(file,verbosity=verbosity)
  else:
    import mout
    mout.errorOut("Unsupported file type for MolParse parsing, using ASE.io.read")
    return read(file,verbosity=verbosity)

def parseMol(mol, verbosity=1):
    
    mol = str(mol)
    assert mol.endswith(".mol"), f"filename doesn't end with '.mol': {mol}"

    if verbosity:
      import mout
      import mcol
      mout.out("parsing "+mcol.file+
               mol+
               mcol.clear+" ... ",
               end='') # user output

    from rdkit.Chem import MolFromMolFile
    mol = MolFromMolFile(mol)

    if verbosity:
        mout.out("Done.") # user output

    return mol

def parsePDB(pdb,
  systemName=None,
  index=1,
  fix_indices=True,
  fix_atomnames=True,
  autoname_chains=False,
  prune_alternative_sites=False,
  verbosity=1,
  debug=False,
  dry=False,
  resfilter=None,
  alternative_site_warnings=False,
  split_chains_by_type=True,
  remove_empty_chains=True,
  num_appended_hydrogen_chains=False,
  detect_misordered_chains=True,
  reordering_warnings=False,
  reordering_summary_warnings=True,
  skipping_residue_warning=True,
  keep_headers=False,
):

  if resfilter is None:
    resfilter = []
  
  try:

    from .system import System
    from .chain import Chain
    from .residue import Residue, RES_TYPES
    import mout
    import mcol

    pdb = str(pdb)
    assert pdb.endswith(".pdb"), f"filename doesn't end with '.pdb': {pdb}"

    try:
      index = int(index)
    except:
      if index == ":":
        if verbosity > 0:
          mout.warningOut("Parsing all models in "+mcol.file+pdb)
      else:
        mout.errorOut("Unsupported index: '"+str(index)+"'",fatal=True)

    if index == ":":

      all_systems = []

      import subprocess
      
      last_model_line = subprocess.check_output("grep MODEL "+pdb+" | tail -n1 ", shell=True)
      
      for i in range(1,int(last_model_line.split()[-1])+1):

        if (verbosity > 0):
          mout.out("\rparsing "+mcol.file+
                   pdb+
                   mcol.clear+" (model "+str(i)+") ... ",
                   end='') # user output

        all_systems.append(parsePDB(pdb,
                                    systemName=systemName,
                                    index=i,
                                    fix_indices=fix_indices,
                                    fix_atomnames=fix_atomnames,
                                    verbosity=verbosity-1,
                                    debug=debug,dry=dry))

      if (verbosity > 0):
        mout.out("\rparsing "+mcol.file+
                   pdb+
                   mcol.clear+" (models 1-"+str(i)+") ... Done.") # user output

      return all_systems

    else:

      if (verbosity > 0):
        mout.out("parsing "+mcol.file+
                 pdb+
                 mcol.clear+" ... ",
                 end='') # user output

      if debug:
        mout.warningOut("Debug mode active!")

      import os

      residue = None
      chain = None
      last_residue_name = None
      last_residue_number = None
      last_chain_name = None
      last_residue_type = None
      res_counter = 0
      chain_counter = 0
      atom_counter = 1
      was_terminal = False

      if systemName is None:
        systemName = os.path.splitext(pdb)[0]
      system = System(name = systemName)

      try: 
        with open(pdb,"r") as input_pdb:
          header_data = []
          searching = True
          for line in input_pdb:

            if keep_headers:
              if line.startswith('TITLE'):
                header_data.append(line)
              elif line.startswith('HEADER'):
                header_data.append(line)
              elif line.startswith('COMPND'):
                header_data.append(line)
              elif line.startswith('REMARK'):
                header_data.append(line)
              elif line.startswith('CRYST'):
                header_data.append(line)
              elif line.startswith('SCALE'):
                header_data.append(line)

            else:

              if line.startswith('TITLE'):
                system.title = line.removeprefix('TITLE ').strip()
                # mout.debug(f'{system.title=}')
                continue
              elif line.startswith('HEADER'):
                system.header = line.removeprefix('HEADER').strip()
                # mout.debug(f'{system.header=}')
                continue

            if searching:
              if line.startswith("MODEL"):

                if index == 1:
                  searching = False
                elif " "+str(index) in line:
                  searching = False
          
              elif (index == 1 and line[0:6] in ['ATOM  ','HETATM']) or not searching:
                searching = False

                if debug:
                  mout.debug(line)

                #### PARSELINE
                try:
                  atom = parsePDBAtomLine(line,res_counter,atom_counter,chain_counter,debug=debug)
                except Exception as e:
                  mout.error(f'{pdb=} {index=}')
                  raise Exception(e)

                chain = Chain(atom.chain)
                residue = new_residue(atom.residue,res_counter,atom.res_number,atom.chain)
                residue.addAtom(atom,copy=False)
                last_residue_name = atom.residue
                last_residue_number = atom.res_number
                last_chain_name = atom.chain
            else:
              if line.startswith("ENDMDL") and not searching:
                break
              elif line.startswith("END") and not searching:
                break
              elif line.startswith("# All scores below are weighted scores, not raw scores.") and not searching:
                break
              elif line.startswith("TER"):
                if verbosity > 1:
                  mout.warningOut("Terminal added to "+mcol.arg+chain.name+":"+residue.name)
                residue.atoms[-1].terminal = True
                was_terminal = True
                continue
              elif line.startswith("CONECT"):
                continue
              elif line.startswith("OMPND"):
                continue
              elif line.startswith("COMPND"):
                continue
              elif line.startswith("AUTHO"):
                continue
              elif line.startswith("MASTER"):
                continue
              elif line.startswith("ANISOU"):
                continue
              elif dry and any(res in line for res in RES_TYPES['SOL']):
                continue
              else:

                if debug:
                  mout.debug(line)

                ### PARSELINE
                try:
                  atom = parsePDBAtomLine(line,res_counter,atom_counter,chain_counter,debug=debug,alternative_site_warnings=alternative_site_warnings)
                except Exception as e:
                  mout.error(f'{pdb=} {index=}')
                  raise Exception(e)

                make_new_res = False
                if residue is None: make_new_res = True
                if last_residue_name != atom.residue: make_new_res = True
                if last_residue_number != atom.res_number: make_new_res = True

                make_new_chain = False
                if residue is None: make_new_chain = True
                if last_chain_name != atom.chain: make_new_chain = True
                if was_terminal: make_new_chain = True

                if make_new_res:
                  if residue is not None: 
                    if not resfilter or residue.name in resfilter:
                      chain.add_residue(residue)
                      res_counter = res_counter+1
                    elif skipping_residue_warning:
                      mout.warning(f'Skipping residue {residue} {residue.number}')
                  residue = new_residue(atom.residue,res_counter,atom.res_number,atom.chain)
                  if split_chains_by_type and last_residue_type != residue.type: 
                    make_new_chain = True
                  if make_new_chain:
                    if chain is not None:
                      system.add_chain(chain)
                      chain_counter += 1
                    chain = Chain(atom.chain)

                residue.addAtom(atom,copy=False)

                atom_counter += 1

                last_residue_type = residue.type
                last_residue_number = atom.res_number
                last_residue_name = atom.residue
                last_chain_name = atom.chain
                was_terminal = False

      except FileNotFoundError:
        mout.errorOut("File "+mcol.file+pdb+mcol.error+" not found",fatal=True)

      if keep_headers:
        system._header_data = header_data

      if not resfilter or residue.name in resfilter:
        chain.add_residue(residue)
      elif skipping_residue_warning:
        mout.warning(f'Skipping residue {residue} {residue.number}')

      system.add_chain(chain)

      if remove_empty_chains:
        system.chains = [c for c in system.chains if c.residues]

      if prune_alternative_sites:
        system.prune_alternative_sites()

      if fix_atomnames:
        system.fix_atomnames(verbosity=1)

      if autoname_chains:
        system.autoname_chains()

      if num_appended_hydrogen_chains:

        if reordering_summary_warnings:
          mout.warning(f'Will try to re-order hydrogens in last {num_appended_hydrogen_chains} chains')

        assert num_appended_hydrogen_chains < len(system.chains), mout.error('num_appended_hydrogen_chains >= len(system.chains)')

        hydrogen_chains = []
        re_add_chains = []
        for i in range(num_appended_hydrogen_chains):
          
          chain = system.chains.pop(-1)

          if any([a.symbol != 'H' for a in chain.atoms]):
            re_add_chains.append(chain)
            continue

          if reordering_summary_warnings:
            mout.warning(f'Re-ordering hydrogens in {chain=}')
          hydrogen_chains.append(chain)

        for chain in re_add_chains:
          system.add_chain(chain)

        for h_chain in hydrogen_chains:
          for atom in h_chain.atoms:

            # chain = system.get_chain(atom.chain)
            chains = [c for c in system.chains if c.name == atom.chain]
            for chain in chains:
              residue = chain.residues[f'n{atom.res_number}']
              if residue is not None:
                break
            
            if reordering_warnings:
              mout.warning(f'Moving {atom} {atom.number} (i) to {residue} {residue.number} (res_indices {atom.res_index}-->{residue.index})')
            
            if residue is not None:
              residue.add_atom(atom)
            else:
              mout.warning(f"Didn't find a home for {atom} {atom.number} ({pdb=})")
              mout.warning(f"Skipping what appears to be a normal chain")
              system.add_chain(h_chain)
              break

      if detect_misordered_chains:
        for chain in system.chains:
          for i,res in enumerate(chain.residues):

            if res.number in [r.number for r in chain.residues[:i]]:
              mout.warning(f'{res} {res.number} has already appeared in chain {chain} {chain.index}')

      if fix_indices:
        system.fix_indices()

      if (verbosity > 0):
        mout.out("Done.") # user output

      return system

  except KeyboardInterrupt:
    mout.errorOut("KeyboardInterrupt")
    exit()

def parsePDBAtomLine(line,res_index,atom_index,chain_counter,debug=False,alternative_site_warnings=True):
  from .atom import Atom
  import mout
  import string
  from .atom import Atom

  if debug:
    mout.debug("Attempting to parse atom with index: "+str(atom_index))
    print(line.strip())

  charge = None

  try:
    atom_name = line[12:16].strip()
    if debug: mout.var(f'{atom_index}.name',atom_name)

    residue = line[17:21].strip()
    if debug: mout.var(f'{atom_index}.residue',residue)

    try:
      pdb_index = int(line[6:12].strip())
    except:
      pdb_index = atom_index
    if debug: mout.var(f'{atom_index}.pdb_index',pdb_index)
    chain = line[21:22]
    if chain == ' ':
      chain = string.ascii_uppercase[chain_counter%26]
    if debug: mout.var(f'{atom_index}.chain',chain)
    res_number = line[22:26].strip()
    if debug: mout.var(f'{atom_index}.res_number',res_number)
      
    alt_site_str = line[16:17]
    if len(line[16:17].strip()) == 0:
      alt_site_str = None

    if alternative_site_warnings and alt_site_str:
      mout.warningOut(f"Alternative site in PDB! res={residue}, atom={atom_name}, res_number={res_number}")
    if debug:
      mout.var(f'{atom_index}.alt_site_str',alt_site_str)

    position = []
    position.append(float(line[30:38].strip()))
    position.append(float(line[38:46].strip()))
    position.append(float(line[46:54].strip()))
    if debug: mout.var(f'{atom_index}.position',position)

    try:
      occupancy = float(line[55:61].strip())
    except:
      occupancy = None
    if debug: mout.var(f'{atom_index}.occupancy',occupancy)

    try:
      temp_factor = float(line[61:67].strip())
    except:
      temp_factor = None
    if debug: mout.var(f'{atom_index}.temp_factor',temp_factor)

    element = line[76:78].strip()
    if len(element) > 1:
      element = f'{element[0]}{element[1:].lower()}'
    if debug: mout.var(f'{atom_index}.element',f'###{element}###')

    chg_str = line[78:80].rstrip('\n')
    if len(chg_str.strip()) < 1:
      chg_str = None
    if debug: mout.var(f'{atom_index}.chg_str',chg_str)

    end = line[80:].strip()

    if line.startswith("HETATM"):
      hetatm=True
    else:
      hetatm=False
    if debug: mout.var(f'{atom_index}.hetatm',hetatm)

    if end:
      charge = float(line[80:89])

    atom = Atom(atom_name,pdb_index,pdb_index,position,residue,chain,res_number,occupancy=occupancy,temp_factor=temp_factor,heterogen=hetatm,charge_str=chg_str,alternative_site=alt_site_str,res_index=res_index,charge=charge,element=element)

    if debug:
      atom.summary()

    if debug:
      raise Exception('Stopping DEBUG')

    return atom

  except Exception as e:
    mout.error(e)
    mout.error("Unsupported PDB line shown below:")
    mout.error(line)
    raise Exception("Unsupported PDB line shown above")

# def parseGRO(gro,systemName=None,fix_indices=True,fix_atomnames=True,verbosity=1,auto_ter=None):
def parseGRO(gro,systemName=None,fix_indices=True,fix_atomnames=True,autoname_chains=True,verbosity=1,auto_ter=["DA3","DT3","DG3","DC3"]):
  import mcol
  import mout
  from .system import System
  from .chain import Chain

  assert gro.endswith(".gro")

  if (verbosity > 0):
    mout.out("parsing "+mcol.file+
             gro+
             mcol.clear+" ... ",
             end='') # user output

  import os

  residue = None
  chain = None
  last_residue_name = None
  last_residue_number = None
  last_chain_name = None
  res_counter = 0
  chain_counter = 0
  atom_counter = 1
  was_terminal = False

  if systemName is None:
    systemName = os.path.splitext(gro)[0]
  system = System(name = systemName)

  try: 
    with open(gro,"r") as input_gro:
      first = True
      line_counter = 0
      # make_new_res = True
      for line in input_gro:
        line_counter += 1

        # parse the header
        if line_counter == 1:
          sys_description = line.strip()
        elif line_counter == 2:
          sys_atom_count = int(line.strip())

        # parse the rest
        else:

          # check if last line
          if len(line.strip()) < 40:
            split_line = line.strip().split()
            system.box = [float(split_line[0]),float(split_line[1]),float(split_line[2])]
            break

          # parse an atom line:
          atom = parseGROAtomLine(line,res_counter,atom_counter,chain_counter)

          # first atom
          if line_counter == 3:
            chain = Chain(atom.chain)
            index = res_counter
            residue = new_residue(atom.residue,index,atom.res_number,atom.chain)
            # residue.addAtom(atom)
            last_residue_name = atom.residue
            last_residue_number = atom.res_number
            last_chain_name = atom.chain

          # check if a new residue is needed
          make_new_res = False
          if residue is None: make_new_res = True
          if last_residue_name != atom.residue: make_new_res = True
          if last_residue_number != atom.res_number: make_new_res = True

          # make the new residue
          if make_new_res:
            if residue is not None: 
              chain.add_residue(residue)
              res_counter += 1
            residue = new_residue(atom.residue,index,atom.res_number,atom.chain)

          # add the atom to the residue
          residue.addAtom(atom)
          residue_type = residue.type
          
          # check if a new chain is needed
          make_new_chain = False
          if residue is None: 
            # print("NONE-TER")
            make_new_chain = True
          if auto_ter is not None and make_new_res and last_residue_name in auto_ter:
            make_new_chain = True
            # print("AUTO-TER",residue.name)
          if make_new_res and residue_type != last_residue_type:
            make_new_chain = True
            # print("TYPE-TER",residue.name)
          
          # make the new chain
          if make_new_res:
            if make_new_chain:
              if chain is not None:
                system.add_chain(chain)
                chain_counter += 1
              chain = Chain(atom.chain)

          # update variables
          atom_counter += 1
          last_residue_type = residue_type
          last_residue_number = atom.res_number
          last_residue_name = atom.residue

          system.description = sys_description

  except FileNotFoundError:
    mout.errorOut("File "+mcol.file+gro+mcol.error+" not found",fatal=True)

  chain.add_residue(residue)
  system.add_chain(chain)

  if fix_indices:
    system.fix_indices()

  if fix_atomnames:
    system.fix_atomnames()

  if autoname_chains:
    system.autoname_chains()

  if (verbosity > 0):
    mout.out("Done.") # user output

  return system

def new_residue(name,index,number,chain):
  from .residue import res_type

  if res_type(name) == "PRO":
    if name in ['HIS','HID','HIE','HSE','HSD','HSP']:
      from .histidine import Histidine
      return Histidine(name,index,number,chain)
    else:
      from .amino import AminoAcid
      return AminoAcid(name,index,number,chain)
  elif res_type(name) == "DNA":
    from .nucleic import NucleicAcid
    return NucleicAcid(name,index,number,chain)
  elif name in ['ATP','GTP','CTP','TTP','OGTP']:
    from .ntp import NucleobaseTriPhosphate
    return NucleobaseTriPhosphate(name,index,number,chain)
  else:
    from .residue import Residue
    return Residue(name,index,number,chain)

def parseGROAtomLine(line,res_index,atom_index,chain_counter):

  import string
  from .atom import Atom

  res_number = line[0:5].strip()
  residue = line[5:10].strip()
  atom_name = line[11:15].strip()
  gro_index = line[15:21].strip()

  chain = string.ascii_uppercase[chain_counter%26]

  position = []
  # position needs conversion nm --> Angstrom
  position.append(10.0*float(line[21:29].strip()))
  position.append(10.0*float(line[29:37].strip()))
  position.append(10.0*float(line[37:45].strip()))
  
  velocity = []
  try:
    # velocity needs conversion: nm/ps --> Angstrom/fs
    velocity.append(0.01*float(line[45:53].strip()))
    velocity.append(0.01*float(line[53:61].strip()))
    velocity.append(0.01*float(line[61:69].strip()))
  except:
    velocity = [0.0,0.0,0.0]

  hetatm=False

  atom = Atom(atom_name,atom_index,gro_index,position,residue,chain,res_number,velocity=velocity,res_index=res_index)

  return atom

def writeCJSON(filename,system,use_atom_types=False,gulp_names=False,noPrime=False,verbosity=1):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             end='') # user output

  # Check that the input is the correct class
  assert isinstance(system,System)

  # Load module
  import json

  if gulp_names:
    use_atom_types = True
    names = []
    temp_names = system.FF_atomtypes
    for name in temp_names:
      name = name[0]+"_"+name[1:]
      names.append(name)
  else:
    if not use_atom_types:
      names = system.atom_names(wRes=True,noPrime=noPrime)
    else:
      names = system.FF_atomtypes

  # Create the dictionary
  data = {}
  data['chemical json'] = 0
  data['name'] = system.name
  data['atoms'] = {
    'names' : names,
    'elements' : {
      'number' : system.atomic_numbers
    },
    'coords' : {
      'unit' : "angstrom",
      '3d' : [item for sublist in system.positions for item in sublist]
    },
    "charges" : system.charges
  }

  # Write the CJSON file
  with open(filename, 'w') as f: json.dump(data, f,indent=4)

  if (verbosity > 0):
    mout.out("Done.") # user outpu

def writePDB(filename,system,verbosity=1,append=False,model=1,header=None,title=None,compound=None,write_model=True, charges=True, shift_name=0):
  import mcol
  import mout
  from .system import System
  from .group import AtomGroup

  if (verbosity > 0):
    mout.out(f"writing {mcol.file}{filename}{mcol.clear} ... ",
             end='') # user output

  # Check that the input is the correct class
  assert isinstance(system,System) or isinstance(system,AtomGroup)

  end = '\n'

  if not append:

    if system._header_data is not None:

      strbuff = "".join(system._header_data)      

    else:

      strbuff = ""

      if not header:
        strbuff += f"HEADER {filename}{end}"
      else:
        strbuff += f"HEADER {header}{end}"

      if not title:
        strbuff += f"TITLE  {system.name}{end}"
      else:
        strbuff += f"TITLE  {title}{end}"

      strbuff += "REMARK "+"generated by molparse.io.writePDB()"+end
      # strbuff += "REMARK "+end
      # strbuff += "REMARK "+"System Summary:"+end
      if isinstance(system, System):
        strbuff += "REMARK "+"# Chains:   "+str(system.num_chains)+end
        strbuff += "REMARK "+"# Residues: "+str(system.num_residues)+end
        if system.remarks:
          for line in system.remarks:
            strbuff += f"REMARK {line}"+end

      strbuff += "REMARK "+"# Atoms:    "+str(system.num_atoms)+end

  else:
    strbuff = ""

  if write_model:
    strbuff += "MODEL "+str(model)+end

  for i,atom in enumerate(system.atoms):

    atomline = constructPDBAtomLine(atom,i, charges=charges, shift_name=shift_name)

    strbuff += atomline

    if atom.terminal:
      # print(atom, atom.number, 'TERMINAL2')
      strbuff += atom.ter_line

  if write_model:
    strbuff += "ENDMDL"+end

  if append:
    out_stream = open(filename,"a")
  else:
    out_stream = open(filename,"w")

  out_stream.write(strbuff)
  out_stream.close()

  if (verbosity > 0):
    mout.out("Done.") # user output    

def modifyPDB(filename,atoms,copy_from=None):

  """ in-place modification of a PDB file for improved performance when only modifying a certain subset of atoms 

  filename: PDB to modify
  atoms: list of mp.Atom objects that will be modified in the PDB
  copy_from: (optional) copy {copy_from} --> {filename} before modifying

  """

  if copy_from:
    import shutil
    shutil.copyfile(copy_from,filename)

  import fileinput

  pdb_indices = [a.pdb_index for a in atoms]

  max_index = max(pdb_indices)
  min_index = min(pdb_indices)

  for line in fileinput.FileInput(filename, inplace=1):

    try:
      pdb_index = int(line[6:12].strip())

    except ValueError:
      print(line,end='')
      continue
      
    # check if the atom index matches a pdb_index in the atoms
    if pdb_index >= min_index and pdb_index <= max_index and pdb_index in pdb_indices:
      line = constructPDBAtomLine(atoms[pdb_indices.index(pdb_index)],pdb_index)

    print(line,end='')

def constructPDBAtomLine(atom,index,charges=True, shift_name=False, alt_sites=True):

  import mout
  import mcol

  end = '\n'

  strlist = []

  atom_serial = atom.number or atom.index
  residue_serial = atom.res_number

  if atom_serial is None:
    mout.warningOut(f"{mcol.varName}atom_serial{mcol.clear}{mcol.warning} is None",code=f'mp.io.constructPDBAtomLine({atom},{index})')
    atom_serial = index
  if residue_serial is None:
    mout.errorOut(f"{mcol.varName}residue_serial{mcol.clear}{mcol.warning} is None",code=f'mp.io.constructPDBAtomLine({atom},{index})')
    residue_serial = "    "

  if not atom.heterogen:
    strlist.append("ATOM  ")
  else:
    strlist.append("HETATM")

  atom_serial_str = str(atom_serial).rjust(5)
  if len(atom_serial_str) > 5: 
    atom_serial_str = "XXXXX"

  residue_serial_str = str(residue_serial).rjust(4)
  if len(residue_serial_str) > 4: 
    residue_serial_str = residue_serial_str[-4:]

  strlist.append(atom_serial_str)
  strlist.append(" ")
  if shift_name:
    strlist.append(f' {atom.name[:3].ljust(3)}')
  else:
    strlist.append(str(atom.name[:4]).ljust(4))
  if alt_sites and atom.alternative_site:
    strlist.append(str(atom.alternative_site))
  else:
    strlist.append(" ")
  # strlist.append(" ")
  if atom.residue is None:
    mout.errorOut(f"{mcol.varName}atom.residue{mcol.clear}{mcol.warning} is None",code=f'mp.io.constructPDBAtomLine({atom},{index})')
    strlist.append("    ".ljust(4))
  else:
    res_str = str(atom.residue).ljust(4)
    if len(res_str) > 4:
      res_str = res_str[:4]
    strlist.append(res_str)
  if atom.chain is None:
    mout.errorOut(f"{mcol.varName}atom.chain{mcol.clear}{mcol.warning} is None",code=f'mp.io.constructPDBAtomLine({atom},{index})')
    strlist.append(" ")
  else:
    strlist.append(str(atom.chain))
  strlist.append(residue_serial_str)
  strlist.append("    ")

  x_str = '{:.3f}'.format(atom.x).rjust(8)
  y_str = '{:.3f}'.format(atom.y).rjust(8)
  z_str = '{:.3f}'.format(atom.z).rjust(8)

  strlist.append(x_str+y_str+z_str)

  if atom.occupancy is not None:
    strlist.append('{:.2f}'.format(atom.occupancy).rjust(6))
  else:
    strlist.append('      ')
  if atom.temp_factor is not None:
    strlist.append('{:.2f}'.format(atom.temp_factor).rjust(6))
  else:
    strlist.append('      ')

  strlist.append("          ")
  strlist.append(atom.species.rjust(2))
  
  if atom.charge_str is not None:
    strlist.append(atom.charge_str)
  # else:
    # strlist.append("  ")

  if charges and atom.charge is not None:
    strlist.append(f'{atom.charge:10.6f}')

  strlist.append(end)

  if atom.terminal:
    # print(atom, atom.number, 'TERMINAL1')
    atom_serial_str = str(atom_serial+1).rjust(5)
    if len(atom_serial_str) > 5: 
      atom_serial_str = "XXXXX"
    
    atom.ter_line =  "TER   "
    atom.ter_line += atom_serial_str
    atom.ter_line += "      "
    atom.ter_line += atom.residue.ljust(4)
    atom.ter_line += atom.chain
    atom.ter_line += str(residue_serial).rjust(4)
    atom.ter_line += end

  return ''.join(strlist)

def writeGRO(filename,system,verbosity=1):
  import mcol
  import mout
  from .system import System

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             end='') # user output

  # Check that the input is the correct class
  assert isinstance(system,System)

  end = '\n'

  strbuff =  system.name+" (amp.io)"+end
  strbuff += str(system.num_atoms)+end

  atom_serial = 1
  residue_serial = 1

  for chain in system.chains:
    for residue in chain.residues:
      for atom in residue.atoms:

        residue_serial_str = str(atom.res_number).rjust(5)
        if len(residue_serial_str) > 5: 
          residue_serial_str = residue_serial_str[-5:]

        index = atom.pdb_index or atom.index
        atom_serial_str = str(index).rjust(4)
        if len(atom_serial_str) > 4: 
          atom_serial_str = atom_serial_str[-4:]

        strbuff += residue_serial_str
        strbuff += str(atom.residue).ljust(4)
        strbuff += " "
        strbuff += str(atom.name).rjust(5)
        strbuff += " "
        strbuff += atom_serial_str
        
        x_str = '{:.3f}'.format(atom.x/10.0).rjust(8)
        y_str = '{:.3f}'.format(atom.y/10.0).rjust(8)
        z_str = '{:.3f}'.format(atom.z/10.0).rjust(8)
        strbuff += x_str+y_str+z_str
        
        if atom.velocity is None:
          x_str = '{:.4f}'.format(0.0).rjust(8)
          y_str = '{:.4f}'.format(0.0).rjust(8)
          z_str = '{:.4f}'.format(0.0).rjust(8)
        else:
          x_str = '{:.4f}'.format(atom.velocity[0]*100.0).rjust(8)
          y_str = '{:.4f}'.format(atom.velocity[1]*100.0).rjust(8)
          z_str = '{:.4f}'.format(atom.velocity[2]*100.0).rjust(8)
        strbuff += x_str+y_str+z_str

        strbuff += end

        atom_serial += 1

      residue_serial += 1

  if system.box is None:
    mout.warningOut("System has no box information. Using bbox")
    bbox = system.bbox
    x_str = '{:.5f}'.format(bbox[0][1]-bbox[0][0]).rjust(10)
    y_str = '{:.5f}'.format(bbox[1][1]-bbox[1][0]).rjust(10)
    z_str = '{:.5f}'.format(bbox[2][1]-bbox[2][0]).rjust(10)
  else:
    x_str = '{:.5f}'.format(system.box[0]).rjust(10)
    y_str = '{:.5f}'.format(system.box[1]).rjust(10)
    z_str = '{:.5f}'.format(system.box[2]).rjust(10)
  strbuff += x_str+y_str+z_str
  strbuff += end

  out_stream = open(filename,"w")

  out_stream.write(strbuff)
  out_stream.close()

  if (verbosity > 0):
    mout.out("Done.") # user output
    
def parseXYZ(xyz,index=":",verbosity=1):

  import mout
  import mcol

  assert xyz.endswith(".xyz")

  try:
    index = int(index)
  except:
    if index == ":":
      if verbosity > 0:
        mout.warningOut("Parsing all models in "+mcol.file+xyz)
    else:
      mout.errorOut("Unsupported index: '"+str(index)+"'",fatal=True)

  if index == ":":

    all_models = []

    import subprocess
    num_lines = int(subprocess.check_output(f"cat {xyz} | wc -l", shell=True))

    with open(xyz,"r") as input_xyz:
      for line in input_xyz:
        num_atoms = int(line)
        break

    if (num_lines/(num_atoms+2)) != (num_lines//(num_atoms+2)):
      mout.errorOut("Wrong number of lines, check XYZ formatting is correct.",fatal=True)

    num_models = num_lines//(num_atoms+2)

    for i in range(num_models):

      all_models.append(parseXYZ(xyz,index=i,verbosity=verbosity-1))

    return all_models

  elif index < 0:

    import subprocess
    num_lines = int(subprocess.check_output(f"cat {xyz} | wc -l", shell=True))

    with open(xyz,"r") as input_xyz:
      for line in input_xyz:
        num_atoms = int(line)
        break

    if (num_lines/(num_atoms+2)) != (num_lines//(num_atoms+2)):
      mout.errorOut("Wrong number of lines, check XYZ formatting is correct.",fatal=True)

    num_models = num_lines//(num_atoms+2)

    return parseXYZ(xyz,index=num_models+index-1,verbosity=verbosity)

  else:

    from .group import AtomGroup

    with open(xyz,"r") as input_xyz:

      start_line = 0

      for j,line in enumerate(input_xyz):

        # print(line)
        
        if j == 0:
          num_atoms = int(line)
          start_line = index*(num_atoms+2)
          continue

        elif j == start_line + 1:

          try:
            # get the model information
            i_str, E_str = line.strip().split(',')
            i = int(i_str.split("=")[-1])
            E = float(E_str.split("=")[-1])
            system = AtomGroup(f"{xyz}, image={i}")
            system._energy = E
            system._traj_index = i
          except ValueError:
            mout.warningOut("XYZ image header string has an unsupported format")
            system = AtomGroup(f"{xyz}, {line.strip()}")

          # create the system object
          atom_counter = 0
          continue

        elif j < start_line + 1:
          continue

        if num_atoms == atom_counter:
          break

        atom = parseXYZAtomLine(line, atom_counter)
        system.add_atom(atom)
        atom_counter += 1

      if j < start_line:
        mout.errorOut(f"XYZ has no image with index {index}")
    
    return system

def parseXYZAtomLine(line,atom_index):

  from .atom import Atom

  split_line = line.strip().split()

  if len(split_line) == 4:
    s,x,y,z = split_line
  else:
    s,x,y,z = split_line[:4]

  x = float(x)
  y = float(y)
  z = float(z)

  atom = Atom(s,index=atom_index,position=[x,y,z])

  return atom

