from ase import io
from ase import atoms as aseatoms

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

import string

# Custom Classes
from .system import System
from .chain import Chain
from .residue import Residue
from .atom import Atom

def write(filename,image,verbosity=1,printScript=False,**parameters):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  # Different behaviour depending on format:
  try:
    # CJSON:
    if filename.endswith(".cjson"):
      writeCJSON(filename,image)
    # PBD and amp.System:
    elif filename.endswith(".pdb") or filename.endswith(".pdb2") and isinstance(image,System):
      writePDB(filename,image,verbosity=verbosity-1)
    # GRO and amp.System:
    elif filename.endswith(".gro") and isinstance(image,System):
      writeGRO(filename,image,verbosity=verbosity-1)
    # PDB and list of amp.System's:
    elif filename.endswith(".pdb") and isinstance(image,list):

      ### List of Systems
      if all([isinstance(i,System) for i in image ]):
        mout.out("Cool")

        for index,frame in enumerate(image):
          if index == 0:
            writePDB(filename,frame,verbosity=verbosity-1)
          else:
            writePDB(filename,frame,verbosity=verbosity-1,append=True,model=index+1)

      else:
        mout.errorOut("Unsupported",fatal=True)

    # Others:
    else:
      io.write(filename,image,**parameters)
  except TypeError:
    mout.out("Fail.")
    mout.errorOut("Could not write empty object to "+mcol.file+filename,code="amp.io.write[1]",end="")
    return None

  if (verbosity > 0):
    mout.out("Done.") # user output

def read(filename,index=None,printScript=False,verbosity=1,tagging=True,tagByResidue=False,**parameters):

  if (verbosity > 0):
    mout.out("reading "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
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

  print(len(taglist))
  print(len(atoms))

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

def parsePDB(pdb,systemName=None,fix_indices=True,fix_atomnames=True,verbosity=1,debug=False):

  if (verbosity > 0):
    mout.out("parsing "+mcol.file+
             pdb+
             mcol.clear+" ... ",
             end='') # user output

  if debug:
    mout.warningOut("Debug mode active!")

  import os

  # file  = open(pdb, 'r').read()
  # if file.count("MODEL") == 0:
  #   mout.errorOut("PDB has no MODELs",fatal=True)
  # file.close()

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
    systemName = os.path.splitext(pdb)[0]
  system = System(name = systemName)

  try: 
    with open(pdb,"r") as input_pdb:
      searching = True
      for line in input_pdb:
        if searching:
          if line.startswith("MODEL"):
            searching = False
          elif line.startswith("ATOM"):
            searching = False
            #### PARSELINE
            atom = parsePDBAtomLine(line,res_counter,atom_counter,chain_counter,debug=debug)
            chain = Chain(atom.chain)
            residue = Residue(atom.residue,res_counter,atom.chain)
            residue.addAtom(atom)
            last_residue_name = atom.residue
            last_residue_number = atom.res_number
            last_chain_name = atom.chain
        else:
          if line.startswith("ENDMDL"):
            break
          if line.startswith("END"):
            break
          if line.startswith("TER"):
            if verbosity > 1:
              mout.warningOut("Terminal added to "+mcol.arg+chain.name+":"+residue.name)
            residue.atoms[-1].terminal = True
            was_terminal = True
            # atom = residue.atoms[-1]
            # # residue.atoms[-1].ter_line=line

            # atom.terminal = True
            # atom.ter_line =  "TER   "
            # atom.ter_line += str(atom.index)
            # atom.ter_line += "\n"

            # residue.atoms[-1].ter_line="TER   "+x_str = '{:.3f}'.format(atom.x).rjust(8)+"\n"
            # # TER    4060      ILE A 509  
            continue
          if line.startswith("CONECT"):
            continue
          if line.startswith("MASTER"):
            continue
          else:
            ### PARSELINE
            atom = parsePDBAtomLine(line,res_counter,atom_counter,chain_counter,debug=debug)
            
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
                chain.add_residue(residue)
                res_counter = res_counter+1
              residue = Residue(atom.residue,res_counter,atom.chain)
              if make_new_chain:
                if chain is not None:
                  system.add_chain(chain)
                  chain_counter = chain_counter+1
                chain = Chain(atom.chain)

            residue.addAtom(atom)

            atom_counter += 1

            last_residue_number = atom.res_number
            last_residue_name = atom.residue
            last_chain_name = atom.chain
            was_terminal = False

  except FileNotFoundError:
    mout.errorOut("File "+mcol.file+pdb+mcol.error+" not found",fatal=True)

  chain.add_residue(residue)
  system.add_chain(chain)

  if fix_indices:
    system.fix_indices()

  if fix_atomnames:
    system.fix_atomnames()

  if (verbosity > 0):
    mout.out("Done.") # user output

  return system

def parsePDBAtomLine(line,res_index,atom_index,chain_counter,debug=False):

  if debug:
    mout.out("Attempting to parse atom with index: "+str(atom_index))

  try:
    atom_name = line[12:17].strip()
    if debug: print(str(atom_index) + ".name: OK")
    residue = line[17:21].strip()
    if debug: print(str(atom_index) + ".residue: OK")
    try:
      pdb_index = int(line[6:12].strip())
    except:
      pdb_index = atom_index
    if debug: print(str(atom_index) + ".index: OK")
    chain = line[21:22]
    if chain == ' ':
      chain = string.ascii_uppercase[chain_counter%26]
    if debug: print(str(atom_index) + ".chain: OK")
    res_number = line[22:26].strip()
    if debug: print(str(atom_index) + ".res_number: OK")

    position = []
    position.append(float(line[31:39].strip()))
    position.append(float(line[39:47].strip()))
    position.append(float(line[47:55].strip()))
    if debug: print(str(atom_index) + ".position: OK")

    try:
      occupancy = float(line[55:61].strip())
    except:
      occupancy = None
    if debug: print(str(atom_index) + ".occupancy: OK")

    try:
      temp_factor = float(line[61:67].strip())
    except:
      temp_factor = None
    if debug: print(str(atom_index) + ".temp_factor: OK")

    chg_str = line[78:80].rstrip('\n')
    if debug: print(str(atom_index) + ".chg_str: OK")

    end = line[80:]

    if line.startswith("HETATM"):
      hetatm=True
    else:
      hetatm=False
    if debug: print(str(atom_index) + ".hetatm: OK")

    if 'QM' in end:
      isQM = True
    else:
      isQM = False
    if debug: print(str(atom_index) + ".isqm: OK")

    atom = Atom(atom_name,pdb_index,pdb_index,position,residue,chain,res_number,QM=isQM,occupancy=occupancy,temp_factor=temp_factor,heterogen=hetatm,charge_str=chg_str)

    return atom

  except:
    mout.errorOut(line)
    mout.errorOut("Unsupported PDB line shown above",fatal=True)

# def parseGRO(gro,systemName=None,fix_indices=True,fix_atomnames=True,verbosity=1,auto_ter=None):
def parseGRO(gro,systemName=None,fix_indices=True,fix_atomnames=True,verbosity=1,auto_ter=["DA3","DT3","DG3","DC3"]):

  if (verbosity > 0):
    mout.out("parsing "+mcol.file+
             gro+
             mcol.clear+" ... ",
             end='') # user output

  import os

  # file  = open(gro, 'r').read()
  # if file.count("MODEL") == 0:
  #   mout.errorOut("gro has no MODELs",fatal=True)
  # file.close()

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
            residue = Residue(atom.residue,res_counter,atom.chain)
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
            residue = Residue(atom.residue,res_counter,atom.chain)

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

  if (verbosity > 0):
    mout.out("Done.") # user output

  return system

def parseGROAtomLine(line,res_index,atom_index,chain_counter):

  # try:
    res_number = line[0:5].strip()
    residue = line[5:10].strip()
    atom_name = line[11:15].strip()
    gro_index = line[15:21].strip()
    # print(res_number,residue,atom_name,atom_index)

    # assert gro_index == atom_index

    chain = string.ascii_uppercase[chain_counter%26]

    position = []
    position.append(float(line[21:29].strip()))
    position.append(float(line[29:37].strip()))
    position.append(float(line[37:45].strip()))
    
    velocity = []
    velocity.append(float(line[45:53].strip()))
    velocity.append(float(line[53:61].strip()))
    velocity.append(float(line[61:69].strip()))

    hetatm=False

    atom = Atom(atom_name,atom_index,atom_index,position,residue,chain,res_number,velocity=velocity)

    return atom

  # except:
  #   mout.errorOut(line)
  #   mout.errorOut("Unsupported GRO line shown above",fatal=True)

def writeCJSON(filename,system,use_atom_types=False,gulp_names=False,noPrime=False,printScript=False,verbosity=1):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
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

def writePDB(filename,system,verbosity=1,printScript=False,append=False,model=1):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  # Check that the input is the correct class
  assert isinstance(system,System)

  end = '\n'

  if not append:
    strbuff =  "HEADER "+filename+end
    strbuff += "TITLE  "+system.name+end
    strbuff += "REMARK "+"generated by asemolplot.io.writePDB()"+end
    # strbuff += "REMARK "+end
    # strbuff += "REMARK "+"System Summary:"+end
    strbuff += "REMARK "+"# Chains:   "+str(system.num_chains)+end
    strbuff += "REMARK "+"# Residues: "+str(system.num_residues)+end
    strbuff += "REMARK "+"# Atoms:    "+str(system.num_atoms)+end
  else:
    strbuff = ""

  strbuff += "MODEL "+str(model)+end

  atom_serial = 1
  residue_serial = 1

  for chain in system.chains:
    for residue in chain.residues:
      for atom in residue.atoms:

        if not atom.heterogen:
          strbuff += "ATOM  "
        else:
          strbuff += "HETATM"

        atom_serial_str = str(atom_serial).rjust(5)
        if len(atom_serial_str) > 5: 
          atom_serial_str = "XXXXX"

        residue_serial_str = str(residue_serial).rjust(4)
        if len(residue_serial_str) > 4: 
          # residue_serial_str = "XXXX"
          # residue_serial_str = "    "
          residue_serial_str = residue_serial_str[-4:]

        strbuff += atom_serial_str
        strbuff += " "
        strbuff += str(atom.name).rjust(4)
        strbuff += " "
        strbuff += str(atom.residue).ljust(4)
        assert len(chain.name) == 1
        strbuff += str(chain.name)
        strbuff += residue_serial_str
        strbuff += "    "

        x_str = '{:.3f}'.format(atom.x).rjust(8)
        y_str = '{:.3f}'.format(atom.y).rjust(8)
        z_str = '{:.3f}'.format(atom.z).rjust(8)

        # x_str = mout.toPrecision(atom.x,3,sf=False).rjust(8)
        # y_str = mout.toPrecision(atom.y,3,sf=False).rjust(8)
        # z_str = mout.toPrecision(atom.z,3,sf=False).rjust(8)

        strbuff += x_str+y_str+z_str

        if atom.occupancy is not None:
          strbuff += '{:.2f}'.format(atom.occupancy).rjust(6)
        else:
          strbuff += '      '
        if atom.temp_factor is not None:
          strbuff += '{:.2f}'.format(atom.temp_factor).rjust(6)
        else:
          strbuff += '      '

        strbuff += "          "
        strbuff += atom.species.rjust(2)
        
        if atom.charge_str is not None:
          strbuff += atom.charge_str

        if atom.QM:
          strbuff += "QM"

        strbuff += end


        if atom.terminal:
          atom.ter_line =  "TER   "
          atom.ter_line += atom_serial_str
          atom.ter_line += "      "
          atom.ter_line += atom.residue.ljust(4)
          atom.ter_line += atom.chain
          atom.ter_line += str(residue_serial).rjust(4)
          atom.ter_line += end
          strbuff += atom.ter_line

        atom_serial += 1

      residue_serial += 1

  strbuff += "ENDMDL"+end

  if append:
    out_stream = open(filename,"a")
  else:
    out_stream = open(filename,"w")

  out_stream.write(strbuff)
  out_stream.close()

  if (verbosity > 0):
    mout.out("Done.") # user output
    

def writeGRO(filename,system,verbosity=1,printScript=False):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
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

        # strbuff += "ATOM  "

        residue_serial_str = str(residue_serial).rjust(5)
        if len(residue_serial_str) > 5: 
          # residue_serial_str = "XXXX"
          # residue_serial_str = "    "
          residue_serial_str = residue_serial_str[-5:]

        atom_serial_str = str(atom_serial).rjust(4)
        if len(atom_serial_str) > 4: 
          atom_serial_str = atom_serial_str[-4:]
          # atom_serial_str = "XXXXX"

        strbuff += residue_serial_str
        strbuff += str(atom.residue).ljust(4)
        strbuff += " "
        strbuff += str(atom.name).rjust(5)
        strbuff += " "
        strbuff += atom_serial_str
        
        x_str = '{:.3f}'.format(atom.x).rjust(8)
        y_str = '{:.3f}'.format(atom.y).rjust(8)
        z_str = '{:.3f}'.format(atom.z).rjust(8)
        strbuff += x_str+y_str+z_str
        
        x_str = '{:.3f}'.format(atom.velocity[0]).rjust(8)
        y_str = '{:.3f}'.format(atom.velocity[1]).rjust(8)
        z_str = '{:.3f}'.format(atom.velocity[2]).rjust(8)
        strbuff += x_str+y_str+z_str

        strbuff += end

        atom_serial += 1

      residue_serial += 1

  x_str = '{:.3f}'.format(system.box[0]).rjust(12)
  y_str = '{:.3f}'.format(system.box[1]).rjust(12)
  z_str = '{:.3f}'.format(system.box[2]).rjust(12)
  strbuff += x_str+y_str+z_str
  strbuff += end

  out_stream = open(filename,"w")

  out_stream.write(strbuff)
  out_stream.close()

  if (verbosity > 0):
    mout.out("Done.") # user output
    