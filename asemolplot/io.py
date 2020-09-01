from ase import io
from ase import atoms as aseatoms

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

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

  if filename.endswith(".cjson"):
    writeCJSON(filename,image)
  elif filename.endswith(".pdb") and isinstance(image,System):
    writePDB(filename,image)
  else:
    io.write(filename,image,**parameters)

  if (verbosity > 0):
    mout.out("Done.") # user output

def read(filename,index=None,printScript=False,verbosity=1,tagging=True,tagByResidue=False,**parameters):

  if (verbosity > 0):
    mout.out("reading "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  atoms = io.read(filename,index,**parameters)

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

def parsePDB(pdb,systemName=None,fix_indices=True):

  import os

  # file  = open(pdb, 'r').read()
  # if file.count("MODEL") == 0:
  #   mout.errorOut("PDB has no MODELs",fatal=True)
  # file.close()

  residue = None
  chain = None
  last_residue_name = None
  last_chain_name = None
  res_counter = 0

  if systemName is None:
    systemName = os.path.splitext(pdb)[0]
  system = System(name = systemName)

  with open(pdb,"r") as input_pdb:
    searching = True
    for line in input_pdb:
      if searching:
        if line.startswith("MODEL"):
          searching = False
        elif line.startswith("ATOM"):
          searching = False
          #### PARSELINE
          atom = parsePDBAtomLine(line,res_counter)
          chain = Chain(atom.chain)
          residue = Residue(atom.residue,res_counter,atom.chain)
          residue.addAtom(atom)
          last_residue_name = atom.residue
          last_chain_name = atom.chain
      else:
        if line.startswith("ENDMDL"):
          break
        if line.startswith("END"):
          break
        if line.startswith("TER"):
          continue
        if line.startswith("CONECT"):
          continue
        if line.startswith("MASTER"):
          continue
        else:
          ### PARSELINE
          atom = parsePDBAtomLine(line,res_counter)
          
          make_new_res = False
          if residue is None: make_new_res = True
          if last_residue_name != atom.residue: make_new_res = True

          make_new_chain = False
          if residue is None: make_new_chain = True
          if last_chain_name != atom.chain: make_new_chain = True

          if make_new_res:
            if residue is not None: 
              chain.addResidue(residue)
              res_counter = res_counter+1
            residue = Residue(atom.residue,res_counter,atom.chain)
            if make_new_chain:
              if chain is not None:
                system.addChain(chain)
              chain = Chain(atom.chain)

          residue.addAtom(atom)

          last_residue_name = atom.residue
          last_chain_name = atom.chain

  chain.addResidue(residue)
  system.addChain(chain)

  if fix_indices:
    system.fix_indices()

  return system

def parsePDBAtomLine(line,index):

  try:
    atom_name = line[12:17].strip()
    residue = line[17:21].strip()
    pdb_index = int(line[6:12].strip())
    chain = line[21:22]
    res_number = line[22:26].strip()
    index = line[7:12].strip()

    position = []

    # print("#"+line[31:39]+"#")


    position.append(float(line[31:39].strip()))
    position.append(float(line[39:47].strip()))
    position.append(float(line[47:55].strip()))

    occupancy = float(line[55:61].strip())
    temp_factor = float(line[61:67].strip())

    chg_str = line[78:80].rstrip('\n')

    end = line[80:]

    if line.startswith("HETATM"):
      hetatm=True
    else:
      hetatm=False

    if 'QM' in end:
      isQM = True
    else:
      isQM = False

    atom = Atom(atom_name,index,pdb_index,position,residue,chain,res_number,QM=isQM,occupancy=occupancy,temp_factor=temp_factor,heterogen=hetatm,charge_str=chg_str)
    
    return atom
  except:
    mout.errorOut("Unsupported PDB line "+str(index)+":")
    print(line)

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

def writePDB(filename,system,verbosity=1,printScript=False):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  # Check that the input is the correct class
  assert isinstance(system,System)

  end = '\n'

  strbuff =  "HEADER "+filename+end
  strbuff += "TITLE  "+system.name+end
  strbuff += "REMARK "+"generated by asemolplot.io.writePDB()"+end
  # strbuff += "REMARK "+end
  # strbuff += "REMARK "+"System Summary:"+end
  strbuff += "REMARK "+"# Chains:   "+str(system.num_chains)+end
  strbuff += "REMARK "+"# Residues: "+str(system.num_residues)+end
  strbuff += "REMARK "+"# Atoms:    "+str(system.num_atoms)+end

  strbuff += "MODEL 1"+end

  atom_serial = 1
  residue_serial = 1

  for chain in system.chains:
    for residue in chain.residues:
      for atom in residue.atoms:

        if not atom.heterogen:
          strbuff += "ATOM  "
        else:
          strbuff += "HETATM"

        strbuff += str(atom_serial).rjust(5)
        strbuff += " "
        strbuff += str(atom.name).rjust(4)
        strbuff += " "
        strbuff += str(atom.residue).ljust(4)
        assert len(chain.name) == 1
        strbuff += str(chain.name)
        strbuff += str(residue_serial).rjust(4)
        strbuff += "    "

        x_str = '{:.3f}'.format(atom.x).rjust(8)
        y_str = '{:.3f}'.format(atom.y).rjust(8)
        z_str = '{:.3f}'.format(atom.z).rjust(8)

        # x_str = mout.toPrecision(atom.x,3,sf=False).rjust(8)
        # y_str = mout.toPrecision(atom.y,3,sf=False).rjust(8)
        # z_str = mout.toPrecision(atom.z,3,sf=False).rjust(8)

        strbuff += x_str+y_str+z_str

        strbuff += '{:.2f}'.format(atom.occupancy).rjust(6)
        strbuff += '{:.2f}'.format(atom.temp_factor).rjust(6)
        strbuff += "          "
        strbuff += atom.species.rjust(2)
        strbuff += atom.charge_str

        if atom.QM:
          strbuff += "QM"

        strbuff += end

        atom_serial += 1

      residue_serial += 1

  strbuff += "ENDMDL"+end

  out_stream = open(filename,"w")
  out_stream.write(strbuff)
  out_stream.close()

  if (verbosity > 0):
    mout.out("Done.") # user outpu