from ase import io
from ase import atoms as aseatoms

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

# Custom Classes
from .chain_res_atoms import System
from .chain_res_atoms import Chain
from .chain_res_atoms import Residue
from .chain_res_atoms import Atom

def write(filename,image,verbosity=1,printScript=False,**parameters):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  if filename.endswith(".cjson"):
    writeCJSON(filename,image)
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

def parsePDB(pdb,systemName=None):

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

  return system

def parsePDBAtomLine(line,index):

  atom_name = line[12:17].strip()
  residue = line[17:21].strip()
  pdb_index = int(line[5:12].strip())
  chain = line[21:22]
  res_number = line[22:26].strip()

  position = []
  position.append(float(line[31:39].strip()))
  position.append(float(line[39:47].strip()))
  position.append(float(line[47:55].strip()))

  atom = Atom(atom_name,index,pdb_index,position,residue,chain,res_number)
  
  return atom

def writeCJSON(filename,system,use_atom_types=False):

  # Check that the input is the correct class
  assert isinstance(system,System)

  # Load module
  import json

  if use_atom_types:
    names = system.atom_names(wRes=True,noPrime=True)
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
