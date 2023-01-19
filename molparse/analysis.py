
"""Undocumented. Speak to Max"""

def getAtomLabel(trajectory,index,full=False):
  atom_symbols = trajectory[0].get_chemical_symbols()
  atom_tags = trajectory[0].get_tags()

  atom_symbol = atom_symbols[index]
  if atom_tags[index] != 0:
    atom_tag = str(atom_tags[index])
  else:
    atom_tag = ""

  atom_label = atom_symbol+str(atom_tag)+"["+str(index)+"]"

  if full:
    return atom_label,atom_symbol,atom_tag
  else:
    return atom_label

def getBondLabel(trajectory,index_pair,full=False):

  index1=index_pair[0]
  index2=index_pair[1]

  atom_label1,atom_symbol1,atom_tag1 = getAtomLabel(trajectory,index1,full=True)
  atom_label2,atom_symbol2,atom_tag2 = getAtomLabel(trajectory,index2,full=True)

  bond_title = atom_symbol1+str(atom_tag1)+"-"+atom_symbol2+str(atom_tag2)+" bond ["+str(index1)+"-"+str(index2)+"]"

  if full:
    return bond_title,atom_label1,atom_label2
  else:
    return bond_title

def bondLengthStats(trajectory,index_pair,printScript=False,verbosity=1,timestep=None,yUnit="Angstroms",fitMin=None,fitMax=None,returnData=False,dataFile=None):
  
  import mcol
  import mout
  import mplot

  index1 = index_pair[0]
  index2 = index_pair[1]

  bond_title = getBondLabel(trajectory,index_pair)

  xdata=[]
  ydata=[]

  if timestep is None:
    xlab = "MD Steps"
    xUnit = "MD Step"
  else:  
    mout.warningOut(f"Scaling xdata by {timestep}x")
    xlab = "Time [fs]"
    xUnit = "femtoseconds"

  if len(trajectory) > 1:

    for n, atoms in enumerate(trajectory):

      if timestep is None:
        xdata.append(n)
      else:
        xdata.append(n*timestep)

      dist = atoms.get_distance(index1,index2)

      ydata.append(dist)

    val,err,fit_func = mplot.fit(xdata,ydata,rank=0,verbosity=verbosity-1,fitMin=fitMin,fitMax=fitMax,title=bond_title,yUnit=yUnit,xUnit=xUnit)

  else:

    val = trajectory[0].get_distance(index1,index2)
    err = None

  mout.varOut(bond_title,val,error=err,unit=yUnit,valCol=mcol.result,dataFile=dataFile,verbosity=verbosity)

  if returnData:
    return val, err, bond_title, xdata, ydata
  else:
    return val, err, bond_title

def bondAngleStats(trajectory,index_triplet,printScript=False,verbosity=1,timestep=None,fitMin=None,fitMax=None,yUnit="degrees",returnData=False,dataFile=None):

  import mcol
  import mout
  import mplot
  from ase.units import fs

  atom_symbols = trajectory[0].get_chemical_symbols()
  atom_tags = trajectory[0].get_tags()

  index1 = index_triplet[0]
  index2 = index_triplet[1]
  index3 = index_triplet[2]

  atom_symbol1 = atom_symbols[index1]
  atom_symbol2 = atom_symbols[index2]
  atom_symbol3 = atom_symbols[index3]

  if atom_tags[index1] != 0:
    atom_tag1 = str(atom_tags[index1])
  else:
    atom_tag1 = ""
  if atom_tags[index2] != 0:
    atom_tag2 = str(atom_tags[index2])
  else:
    atom_tag2 = ""
  if atom_tags[index3] != 0:
    atom_tag3 = str(atom_tags[index3])
  else:
    atom_tag3 = ""

  atom_label1 = atom_symbol1+str(atom_tag1)+"["+str(index1)+"]"
  atom_label2 = atom_symbol2+str(atom_tag2)+"["+str(index2)+"]"
  atom_label3 = atom_symbol3+str(atom_tag3)+"["+str(index3)+"]"
  labels = [atom_label1, atom_label2, atom_label3]
  bond_title = atom_symbol1+str(atom_tag1)+"-"+atom_symbol2+str(atom_tag2)+"-"+atom_symbol3+str(atom_tag3)+" bond ["+str(index1)+"-"+str(index2)+"-"+str(index3)+"]"
  
  xdata=[]
  ydata=[]

  if timestep is None:
    xlab = "MD Steps"
    xUnit = "MD Step"
  else:  
    xlab = "Time [fs]"
    xUnit = "picosecond"

  if len(trajectory) > 1:

    for n, atoms in enumerate(trajectory):
      if timestep is None:
        xdata.append(n)
      else:
        xdata.append(n*timestep/fs)

      ang = atoms.get_angle(index1,index2,index3)

      ydata.append(ang)

    val,err,fit_func = mplot.fit(xdata,ydata,rank=0,verbosity=verbosity-1,fitMin=fitMin,fitMax=fitMax,title=bond_title,yUnit=yUnit,xUnit=xUnit)

  else:
    val = trajectory[0].get_angle(index1,index2,index3)
    err = None

  mout.varOut(bond_title,val,error=err,unit=yUnit,valCol=mcol.result,dataFile=dataFile,verbosity=verbosity)

  if returnData:
    return val, err, bond_title, xdata, ydata
  else:
    return val, err, bond_title

def getCentreOfMass(atoms,verbosity=1):

  positions = getPositions(atoms)

  this_com = [sum([pos[0] for pos in positions])/len(positions),
              sum([pos[1] for pos in positions])/len(positions),
              sum([pos[2] for pos in positions])/len(positions)]

  if verbosity > 0:
    mout.varOut("CoM ("+str(atoms)+")",this_com,unit="Å")

  return this_com

def getRMSD(atoms,reference=None,verbosity=1):

  import mout

  if reference is None:
    CoM = getCentreOfMass(atoms,verbosity=verbosity-1)
  else:
    ref_positions = getPositions(reference)

  positions = getPositions(atoms)

  num_atoms = len(positions)

  big_sum = 0.0

  for index,pos in enumerate(positions):
    if reference is None:
      d = getDisplacement(pos,CoM)
    else:
      d = getDisplacement(pos,ref_positions[index])

    big_sum += pow(d,2)

  rmsd = pow(big_sum/num_atoms,0.5)

  mout.varOut("RMSD ("+str(atoms)+")",rmsd,unit="Å")

  return rmsd

def getDisplacement(position1,position2):
  x = position1[0] - position2[0] 
  y = position1[1] - position2[1] 
  z = position1[2] - position2[2] 
  d = pow(pow(x,2) + pow(y,2) + pow(z,2),0.5)
  return d

def getAngle(position1,position2,position3):
  import numpy as np

  a = np.array(position1)
  b = np.array(position2)
  c = np.array(position3)

  ba = a - b
  bc = c - b

  cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
  angle = np.arccos(cosine_angle)

  return np.degrees(angle)

def getPositions(atoms):

  from ase import Atoms

  if isinstance(atoms,Atoms):
    # ASE ATOMS
    positions = atoms.get_positions()
  elif isinstance(atoms,list):
    positions = []
    for thing in atoms: 
      positions += thing.positions
  else:
    # AMP System/Chain/Residue
    positions = atoms.positions

  return positions
  
  