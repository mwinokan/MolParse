
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

def bondLengthStats(trajectory,index_pair,printScript=False,verbosity=1,timestep=None,yUnit="Angstroms",returnData=False):

  atom_symbols = trajectory[0].get_chemical_symbols()
  atom_tags = trajectory[0].get_tags()

  index1 = index_pair[0]
  index2 = index_pair[1]

  atom_symbol1 = atom_symbols[index1]
  atom_symbol2 = atom_symbols[index2]
  if atom_tags[index1] != 0:
    atom_tag1 = str(atom_tags[index1])
  else:
    atom_tag1 = ""
  if atom_tags[index2] != 0:
    atom_tag2 = str(atom_tags[index2])
  else:
    atom_tag2 = ""

  atom_label1 = atom_symbol1+str(atom_tag1)+"["+str(index1)+"]"
  atom_label2 = atom_symbol2+str(atom_tag2)+"["+str(index2)+"]"
  labels = [atom_label1, atom_label2]
  bond_title = atom_symbol1+str(atom_tag1)+"-"+atom_symbol2+str(atom_tag2)+" bond ["+str(index1)+"-"+str(index2)+"]"

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
        xdata.append(n*timestep/units.fs)

      dist = atoms.get_distance(index1,index2)

      ydata.append(dist)

    val,err,fit_func = mplot.fit(xdata,ydata,rank=0,verbosity=verbosity-1,title=bond_title,yUnit=yUnit,xUnit=xUnit)

  else:

    val = trajectory[0].get_distance(index1,index2)
    err = None

  if verbosity > 0:
    mout.varOut(bond_title,val,error=err,unit=yUnit,valCol=mcol.result)

  if returnData:
    return val, err, bond_title, xdata, ydata
  else:
    return val, err, bond_title


def bondAngleStats(trajectory,index_triplet,printScript=False,verbosity=1,timestep=None,yUnit="degrees",returnData=False):

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
        xdata.append(n*timestep/units.fs)

      ang = atoms.get_angle(index1,index2,index3)

      ydata.append(ang)

    val,err,fit_func = mplot.fit(xdata,ydata,rank=0,verbosity=verbosity-1,title=bond_title,yUnit=yUnit,xUnit=xUnit)

  else:
    val = trajectory[0].get_angle(index1,index2,index3)
    err = None

  mout.varOut(bond_title,val,error=err,unit=yUnit,valCol=mcol.result)

  if returnData:
    return val, err, bond_title, xdata, ydata
  else:
    return val, err, bond_title

def getCentreOfMass(atoms):

  positions = atoms.get_positions()

  this_com = [sum([pos[0] for pos in positions])/len(positions),
              sum([pos[1] for pos in positions])/len(positions),
              sum([pos[2] for pos in positions])/len(positions)]

  print(this_com)
