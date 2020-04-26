
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

import numpy as np

from ase import units

"""

  To-Do's
    * If show = True close all other figures?
    * graphVelocity
    * graphTemperature
    * graphGyration

"""

def graphEnergy(trajectory,perAtom=True,filename=None,show=True,verbosity=2,kJpermol=False,xlab=None,timestep=None):

  if (verbosity > 0):
    mout.out("graphing "+mcol.varName+
             "Energy"+
             mcol.clear+" ... ",
             printScript=True,
             end='') # user output

  xdata=[]
  ekins=[]
  epots=[]
  etots=[]

  for n, atoms in enumerate(trajectory):
    if timestep is None:
      if xlab is None:
        xlab = "MD Steps"
      xdata.append(n)
    else:
      xlab = "Time [fs]"
      xdata.append(n*timestep/units.fs)
    epot = atoms.get_potential_energy()
    ekin = atoms.get_kinetic_energy()

    if not kJpermol:
      if perAtom:
        epot /= len(atoms)
        ekin /= len(atoms)
        ylab = "Energy [eV/atom]"
      else:
        ylab = "Energy [eV]"
    else:
      epot *= units.eV / units.kJ * units.mol
      ekin *= units.eV / units.kJ * units.mol
      ylab = "Energy [kJ/mol]"

    epots.append(epot)
    ekins.append(ekin)
    etots.append(epot+ekin)

  mplot.graph2D(xdata,[epots,ekins,etots],ytitles=["Potential","Kinetic","Total"],show=show,xlab=xlab,ylab=ylab,filename=filename,verbosity=verbosity-1)

  if (verbosity > 0):
    mout.out("Done.") # user output

def graphDisplacement(trajectory,show=True,filename=None,relative=True,verbosity=2,timestep=None):
  """
    Root Mean Square Displacement 

  """

  if (verbosity > 0):
    mout.out("graphing "+mcol.varName+
             "RMSD"+
             mcol.clear+" ... ",
             printScript=True,
             end='') # user output

  xdata=[]
  rmsd=[]

  for n, atoms in enumerate(trajectory):

    if timestep is None:
      xlab = "MD Steps"
      xdata.append(n)
    else:
      xlab = "Time [fs]"
      xdata.append(n*timestep/units.fs)

    positions = atoms.get_positions()
    
    if relative:
      if (n==0):
        reference=positions.copy()
      positions -= reference

    rmsd.append(np.sqrt(np.mean(positions**2)))

  mplot.graph2D(xdata,rmsd,show=show,xlab="MD Steps",ylab="RMS Displacement",filename=filename,verbosity=verbosity-1)

  if (verbosity > 0):
    mout.out("Done.") # user output

def graphBondLength(trajectory,indices,printScript=False,show=True,filename=None,verbosity=2,timestep=None,title=None,fitOrder=None):
  """
    Graph the bond lengths (displacement) between atoms.

  """

  if (verbosity > 0):
    mout.out("graphing "+mcol.varName+
             title+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  many = any(isinstance(el,list) for el in indices)

  xdata=[]
  ydata=[]
  labels=[]

  if many:
    
    atom_symbols = trajectory[0].get_chemical_symbols()
    atom_tags = trajectory[0].get_tags()

    for i, pair in enumerate(indices):
      this_data=[]

      index1 = pair[0]
      index2 = pair[1]

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

      labels.append(atom_symbol1+str(atom_tag1)+"-"+
                    atom_symbol2+str(atom_tag2)+
                    " bond ["+str(index1)+"-"+
                    str(index2)+"]")

      for n,atoms in enumerate(trajectory):
        if (i==0):
          if timestep is None:
            xlab = "MD Steps"
            xUnit = "MD Step"
            xdata.append(n)
          else:
            xlab = "Time [ps]"
            xUnit = "picosecond"
            xdata.append(n*timestep)
            # xdata.append(n*timestep/units.fs)

        dist = atoms.get_distance(index1,index2)

        this_data.append(dist)

      ydata.append(this_data)

    if len(xdata) < fitOrder+1:
      mout.warningOut("Not enough data-points for fitting.")
      fitOrder = None

    if fitOrder is None:
      mplot.graph2D(xdata,ydata,ytitles=labels,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=title,verbosity=verbosity-1)
    else:
      if verbosity > 1:
        print("")
      val,err,fit_func = mplot.fit(xdata,ydata,rank=fitOrder,verbosity=verbosity-1,title=title,yUnit="Angstroms",xUnit=xUnit)
      text = mplot.getCoeffStr(val,err,1,yUnit="Angstroms",xUnit=xUnit)
      mplot.graph2D(xdata,ydata,fitFunc=fit_func,ytitles=labels,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=title,verbosity=verbosity-1,subtitle=text)

  else:
    index1 = indices[0]
    index2 = indices[1]

    atom_symbols = trajectory[0].get_chemical_symbols()
    atom_symbol1 = atom_symbols[index1]
    atom_symbol2 = atom_symbols[index2]
    label= atom_symbol1+atom_symbol2+" bond ["+str(index1)+"-"+str(index2)+"]"

    ydata=[]

    if timestep is None:
      xlab = "MD Steps"
      xUnit = "MD Step"
    else:  
      xlab = "Time [fs]"
      xUnit = "picosecond"

    for n, atoms in enumerate(trajectory):
      if timestep is None:
        xdata.append(n)
      else:
        xdata.append(n*timestep/units.fs)

      dist = atoms.get_distance(index1,index2)

      ydata.append(dist)

    if fitOrder is None:
      mplot.graph2D(xdata,ydata,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=label,verbosity=verbosity-1)
    else:
      if verbosity > 1:
        print("")
      val,err,fit_func = mplot.fit(xdata,ydata,rank=fitOrder,verbosity=verbosity-1,title=title,yUnit="Angstroms",xUnit=xUnit)
      text = mplot.getCoeffStr(val,err,1,yUnit="Angstroms",xUnit=xUnit)
      mplot.graph2D(xdata,ydata,fitFunc=fit_func,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=label,verbosity=verbosity-1,subtitle=text)

  if (verbosity > 1):
    mout.out("") # user output
  if (verbosity == 1):
    mout.out("Done.") # user output

  if fitOrder is not None:
    return val, err, fit_func
  else:
    return None, None, None

def graphBondVibSpec(trajectory,indices,printScript=False,show=True,filename=None,verbosity=2,timestep=None,title=None):
  """
    Graph the fourier transfort of bond lengths (displacement) between atoms.

  """

  # if (verbosity > 0):
  #   mout.out("graphing "+mcol.varName+
  #            title+
  #            mcol.clear+" ... ",
  #            printScript=printScript,
  #            end='') # user output

  many = any(isinstance(el,list) for el in indices)

  xdata=[]
  ydata=[]
  labels=[]

  index1 = indices[0]
  index2 = indices[1]

  for n, atoms in enumerate(trajectory):
    if timestep is None:
      xdata.append(n)

    dist = atoms.get_distance(index1,index2)
    ydata.append(dist)

  ydata = np.fft.fft(ydata)

  mplot.graph2D(xdata,ydata,show=show,filename=filename,verbosity=verbosity-1)

# just a wrapper for mplot.show()
def showFigs(verbosity=1):
  mplot.show(verbosity=verbosity)
