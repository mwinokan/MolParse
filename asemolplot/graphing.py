
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

import numpy as np

"""

  To-Do's
    * Time on x-axis
    * graphVelocity
    * graphBondLength

"""

def graphEnergy(trajectory,perAtom=True,filename=None,show=True,verbosity=2):

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
    xdata.append(n)
    epot = atoms.get_potential_energy()
    ekin = atoms.get_kinetic_energy()

    if perAtom:
      epot /= len(atoms)
      ekin /= len(atoms)
      ylab = "Energy eV/atom"
    else:
      ylab = "Energy eV"

    epots.append(epot)
    ekins.append(ekin)
    etots.append(epot+ekin)

  mplot.graph2D(xdata,[epots,ekins,etots],ytitles=["Potential","Kinetic","Total"],show=show,xlab="MD Steps",ylab=ylab,filename=filename,verbosity=verbosity-1)

  if (verbosity > 0):
    mout.out("Done.") # user output

def graphDisplacement(trajectory,show=True,filename=None,relative=True,verbosity=2):
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

    xdata.append(n) # set x-axis to array index

    positions = atoms.get_positions()
    
    if relative:
      if (n==0):
        reference=positions.copy()
      positions -= reference

    rmsd.append(np.sqrt(np.mean(positions**2)))

  mplot.graph2D(xdata,rmsd,show=show,xlab="MD Steps",ylab="RMS Displacement",filename=filename,verbosity=verbosity-1)

  if (verbosity > 0):
    mout.out("Done.") # user output

def showFigs(verbosity=1):
  mplot.show(verbosity=verbosity)
