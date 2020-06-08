
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

import numpy as np

from ase import units

from .analysis import bondLengthStats
from .analysis import bondAngleStats

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

def graphBondLength(trajectory,indices,printScript=False,show=True,filename=None,fitMin=None,fitMax=None,verbosity=2,timestep=None,title=None,fitOrder=None,yUnit="Angstroms",dataFile=None,ymin=0,ymax=None,xmin=None,xmax=None):
  """
    Graph the bond lengths (displacement) between atoms.

  """

  many = any(isinstance(el,list) for el in indices)

  ydata=[]
  labels=[]

  if timestep is None:
    xlab = "MD Steps"
    xUnit = "MD Step"
  else:  
    xlab = "Time [fs]"
    xUnit = "picosecond"

  if many:
    
    for i, pair in enumerate(indices):
      
      if verbosity > 2:
        print("")
      val, err, label, xdata, this_data = bondLengthStats(trajectory,pair,printScript=printScript,verbosity=verbosity-2,timestep=timestep,yUnit=yUnit,returnData=True)

      labels.append(label)
      ydata.append(this_data)

    if fitOrder is None:
      mplot.graph2D(xdata,ydata,ytitles=labels,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=title,verbosity=verbosity-1)
    else:
      if verbosity > 1:
        print("")
      if len(xdata) < fitOrder+1:
        mout.warningOut("Not enough data-points for fitting.")
        fitOrder = None
      if title is not None:
        dataFile.write('# ')
        dataFile.write(title)
        dataFile.write('\n')
      val,err,fit_func = mplot.fit(xdata,ydata,rank=fitOrder,verbosity=verbosity-1,title=title,fitMin=fitMin,fitMax=fitMax,yUnit=yUnit,xUnit=xUnit,dataFile=dataFile)
      text = mplot.getCoeffStr(val,err,1,yUnit=yUnit,xUnit=xUnit)
      mplot.graph2D(xdata,ydata,fitFunc=fit_func,ytitles=labels,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=title,verbosity=verbosity,subtitle=text,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax)

  else:

    if verbosity > 2:
      print("")
    val, err, label, xdata, ydata = bondLengthStats(trajectory,indices,printScript=printScript,verbosity=verbosity-2,timestep=timestep,yUnit=yUnit,returnData=True)

    if fitOrder is None:
      mplot.graph2D(xdata,ydata,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=label,verbosity=verbosity-1)
    else:
      if verbosity > 1:
        print("")
      if len(xdata) < fitOrder+1:
        mout.warningOut("Not enough data-points for fitting.")
        fitOrder = None
      if title is not None:
        dataFile.write('#')
        dataFile.write(title)
        dataFile.write('\n')
      val,err,fit_func = mplot.fit(xdata,ydata,rank=fitOrder,verbosity=verbosity-1,fitMin=fitMin,fitMax=fitMax,title=title,yUnit=yUnit,xUnit=xUnit,dataFile=dataFile)
      text = mplot.getCoeffStr(val,err,1,yUnit=yUnit,xUnit=xUnit)
      mplot.graph2D(xdata,ydata,fitFunc=fit_func,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=label,verbosity=verbosity,subtitle=text,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax)

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

def graphBondAngle(trajectory,indices,printScript=False,show=True,filename=None,fitMin=None,fitMax=None,verbosity=2,timestep=None,title=None,fitOrder=None,yUnit="Degrees",dataFile=None,ymin=0,ymax=180,xmin=None,xmax=None):
  """
    Graph the bond angle (acute) between atoms.

  """

  many = any(isinstance(el,list) for el in indices)

  ydata=[]
  labels=[]

  if timestep is None:
    xlab = "MD Steps"
    xUnit = "MD Step"
  else:  
    xlab = "Time [fs]"
    xUnit = "picosecond"

  if many:
    
    for i, triple in enumerate(indices):
      
      if verbosity > 2:
        print("")
      val, err, label, xdata, this_data = bondAngleStats(trajectory,triple,printScript=printScript,verbosity=verbosity-2,timestep=timestep,yUnit=yUnit,returnData=True)

      labels.append(label)
      ydata.append(this_data)

    if fitOrder is None:
      mplot.graph2D(xdata,ydata,ytitles=labels,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=title,verbosity=verbosity-1)
    else:
      if verbosity > 1:
        print("")
      if len(xdata) < fitOrder+1:
        mout.warningOut("Not enough data-points for fitting.")
        fitOrder = None
      if title is not None:
        dataFile.write('#')
        dataFile.write(title)
        dataFile.write('\n')
      val,err,fit_func = mplot.fit(xdata,ydata,rank=fitOrder,verbosity=verbosity-1,title=title,fitMin=fitMin,fitMax=fitMax,yUnit=yUnit,xUnit=xUnit,dataFile=dataFile)
      text = mplot.getCoeffStr(val,err,1,yUnit=yUnit,xUnit=xUnit)
      mplot.graph2D(xdata,ydata,fitFunc=fit_func,ytitles=labels,show=show,xlab=xlab,ylab="Angle [Degrees]",filename=filename,title=title,verbosity=verbosity,subtitle=text,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax)

  else:

    if verbosity > 2:
      print("")
    val, err, label, xdata, ydata = bondAngleStats(trajectory,indices,printScript=printScript,verbosity=verbosity-2,timestep=timestep,yUnit=yUnit,returnData=True)

    if fitOrder is None:
      mplot.graph2D(xdata,ydata,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=label,verbosity=verbosity-1)
    else:
      if verbosity > 1:
        print("")
      if len(xdata) < fitOrder+1:
        mout.warningOut("Not enough data-points for fitting.")
        fitOrder = None
      if title is not None:
        dataFile.write('# ')
        dataFile.write(title)
        dataFile.write('\n')
      val,err,fit_func = mplot.fit(xdata,ydata,rank=fitOrder,verbosity=verbosity-1,fitMin=fitMin,fitMax=fitMax,title=title,yUnit=yUnit,xUnit=xUnit,dataFile=dataFile)
      text = mplot.getCoeffStr(val,err,1,yUnit=yUnit,xUnit=xUnit)
      mplot.graph2D(xdata,ydata,fitFunc=fit_func,show=show,xlab=xlab,ylab="Distance [Angstrom]",filename=filename,title=label,verbosity=verbosity,subtitle=text,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax)

  if fitOrder is not None:
    return val, err, fit_func
  else:
    return None, None, None
