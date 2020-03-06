
import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt
# from ase.io.trajectory import Trajectory

import numpy as np

"""

  To-Do's
    * graph2D (include the graph2DMany functionality and check the shape of ydata to allow lists and lists of lists.)
    * graph2D (check if ytitles is passed)
    * Verbosity

"""

def graphEnergy(trajectory,perAtom=True,filename=None,verbosity=1):

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

  # graph2D(xdata,epots,xlab='index',ylab='Potential Energy eV',filename=filename)

  graph2DMany(xdata,[epots,ekins,etots],ytitles=["Potential","Kinetic","Total"],xlab="MD Steps",ylab=ylab,filename=filename,verbosity=verbosity-1)

  if (verbosity > 0):
    mout.out("Done.") # user output


# def graph2D(xdata,ydata,show=True,xmin=None,xmax=None,ymin=None,ymax=None):
def graph2D(xdata,ydata,filename=None,show=True,xmin=None,xmax=None,ymin=None,ymax=None,xlab='x',ylab='y',title=None):
  plt.plot(xdata,ydata)
  plt.axis([xmin,xmax,ymin,ymax])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.suptitle(title)
  if show:
    plt.show()
  if filename is not None:
    plt.savefig(filename)

def graph2DMany(xdata,ydata,ytitles=None,filename=None,show=True,xmin=None,xmax=None,ymin=None,ymax=None,xlab='x',ylab='y',title=None,verbosity=1):

  for curve, label in zip(ydata,ytitles):
    plt.plot(xdata,curve,label=label)

  plt.axis([xmin,xmax,ymin,ymax])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.suptitle(title)
  plt.legend()
  if show:
    if (verbosity > 0):
      mout.out("showing ... ",end='')
    plt.show()
  if filename is not None:
    if (verbosity > 0):
      mout.out("saving as " + mcol.file + filename + mcol.clear + " ... ",end='')
    plt.savefig(filename)
