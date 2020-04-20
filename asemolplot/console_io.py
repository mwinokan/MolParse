
import mcol # https://github.com/mwinokan/MPyTools
from .output import varOut

def printEnergy(atoms,perAtom=True,precision=8,printScript=False):
  if (perAtom):
    epot = atoms.get_potential_energy() / len(atoms)
    ekin = atoms.get_kinetic_energy() / len(atoms)
  else:
    epot = atoms.get_potential_energy()
    ekin = atoms.get_kinetic_energy()
  
  varOut("E_pot",epot,unit="eV/atom",valCol=mcol.result,precision=precision,printScript=printScript,end=", ")
  varOut("E_kin",ekin,unit="eV/atom",valCol=mcol.result,precision=precision,printScript=False)
