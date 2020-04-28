
import mout # https://github.com/mwinokan/MPyTools
import mcol # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

from ase.io.trajectory import Trajectory

from .io import read,write

def pdb2traj(input,output,verbosity=1,printScript=False,tagging=True): # move to AMP
  if not input.endswith(".pdb"):
    mout.errorOut("Input is not a PDB!")
    return None

  if (verbosity > 0):
    mout.out("Converting "+
             mcol.file+input+
             mcol.clear+" to "+
             mcol.file+output+
             mcol.clear+" ... ",
             printScript=printScript,end='')
    if (verbosity > 1):
      mout.out("")

    if tagging:
      # initialise arrays
      taglist=[]

      if (verbosity > 1):
        mout.out("parsing PDB for tags ... ",printScript=printScript,end='')

      # Parse the first model in the PDB to get the tags
      with open(input,"r") as input_pdb:
        searching = True
        for line in input_pdb:
          if searching:
            if line.startswith("MODEL"):
              searching = False
          else:
            if line.startswith("ENDMDL"):
              break
            if line.startswith("TER"):
              continue
            else:
              # get the tags:
              taglist.append(int(''.join(filter(lambda i: i.isdigit(), line.strip().split()[2]))))

      if (verbosity > 1):
        mout.out("Done.")

    in_traj = read(input,index=":",verbosity=verbosity-1)

    if tagging:
      # set the tags
      for atoms in in_traj:
        for index,tag in enumerate(taglist):
          atoms[index].tag = tag

    write(output,in_traj,verbosity=verbosity-1)

    if (verbosity == 1):
      mout.out("Done.")

    return in_traj

def gro2traj(input,output,verbosity=1,printScript=False,tagging=True):
  if not input.endswith(".gro"):
    mout.errorOut("Input is not a .gro file!")
    return None

  if (verbosity > 0):
    mout.out("Converting "+
             mcol.file+input+
             mcol.clear+" to "+
             mcol.file+output+
             mcol.clear+" ... ",
             printScript=printScript,end='')
    if (verbosity > 1):
      mout.out("")

  in_traj = read(input,index=":",verbosity=verbosity-1)

  write(output,in_traj,verbosity=verbosity-1)

  if (verbosity == 1):
    mout.out("Done.")

  return in_traj

def fileLength(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1
  