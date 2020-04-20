
import mout # https://github.com/mwinokan/MPyTools
import mcol # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

from ase.io.trajectory import Trajectory

def pdb2traj(input,output,verbosity=1,printScript=False): # move to AMP
  if (verbosity > 0):
    mout.out("Converting "+
             mcol.file+input+
             mcol.clear+" to "+
             mcol.file+output+
             mcol.clear+"...",
             printScript=printScript,end='')
    if (verbosity > 1):
      mout.out("")

  in_traj = amp.read(input,index=":",verbosity=verbosity-1)
  amp.write(output,in_traj)

  if (verbosity > 0):
    mout.out("Done.")

def fileLength(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1
  