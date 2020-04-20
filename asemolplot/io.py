from ase import io

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

def write(filename,image,verbosity=1,printScript=False,**parameters):

  if (verbosity > 0):
    mout.out("writing "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  io.write(filename,image,**parameters)

  if (verbosity > 0):
    mout.out("Done.") # user output

def read(filename,index=None,printScript=False,verbosity=1,**parameters):

  if (verbosity > 0):
    mout.out("reading "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  atoms = io.read(filename,index,**parameters)

  if (verbosity > 0):
    mout.out("Done.") # user output

  return atoms
