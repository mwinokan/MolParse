from ase import io

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

def makePovImage(filename,image,**parameters):
  mout.out("processing "+mcol.file+
           filename+".pov"+
           mcol.clear+" ... ",
           printScript=True,
           end='') # user output
  io.write(filename+'.pov',image,
    run_povray=True,
    transparent=True,
    camera_type='perspective',
    **parameters)
  mout.out("Done.") # user output
