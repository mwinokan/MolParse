from ase import io

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

import os

def makePovImage(filename,image,**style):
  mout.out("processing "+mcol.file+
           filename+".pov"+
           mcol.clear+" ... ",
           printScript=True,
           end='') # user output
  if (not style['drawCell']):
      image.set_cell([0,0,0])
  del style['drawCell']
  io.write(filename+'.pov',image,
    run_povray=True,
    transparent=True,
    camera_type='perspective',
    **style)
  mout.out("Done.") # user output

def makePovImages(filename,subdirectory="pov",interval=1,**style):
  os.system("mkdir -p "+subdirectory)
  os.system("rm "+subdirectory+"/*")

  for n, image in enumerate(io.read(filename,index=":")):
    if (n % interval != 0 and n != 100):
      continue

    makePovImage(subdirectory+"/"+str(n).zfill(4),image,**style)

  os.system("rm "+subdirectory+"/*.ini")
  os.system("rm "+subdirectory+"/*.pov")

def makePovAnimation(filename,subdirectory="pov",interval=1,**style):

  import module # https://github.com/mwinokan/MPyTools

  makePovImages(filename,subdirectory=subdirectory,interval=interval,**style)

  module.module('--expert','load','ImageMagick/7.0.3-1-intel-2016a')
  mout.out("ImageMagick loaded.",printScript=True)

  mout.out("creating "+mcol.file+
         "animation.gif"+
         mcol.clear+" ... ",
         printScript=True,
         end='') # user output
  os.system("convert -delay 10 "+subdirectory+"/*.png -fill white -opaque none -loop 1 "+subdirectory+".gif")
  mout.out("Done.") # user output
