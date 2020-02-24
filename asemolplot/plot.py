from ase import io

import mcol # https://github.com/mwinokan/MPyTools
import mout # https://github.com/mwinokan/MPyTools

from . import styles

def makeImage(filename,image,verbosity=1,**style):
  if (verbosity > 0):
    mout.out("creating "+mcol.file+
             filename+".png"+
             mcol.clear+" ... ",
             printScript=True,
             end='') # user output

  if (not style['drawCell']):
      image.set_cell([0,0,0])
      style['show_unit_cell'] = 0
  del style['drawCell']

  # delete povray specific style parameters
  if 'canvas_width' in style:
    style['maxwidth'] = style['canvas_width']
    del style['canvas_width']
  if 'canvas_height' in style:
    del style['canvas_height']
  if 'crop_xshift' in style:
    del style['crop_xshift']
  if 'crop_yshift' in style:
    del style['crop_yshift']
  if 'celllinewidth' in style:
    del style['celllinewidth']
  if 'transparent' in style:
    del style['transparent']

  io.write(filename+'.png',image,
    **style)
  
  if (verbosity > 0):
    mout.out("Done.") # user output

def makeImages(filename,subdirectory="amp",interval=1,verbosity=1,**style):
  
  import os
  
  os.system("mkdir -p "+subdirectory)
  os.system("rm "+subdirectory+"/* 2> /dev/null")

  for n, image in enumerate(io.read(filename,index=":")):
    if (n % interval != 0 and n != 100):
      continue

    makeImage(subdirectory+"/"+str(n).zfill(4),image,verbosity=verbosity-1,**style)

# Using imageio
# https://stackoverflow.com/questions/753190/programmatically-generate-video-or-animated-gif-in-python
def makeAnimation(filename,subdirectory="amp",interval=1,gifstyle=styles.gif_standard,verbosity=1,**plotstyle):
  
  import os
  import imageio

  # Check if a crop is desired
  if "canvas_height" in plotstyle:
    cropping=True
    crop_w = plotstyle["canvas_width"]
    crop_h = plotstyle["canvas_height"]
    del plotstyle["canvas_height"]
  
  if "crop_xshift" in plotstyle:
    crop_x = plotstyle["crop_xshift"]
    del plotstyle["crop_xshift"]
  else:
    crop_x = 0

  if "crop_yshift" in plotstyle:
    crop_y = plotstyle["crop_yshift"]
    del plotstyle["crop_yshift"]
  else:
    crop_y = 0
  
  if "background" in gifstyle:
    if (gifstyle["background"] == "white"):
      backwhite=True
    del gifstyle["background"]

  # Generate the PNG's
  makeImages(filename,subdirectory=subdirectory,interval=interval,verbosity=verbosity-1,**plotstyle)

  # Load ImageMagick
  if cropping:
    import module # https://github.com/mwinokan/MPyTools
    module.module('--expert','load','ImageMagick/7.0.3-1-intel-2016a')
    if (verbosity > 0):
      mout.out("ImageMagick loaded.",printScript=True)

  # Combine the images
  if (verbosity > 0):
    mout.out("loading "+mcol.file+
           subdirectory+"/*.png"+
           mcol.clear+" ... ",
           printScript=True,
           end='') # user output
  images = []
  for file in os.listdir(subdirectory):
    filename = subdirectory+"/"+file
    if file.endswith(".png"):
      if cropping:
        os.system("convert "+filename+
            " -crop "+str(crop_w)+
            "x"+str(crop_h)+
            "+"+str(crop_x)+
            "+"+str(crop_y)+
            " "+filename)
      if backwhite:
        os.system("convert "+filename+
                  " -fill white -opaque none "+
                  filename)
      image = imageio.imread(filename)
      images.append(image)
  if (verbosity > 0):
    mout.out("Done.") # user output

  if (verbosity > 0):
    mout.out("creating "+mcol.file+
           subdirectory+".gif"+
           mcol.clear+" ... ",
           printScript=True,
           end='') # user output
  imageio.mimsave(subdirectory+".gif",images,**gifstyle)

  if (verbosity > 0):
    mout.out("Done.") # user output
