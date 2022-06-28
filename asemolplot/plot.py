
def makeImage(filename,image,filter=None,verbosity=1,printScript=False,**style):
  import mcol
  import mout
  from ase import io

  if not filename.endswith(".png"):
    filename = filename+".png"

  if (verbosity > 0):
    mout.out("creating "+mcol.file+
             filename+
             mcol.clear+" ... ",
             printScript=printScript,
             end='') # user output

  if not 'drawCell' in style:
    style['drawCell'] = False
  if not style['drawCell']:
    image.set_cell([0,0,0])
    style['show_unit_cell'] = 0
    # style['celllinewidth'] = 0
  if 'drawCell' in style:
    del style['drawCell']

  # print(style)

  # delete povray specific style parameters
  if 'canvas_width' in style:
    style['maxwidth'] = style['canvas_width']
    del style['canvas_width']
  if 'canvas_height' in style:
    del style['canvas_height']
  # if 'crop_xshift' in style:
  #   del style['crop_xshift']
  # if 'crop_yshift' in style:
  #   del style['crop_yshift']
  if 'celllinewidth' in style:
    del style['celllinewidth']
  # if 'transparent' in style:
  #   del style['transparent']

  io.write(filename,image,
    **style)
  
  if (verbosity > 0):
    mout.out("Done.") # user output

def makeImages(filename,subdirectory="amp",interval=1,verbosity=1,filenamePadding=4,printScript=False,index=":",**style):
  from ase import io
  import os
  import math
  import mcol
  import mout
  
  os.system("mkdir -p "+subdirectory)
  os.system("rm "+subdirectory+"/* 2> /dev/null")

  if index != ":":
    image = io.read(filename,index=index)
    makeImage(subdirectory+"/"+str(index).zfill(filenamePadding),image,verbosity=verbosity-1,**style)
  else:
    traj = io.read(filename,index=index)

    global num_traj_images
    global num_frames
    num_traj_images = len(traj)
    num_frames = math.ceil(num_traj_images/interval)
    if verbosity > 0:
      mout.varOut("Trajectory",num_traj_images,unit="images")
      mout.varOut("Animation",num_frames,unit="frames")
    
    for n, image in enumerate(traj):
      if (n % interval != 0 and n != 100):
        continue

      if verbosity == 1:
        mout.progress(n+1,num_traj_images,prepend="Creating images",printScript=printScript)
      
      makeImage(subdirectory+"/"+str(n).zfill(filenamePadding),image,verbosity=verbosity-1,printScript=printScript,**style)

# Using imageio
# https://stackoverflow.com/questions/753190/programmatically-generate-video-or-animated-gif-in-python
def makeAnimation(filename,subdirectory="amp",interval=1,gifstyle=None,verbosity=1,printScript=False,useExisting=False,dryRun=False,**plotstyle):
  import mcol
  import mout

  from . import styles
  if gifstyle is None:
    gifstyle = styles.gif_standard

  if plotstyle == {}:
    plotstyle=styles.standard

  global num_traj_images
  global num_frames

  import os
  import imageio

  cropping = False
  shifting = False

  # Set canvas sizes:
  canv_w = plotstyle["canvas_width"]
  if "canvas_height" in plotstyle:
    canv_h = plotstyle["canvas_height"]
    del plotstyle["canvas_height"]
  else:
    mout.errorOut("No canvas_height specified in plotstyle!",fatal=True)

  # Check if cropping:
  if "crop_w" in plotstyle and "crop_w" in plotstyle:
    cropping = True
    crop_w = plotstyle["crop_w"]
    crop_h = plotstyle["crop_h"]
  if "crop_w" in plotstyle: del plotstyle["crop_w"]
  if "crop_h" in plotstyle: del plotstyle["crop_h"]

  # Check if crop offset:
  if "crop_x" in plotstyle:
    shifting = True
    crop_x = plotstyle["crop_x"]
    crop_y = plotstyle["crop_y"]
  if "crop_x" in plotstyle: del plotstyle["crop_x"]
  if "crop_y" in plotstyle: del plotstyle["crop_y"]

  # Generate the PNG's
  if not useExisting:
    if (verbosity > 0):
      mout.out("generating "+mcol.file+
             subdirectory+"/*.png"+
             mcol.clear+" ... ",
             printScript=printScript,end='') # user output
    if (verbosity > 1):
      mout.out(" ")

    if not dryRun:
      # Generate all the images
      makeImages(filename,subdirectory=subdirectory,interval=interval,verbosity=verbosity-1,**plotstyle)
    else:
      # Generate just the first image
      makeImages(filename,subdirectory=subdirectory,interval=interval,verbosity=verbosity-1,index=0,**plotstyle)
    
    if (verbosity == 1):
      mout.out("Done.")

  # Load ImageMagick
  # if cropping or backwhite:
  import module # https://github.com/mwinokan/MPyTools
  ret = module.module('--expert','load','ImageMagick/7.0.3-1-intel-2016a')
  if ret == 0 and verbosity > 0:
    mout.out("ImageMagick loaded.",printScript=printScript)

  # Combine the images
  if (verbosity > 0):
    mout.out("loading "+mcol.file+
           subdirectory+"/*.png"+
           mcol.clear+" ... ",
           printScript=printScript,
           end='') # user output

  images = []

  # print(os.listdir(subdirectory))

  # loop over all files in the subdirectory:
  for file in sorted(os.listdir(subdirectory)):

    # get the relative path to the file:
    filename = subdirectory+"/"+file

    # check if the file is a PNG:
    if file.endswith(".png"):

      # run different IM commands depending on cropping and shifting:
      if not cropping and not shifting:
        os.system("convert "+filename+
                  " -background white -extent "+
                  str(canv_w)+"x"+
                  str(canv_h)+" "+
                  filename)
      elif cropping and not shifting:
        os.system("convert "+filename+
                  " -crop "+str(crop_w)+"x"+str(crop_h)+
                  " -background white -extent "+
                  str(crop_w)+"x"+
                  str(crop_h)+" "+
                  filename)
        # print("convert "+filename+
        #           " -crop "+str(crop_w)+"x"+str(crop_h)+
        #           " -background white -extent "+
        #           str(crop_w)+"x"+
        #           str(crop_h)+" "+
        #           filename)
      elif shifting and not cropping:
        os.system("convert "+filename+
                  " -crop +"+str(crop_x)+"+"+str(crop_y)+
                  " -background white -extent "+
                  str(canv_w)+"x"+
                  str(canv_h)+" "+
                  filename)
      else:
        os.system("convert "+filename+
                  " -crop "+str(crop_w)+"x"+str(crop_h)+
                  "+"+str(crop_x)+"+"+str(crop_y)+
                  " -background white -extent "+
                  str(crop_w)+"x"+
                  str(crop_h)+" "+
                  filename)

      # Read in the image and append to the image array
      image = imageio.imread(filename)
      images.append(image)

      # print(len(images),filename)

      if verbosity > 0 and not dryRun and not useExisting:
        mout.progress(len(images),num_frames,prepend="Cropping & loading images",printScript=printScript)

  if (verbosity > 0) and not useExisting:
    mout.out("Done.") # user output

  if (verbosity > 0):
    mout.out("creating "+mcol.file+
           subdirectory+".gif"+
           mcol.clear+" ... ",
           printScript=printScript,
           end='') # user output

  # Generate the animated GIF:
  imageio.mimsave(subdirectory+".gif",images,**gifstyle)

  if (verbosity > 0):
    mout.out("Done.") # user output
