#!/usr/bin/env python3

"""

amp_anim8.py
------------

Animation with ASE & AMP

Help & Usage: python amp_anim8.py -h

- Max Winokan

[ Part of AseMolPlot 
  https://github.com/mwinokan/AseMolPlot ]

"""

import argparse

import mout # https://github.com/mwinokan/MPyTools
import mcol # https://github.com/mwinokan/MPyTools
import mplot # https://github.com/mwinokan/MPyTools

import asemolplot as amp # https://github.com/mwinokan/AseMolPlot

# from ase.io.trajectory import Trajectory

##########################################################################

argparser = argparse.ArgumentParser(description="Animate multi-model PDB's, and ASE trajectories")

argparser.add_argument("input",metavar="INPUT",help="Input PDB/TRAJ")
argparser.add_argument("-pov","--povray",help="Use PoV-Ray",default=False,action='store_true')
argparser.add_argument("-dr","--dry-run",help="Only generate first image",default=False,action='store_true')
argparser.add_argument("-o","--output",help="Output keyword",type=str)
argparser.add_argument("-ps","--print-script", type=mout.str2bool,nargs='?',const=True,default=False,help="Print the script name in console output.")
argparser.add_argument("-i","--interval", type=int,default=1,help="Animation frame interval")
argparser.add_argument("-rx","--rotate-x",type=float,help="Rotation x")
argparser.add_argument("-ry","--rotate-y",type=float,help="Rotation x")
argparser.add_argument("-rz","--rotate-z",type=float,help="Rotation x")

argparser.add_argument("-cw","--width",type=int,help="Canvas Width")
argparser.add_argument("-ch","--height",type=int,help="Canvas Height")
argparser.add_argument("-c","--crop",nargs=2,help="Crop width and height")
argparser.add_argument("-cs","--crop-shift",nargs=2,help="Crop-shift x and y")

args = argparser.parse_args()

##########################################################################

# process the arguments
infile = args.input
printScript = args.print_script

verbosity = 2

if args.output is None:
  out_prefix = "amp_out"
  mout.warningOut("Defaulted to output keyword 'amp_out'.",printScript=printScript,code=3)
else:
  out_prefix = args.output

# custom_style = amp.styles.standard.copy()
custom_style = amp.styles.standard_cell.copy()

# print(args.rotate_x)
# print(args.rotate_y)
# print(args.rotate_z)

if args.rotate_x is not None or args.rotate_y is not None or args.rotate_z is not None:
  custom_style["rotation"] = ""
  if args.rotate_x is not None:
    custom_style["rotation"] = str(args.rotate_x)+"x"
  if args.rotate_y is not None:
    if custom_style["rotation"] == "":
      custom_style["rotation"] = str(args.rotate_y)+"y"
    else:
      custom_style["rotation"] = custom_style["rotation"]+","+str(args.rotate_y)+"y"
  if args.rotate_z is not None:
    if custom_style["rotation"] == "":
      custom_style["rotation"] = str(args.rotate_z)+"z"
    else:
      custom_style["rotation"] = custom_style["rotation"]+","+str(args.rotate_z)+"z"

custom_style["canvas_width"] = 1200
custom_style["canvas_height"] = 1200

if args.width is not None:
  custom_style["canvas_width"] = args.width
if args.height is not None:
  custom_style["canvas_height"] = args.height

if args.crop is not None:
  custom_style["crop_w"] = args.crop[0]
  custom_style["crop_h"] = args.crop[1]
if args.crop_shift is not None:
  custom_style["crop_x"] = args.crop_shift[0]
  custom_style["crop_y"] = args.crop_shift[1]

if not args.povray:
  custom_style["scale"] = custom_style["canvas_width"]/500 * 20

# print(custom_style["rotation"])

if not infile.endswith(".pdb"):
  mout.warningOut("Unrecognised file type, attempting conversion!")
  thing = amp.read(infile)  
  amp.write("__temp__.pdb",thing)
  infile = "__temp__.pdb"

infile_read = open(infile,'r').read()

model_count = infile_read.count("MODEL")

mout.varOut("Models in PDB",model_count)

if model_count > 1:
  mout.headerOut("Animating...")

  if args.povray:
    amp.makePovAnimation(infile,subdirectory=out_prefix,verbosity=verbosity,interval=args.interval,dryRun=args.dry_run,**custom_style)
  else:
    amp.makeAnimation(infile,subdirectory=out_prefix,verbosity=verbosity,interval=args.interval,dryRun=args.dry_run,**custom_style)

else:
  mout.headerOut("Rendering...")

  atoms = amp.read(infile)

  if args.povray:
    amp.makePovImage(out_prefix,atoms,verbosity=verbosity,**custom_style)
  else:
    amp.makeImage(out_prefix,atoms,verbosity=verbosity,**custom_style)

exit()
