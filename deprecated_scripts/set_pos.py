#!/usr/bin/env python3

import molparse as amp
import mout
import mcol

##########################################################################

import argparse

argparser = argparse.ArgumentParser(description="Shift the z-coordinate of a single atom")

argparser.add_argument("input", metavar="INPUT", help="Input .gro or .pdb file")
argparser.add_argument("index", metavar="INDEX", help="Atom index to be shifted [0:N-1]")
argparser.add_argument("xpos", metavar="X", help="X-coordinate [nanometres]")
argparser.add_argument("ypos", metavar="Y", help="Y-ccordinate [nanometres]")
argparser.add_argument("zpos", metavar="Z", help="Z-ccordinate [nanometres]")
args = argparser.parse_args()

##########################################################################

if args.input.endswith(".gro"):
    system = amp.parseGRO(args.input)
if args.input.endswith(".pdb"):
    system = amp.parsePDB(args.input)

atom = system.atoms[int(args.index)]

# print(atom.position)
if args.xpos != " ":
    atom.position[0] = float(args.xpos)
if args.ypos != " ":
    atom.position[1] = float(args.ypos)
if args.zpos != " ":
    atom.position[2] = float(args.zpos)

# print(atom.position)

mout.headerOut("Set coordinates of atom " +
               mcol.result + atom.name +
               mcol.clear + " (" +
               mcol.arg + str(args.index) +
               mcol.clear + ") of residue " +
               mcol.result + atom.residue +
               mcol.clear + " to " +
               mcol.arg + str(atom.position) +
               mcol.varType + " nanometres")

# system.summary()

amp.write(args.input, system)
