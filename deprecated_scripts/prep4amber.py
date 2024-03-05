#!/usr/bin/env python3

import mout
import argparse
import molparse as amp

argparser = argparse.ArgumentParser(description='Prepare a PDB file for an Amber pipeline')
argparser.add_argument("-i", "--input", help="Input filename")
argparser.add_argument("-o", "--output", help="Output filename")
args = argparser.parse_args()

if args.input is None:
    mout.errorOut("No input given", fatal=True)
if args.output is None:
    mout.errorOut("No output given", fatal=True)

system = amp.parsePDB(args.input)
amp.prep4amber(system)
amp.write(args.output, system)
