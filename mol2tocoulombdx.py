import os
import sys
import glob
import numpy
import argparse
import subprocess

from gridData import Grid

import mendeleev

sys.path.append("./common")
import carbo

psizepath = "./psize.py"
apbspath = "/usr/bin/apbs"

if __name__ == "__main__":

    STEPVAL = 1.4
    DELTAVAL = 20.0
    coulombconst = 1.0
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input the mol2 file", \
        required=True, default="", type=str)
    parser.add_argument("-s","--stepvalue", help="input step value defaul="+str(STEPVAL), \
      required=False, default=STEPVAL, type=float)
    parser.add_argument("-d","--deltaval", help="input delta value defaul="+str(DELTAVAL), \
      required=False, default=DELTAVAL, type=float)
    parser.add_argument("-v","--verbose", help="Verbose mode on", \
      required=False, default=False, action="store_true")
    parser.add_argument("--show", help="Show plot", \
      required=False, default=False, action="store_true")

    args = parser.parse_args()

    basename = os.path.splitext(args.file)[0] 

    mols = carbo.mol2atomextractor(args.file, False)

    cfields = carbo.get_cfields(mols, args.stepvalue, args.deltaval, \
      coulombconst, args.verbose)

    mep = numpy.zeros(cfields[0][1].shape)
    coords = cfields[0][0]
    for field in cfields:
      ep = field[1]
      mep += ep

    mep = mep/float(len(cfields))
    g = Grid(mep, origin=coords[0], \
      delta=[args.stepvalue, args.stepvalue, args.stepvalue])
    name = basename + "_coulomb_mean.dx"
    g.export(name)

    for i, field in enumerate(cfields):
      coords = field[0]
      ep = field[1]
      g = Grid(numpy.asarray(ep), origin=coords[0], \
        delta=[args.stepvalue, args.stepvalue, args.stepvalue])
      
      name = basename + "_coulomb_" + str(i+1) + ".dx"

      g.export(name)