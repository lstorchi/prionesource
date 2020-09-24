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
    parser.add_argument("-f","--file", help="input the mol2 list and weights file", \
        required=True, default="", type=str)
    parser.add_argument("-s","--stepvalue", help="input step value defaul="+str(STEPVAL), \
      required=False, default=STEPVAL, type=float)
    parser.add_argument("-d","--deltaval", help="input delta value defaul="+str(DELTAVAL), \
      required=False, default=DELTAVAL, type=float)
    parser.add_argument("-v","--verbose", help="Verbose mode on", \
      required=False, default=False, action="store_true")
    parser.add_argument("--show", help="Show plot", \
      required=False, default=False, action="store_true")
    parser.add_argument("--ddielectric", help="Enable Dielectric damping", \
      required=False, default=False, action="store_true")
    parser.add_argument("--interpolate", help="Use input DX file to interpolte inm the same grid", \
      required=False, default="", type=str)

    args = parser.parse_args()

    basename = os.path.splitext(args.file)[0] 

    mols = []
    names = []
    weights = []
    fp = open(args.file, "r")
    for l in fp:
        name, ws = l.split()
        w = float(ws)
        mol = carbo.mol2atomextractor(name, False)
        if len(mol) != 1:
            print("Error each file should have 1 molecule")
            exit(1)
        mols.append(mol[0])
        weights.append(w)
        names.append(name)
    fp.close()
    sum = numpy.sum(weights)
    weights /= sum
    #print(numpy.sum(weights))
    #mols = carbo.mol2atomextractor(args.file, False)

    cfields = carbo.get_cfields(mols, args.stepvalue, args.deltaval, \
      coulombconst, args.verbose, args.ddielectric)

    mep = numpy.zeros(cfields[0][1].shape)
    coords = cfields[0][0]
    for i, field in enumerate(cfields):
      ep = field[1]
      mep += ep * weights[i]

    mep = mep/float(len(cfields))
    g = Grid(mep, origin=coords[0], \
      delta=[args.stepvalue, args.stepvalue, args.stepvalue])
    name = basename + "_coulomb_mean.dx"
    if args.ddielectric:
      name = basename + "_coulomb_ddieletric_mean.dx"
    g.export(name)

    tofitw = None
    if args.interpolate != "":
      tofitw = Grid(args.interpolate)

    for i, field in enumerate(cfields):
      coords = field[0]
      ep = field[1]
      g = Grid(numpy.asarray(ep), origin=coords[0], \
        delta=[args.stepvalue, args.stepvalue, args.stepvalue])
      
      basename = os.path.splitext(names[i])[0]
      name = basename + "_coulomb.dx"
      if args.ddielectric:
        name = basename + "_coulomb_ddieletric.dx"
      print (name + " %8.3f"%(weights[i]), file=sys.stderr)

      #name = basename + "_coulomb_" + str(i+1) + ".dx"
      #if args.ddielectric:
      #  name = basename + "_coulomb_ddieletric_" + str(i+1) + ".dx"
      
      if tofitw != None:
        g_on = g.resample(tofitw)
        g_on.export(name)
      else:
        g.export(name)
