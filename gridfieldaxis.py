import re
import os
import sys
import math
import glob
import numpy
import os.path
import warnings
import argparse
import subprocess

from scipy import cluster

from sklearn.cluster import KMeans

sys.path.append("./common")
import gridfield

CONVERTER = 1.889716
MINDIM = 20

###############################################################################

def bufcount(filename):
    f = open(filename)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

###############################################################################

def remove_equal(input):
  output = []
  for x in input:
    if x not in output:
      output.append(x)
                        
  return output

###############################################################################

def write_to_cube (mol1, mol1field, fname, xnstep, ynstep, znstep,\
    step, xmin, ymin, zmin):

  opf = open(fname, "w")

  print("Writing final cube... ")

  zero = 0.0
  opf.write("El field\n")
  opf.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (len(mol1.atoms), CONVERTER*xmin, \
      CONVERTER*ymin, CONVERTER*zmin))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (xnstep, CONVERTER*step, zero, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (ynstep, zero, CONVERTER*step, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (znstep, zero, zero, CONVERTER*step))

  for atom in mol1:
    x = CONVERTER*atom.coords[0]
    y = CONVERTER*atom.coords[1]
    z = CONVERTER*atom.coords[2]
    c = atom.partialcharge
    anu = atom.atomicnum

    opf.write("%4d %11.6f %11.6f %11.6f %11.6f\n" % (anu, c, x, y, z))

  for ix in range(xnstep):
    for iy in range(ynstep):
      for iz in range(znstep):
        opf.write("%g "%mol1field[ix,iy,iz])
        if (iz % 6 == 5):
          opf.write("\n")
      opf.write("\n")

  opf.close()


###############################################################################

def get_points(energy, STEPVAL, xmin, ymin, zmin, axis="x", \
  minimaselection=0.0):

  xsets = set()
  ysets = set()
  zsets = set()

  evalset = []

  xvals = []
  yvals = []
  zvals = []

  ixyz_to_xyzval_map = {}
  xyzval_to_ixyz_map = {}

  nx = energy.shape[0]
  ny = energy.shape[1]
  nz = energy.shape[2]

  for iz in range(0, nz):
    z = zmin + float(iz) * (STEPVAL)
    for ix in range(0, nx):
      x = xmin + float(ix) * (STEPVAL)
      for iy in range(0, ny):
        y = ymin + float(iy) * (STEPVAL)

        xsets.add(x)
        ysets.add(y)
        zsets.add(z)

        xvals.append(x)
        yvals.append(y)
        zvals.append(z)

        e = energy[ix, iy, iz]

        if not ( e in evalset):
          evalset.append(e)

        ixyzstr = str(ix)+'_'+str(iy)+'_'+str(iz)
        xstr = "{:.3f}".format(x)
        ystr = "{:.3f}".format(y)
        zstr = "{:.3f}".format(z)
        xyzstr = xstr+'_'+ystr+'_'+zstr
        #print ixyzstr , ' ==> ', xyzstr
        ixyz_to_xyzval_map.update({ixyzstr: xyzstr})
        xyzval_to_ixyz_map.update({xyzstr: ixyzstr})

  dx = sorted(xsets)[1] - sorted(xsets)[0]
  dy = sorted(ysets)[1] - sorted(ysets)[0]
  dz = sorted(zsets)[1] - sorted(zsets)[0]

  botx = min(list(xsets))
  boty = min(list(ysets))
  botz = min(list(zsets))

  topx = max(list(xsets))
  topy = max(list(ysets))
  topz = max(list(zsets))

  print("  Box ")
  print("           x: {:.3f}".format(botx), " {:.3f}".format(topx))
  print("           y: {:.3f}".format(boty), " {:.3f}".format(topy))
  print("           z: {:.3f}".format(botz), " {:.3f}".format(topz))

  print("  field sizes: ", len(xsets), len(ysets), len(zsets))
  print("           dx: {:.3f}".format(dx))
  print("           dy: {:.3f}".format(dy))
  print("           dz: {:.3f}".format(dz))

###############################################################################

if __name__ == "__main__":

  probe = "DRY"
  STEPVAL = 2.0
  DELTAVAL = 10.0
  axis = "x"

  parser = argparse.ArgumentParser()
  parser.add_argument("-f","--file", help="input the mol2 list and weights file", \
      required=True, default="", type=str)
  parser.add_argument("-o","--output", help="output filename ", \
      required=False, default="mean.kont", type=str)
  parser.add_argument("-p","--probe", help="the probe to be used [default="+probe+"]", \
    required=False, default=probe, type=str)
  parser.add_argument("-s","--stepval", help="input stepval value [defaul="+str(STEPVAL)+"]", \
    required=False, default=STEPVAL, type=float)
  parser.add_argument("-d","--deltaval", help="input deltaval value [defaul="+str(DELTAVAL)+"]", \
    required=False, default=DELTAVAL, type=float)
  parser.add_argument("-a","--axis", help="Specify the axis [defaul="+ \
    axis+"]", required=False, default=axis, type=str)

  args = parser.parse_args()

  probe = args.probe
  STEPVAL = args.stepval
  DELLTAVAL = args.deltaval
  axis = args.axis

  if not os.path.isfile("./grid"):
      print("we need grid executable in the current dir")
      exit(1)

  if not os.path.isfile("./fixpdb"):
      print("we need fixpdb executable in the current dir")
      exit(1)

  print("Computing grid fields ...")

  energy, xmin, ymin, zmin = gridfield.compute_grid_mean_field (args.file, \
          STEPVAL, DELTAVAL, probe)

  gridfield.energytofile (energy, args.output, xmin, ymin, zmin, STEPVAL)

  get_points(energy, STEPVAL, xmin, ymin, zmin, args.axis, \
    minimaselection=0.0)

  print("Done ")
