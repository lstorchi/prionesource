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

probe = "DRY"
STEPVAL = 0.5
DELTAVAL = 10.0
MINDIM = 0
NUMOFCLUST = 0

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="input the mol2 list and weights file", \
    required=True, default="", type=str)
parser.add_argument("-p","--probe", help="the probe to be used [default="+probe+"]", \
  required=False, default=probe, type=str)
parser.add_argument("-s","--stepval", help="input stepval value [defaul="+str(STEPVAL)+"]", \
  required=False, default=STEPVAL, type=float)
parser.add_argument("-d","--deltaval", help="input deltaval value [defaul="+str(DELTAVAL)+"]", \
  required=False, default=DELTAVAL, type=float)
parser.add_argument("-m","--numofdim", help="input number of minima to be used [defaul="+ \
  str(MINDIM)+"]", required=False, default=MINDIM, type=int)
parser.add_argument("-c","--numofcluster", help="input number of clusters to be used [defaul=" + \
  str(NUMOFCLUST)+"]", required=False, default=DELTAVAL, type=float)

args = parser.parse_args()

probe = args.probe
MINDIM = args.numofdim
NUMOFCLUST = args.numofcluster
STEPVAL = args.stepval
DELTAVAL = args.stepval


if not os.path.isfile("./grid"):
    print("we nee grid executable in the current dir")
    exit(1)

if not os.path.isfile("./fixpdb"):
    print("we nee fixpdb executable in the current dir")
    exit(1)

print("Computing grid fields ...")

energy, xmin, ymin, zmin = gridfield.compute_grid_avg_field (filename, \
        weightfile, STEPVAL, DELTAVAL, probe)

gridfield.energytofile (energy, "mean.kont", xmin, ymin, zmin, STEPVAL)

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

print("nx: ", nx, " ny: ", ny, " nz: ", nz)

for iz in range(0, nz):
  z = zmin + float(iz) * (1.0/STEPVAL)
  for ix in range(0, nx):
    x = xmin + float(ix) * (1.0/STEPVAL)
    for iy in range(0, ny):
      y = ymin + float(iy) * (1.0/STEPVAL)

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

print("Selecting min values ...")

mineval = []
sortedeval = numpy.sort(evalset)

for i in range(0, min(MINDIM, len(sortedeval))):
  mineval.append(sortedeval[i])

print("Min values selected: ")
i = 1
for m in mineval:
    print(i , " ==> " , m)
    i = i + 1

print("Done...")

print("Start computing ... ")

dimcube = 5

min_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)

for ix in range(0,nx):
  for iy in range(0,ny):
    for iz in range(0,nz):
      min_energy[ix, iy, iz] = 0.0

xmin = []
ymin = []
zmin = []
minvals = []

for ix in range(dimcube,nx-dimcube):
  for iy in range(dimcube,ny-dimcube):
    for iz in range(dimcube,nz-dimcube):

      x = botx + ix * dx
      y = boty + iy * dy
      z = botz + iz * dz

      eref = energy[ix, iy, iz]

      notminima = False
      
      if (eref in mineval):

        for ix_near in range(-dimcube,dimcube):
          for iy_near in range(-dimcube,dimcube):
            for iz_near in range(-dimcube,dimcube):
        
              if not ((ix_near == 0) and (iy_near == 0) and (iz_near == 0)):
        
                x_near = x + (ix_near * dx)
                y_near = y + (iy_near * dy)
                z_near = z + (iz_near * dz)
        
                xstr = "{:.3f}".format(x)
                ystr = "{:.3f}".format(y)
                zstr = "{:.3f}".format(z)
                xyzstr = xstr+'_'+ystr+'_'+zstr
                ixyz = xyzval_to_ixyz_map[xyzstr].split("_")
                near_ix = int(ixyz[0])
                near_iy = int(ixyz[1])
                near_iz = int(ixyz[2])

                e = energy[near_ix, near_iy, near_iz]
                
                if (e < eref):
                  notminima = True
        
        if (notminima):
          min_energy[ix, iy, iz] = 0.0
        else:
          #print "{:.3f}".format(x), " {:.3f}".format(y), " {:.3f}".format(z), " {:.3f}".format(eref)
          min_energy[ix, iy, iz] = eref
          minvals.append(eref)
          xmin.append(x)
          ymin.append(y)
          zmin.append(z)

print("Done ")

print("Start clustering ...")

pointstocluster = numpy.zeros((len(minvals), 3))

for i in range(0,len(minvals)):
  pointstocluster[i,0] = xmin[i]
  pointstocluster[i,1] = ymin[i]
  pointstocluster[i,2] = zmin[i]

with warnings.catch_warnings(record=True) as w:
    centroids, selected = cluster.vq.kmeans2 (pointstocluster, NUMOFCLUST)
    if len(w) != 0:
        print(w[-1].message)
        #exit(1)

numofcluster = 0
for j in range(0, NUMOFCLUST):
    for k in range(len(selected)):
        if (selected[k] == j):
            numofcluster = numofcluster + 1
            break

print("Final num of cluster: ", numofcluster)

centroidvals = []

for j in range(0, NUMOFCLUST):
  centroidvals.append(0.0)

for j in range(0,len(selected)):
  idx = selected[j]
  centroidvals[idx] += minvals[j] 

gridfield.ifextrm("./centroid.xyz")

xyzf = open("./centroid.xyz", "w")

xyzf.write(str(numofcluster)+"\n")
xyzf.write(" \n")

print("Controid")
print("X          Y          Z          EMIN       RMIN       RMAX       RAVG")
for j in range(0, NUMOFCLUST):
  rmin = float('inf')
  rmax = -1.0
  ravg = 0.0
  dim = 0.0

  for k in range(len(selected)):
    if (selected[k] == j):
      x = centroids[j][0] - pointstocluster[k][0]
      y = centroids[j][1] - pointstocluster[k][1]
      z = centroids[j][2] - pointstocluster[k][2]
      r = math.sqrt(x*x + y*y + z*z)

      if r < rmin :
        rmin = r
      if r > rmax :
        rmax = r
     
      ravg = ravg + r
      dim = dim + 1.0

  if dim != 0.0:
    ravg = ravg / dim
            
    xyzf.write(" H %10.5f %10.5f %10.5f\n"%(\
            centroids[j][0], centroids[j][1], \
            centroids[j][2]))
    
    print("%10.5f"%centroids[j][0], \
            "%10.5f"%centroids[j][1], \
            "%10.5f"%centroids[j][2], \
            "%10.5f"%centroidvals[j], \
            "%10.5f"%rmin, \
            "%10.5f"%rmax, \
            "%10.5f"%ravg)

xyzf.close()

"""
centroids_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)

for ix in range(0,nx):
  for iy in range(0,ny):
    for iz in range(0,nz):
      centroids_energy[ix, iy, iz] = 0.0

for j in range(0, NUMOFCLUST):
  ixc = int((centroids[j, 0] - botx) / dx)
  iyc = int((centroids[j, 1] - boty) / dy)
  izc = int((centroids[j, 2] - botz) / dz)

  centroids_energy[ixc, iyc, izc] = centroidvals[j]

  #print ixc, iyc, izc, centroidvals[j]

print "Done "

print "Writing results ..."

mol1 = pybel.readfile("mol2", filenamemol2).next()

write_to_cube (mol1, min_energy, "minim.cube", nx, ny, nz, \
    dx, botx, boty, botz)

write_to_cube (mol1, centroids_energy, "centroid.cube", nx, ny, nz, \
    dx, botx, boty, botz)
"""

print("Done ")
