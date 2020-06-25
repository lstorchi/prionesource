import scipy.spatial
import subprocess
import numpy 
import pybel
import sets
import math
import sys
import re
import os

MAXSTEPVAL = -0.5
MINSTEPVAL = -5.0

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

    f.close()

    return lines

###############################################################################

def readkontfile (kontname):

  lineamnt = bufcount(kontname)
  
  dim = (lineamnt - 1)/2

  energy = numpy.empty([1,1,1], float)

  if ((dim * 2) + 1) != lineamnt :
    print "Maybe invalid kont file"
    exit(1)

  fk = open(kontname)

  xsets = sets.Set()
  ysets = sets.Set()
  zsets = sets.Set()
  switchtofieldm = False

  nx = ny = nz = 0
  ix = iy = iz = 0
  for l in fk:

    if (l[:5] == "Probe"):
      switchtofieldm = True 
      nx = len(xsets)
      ny = len(ysets)
      nz = len(zsets)
      energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)

    else:
      if switchtofieldm:
        p = re.compile(r'\s+')
        line = p.sub(' ', l)
        line = line.lstrip()
        line = line.rstrip()

        e = float(line)
        energy[ix, iy, iz] = e
        #print ix, iy, iz, e

        # seguo la logica con cui sono scritti i kont ascii senza fare deduzioni
        # ovviamente va migliorato
        iy = iy + 1
        if (iy == ny):
          iy = 0
          ix = ix + 1
        
        if (ix == nx):
          ix = 0
          iy = 0
          iz = iz + 1

        if (iz == nz):
          ix = 0
          iy = 0
          iz = 0

      else:
        p = re.compile(r'\s+')
        line = p.sub(' ', l)
        line = line.lstrip()
        line = line.rstrip()
        n, x, y, z = line.split(" ")

        xsets.add(float(x))
        ysets.add(float(y))
        zsets.add(float(z))

  fk.close()

  dltx = max(xsets) - min(xsets)
  dlty = max(ysets) - min(ysets)
  dltz = max(zsets) - min(zsets)

  return energy, min(xsets), dltx/(nx-1), \
      min(ysets), dlty/(ny-1), min(zsets), dltz/(nz-1)

###############################################################################

filename = ""

if (len(sys.argv)) == 2:
  filename = sys.argv[1]
else:
  print "usage :", sys.argv[0] , " filename.kont"
  exit(1)

energy, xmin, xs, ymin, ys, zmin, zs = readkontfile(filename)

nx = energy.shape[0]
ny = energy.shape[1]
nz = energy.shape[2]

x = xmin
for j in range(0, nx):
  conut = 0

  for i in range(0, nz):
    for k in range(0,ny):
      if energy[j,k,i] <= MAXSTEPVAL:
        conut = conut + 1
  
  print x, " " ,  conut

  x = x + xs
