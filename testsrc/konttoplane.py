import re
import os
import sys
import sets
import math
import glob
import numpy
import subprocess


CUBESTEP = 5
MINVAL = -2.0

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

print "Read file... "

dx = 0.0
dy = 0.0
dz = 0.0

botx = 0.0
boty = 0.0
botz = 0.0
topx = 0.0
topy = 0.0
topz = 0.0

kontname = ""

if (len(sys.argv)) == 2: 
  kontname = sys.argv[1]
else:
  print "usage :", sys.argv[0] , " filename.kont"
  exit(1)

print "reading \"", kontname, "\""

lineamnt = bufcount(kontname)

dim = (lineamnt - 1)/2

if ((dim * 2) + 1) != lineamnt :
  print "Maybe invalid kont file"
  exit(1)

fk = open(kontname)

xsets = sets.Set()
ysets = sets.Set()
zsets = sets.Set()
switchtofieldm = False

energy = numpy.empty([1,1,1], float)
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

dx = sorted(xsets)[1] - sorted(xsets)[0]
dy = sorted(ysets)[1] - sorted(ysets)[0]
dz = sorted(zsets)[1] - sorted(zsets)[0]

botx = min(list(xsets))
boty = min(list(ysets))
botz = min(list(zsets))

topx = max(list(xsets))
topy = max(list(ysets))
topz = max(list(zsets))

print "  Box "
print "           x: ", botx, " ", topx
print "           y: ", boty, " ", topy
print "           z: ", botz, " ", topz

print "  field sizes: ", len(xsets), len(ysets), len(zsets)
print "           dx: ", dx
print "           dy: ", dy
print "           dz: ", dz

print "Done..."

print "Start computing ... "

genc = 0
for ix in range(0,nx,CUBESTEP):
  for iy in range(0,ny,CUBESTEP):
    for iz in range(0,nz,CUBESTEP):
      xstart = botx + ix * dx
      ystart = boty + iy * dy
      zstart = botz + iz * dz

      xstop = botx + (ix+5) * dx
      ystop = boty + (iy+5) * dy
      zstop = botz + (iz+5) * dz

      #print xstart, xstop, ystart, ystop, zstart, zstop
      #print "start : ", iz

      couter = 0
      for i in range(0,CUBESTEP):
        if (ix+i < nx-1):
          for j in range(0,CUBESTEP):
            if (iy+j < ny-1):
              for k in range(0,CUBESTEP):
                if (iz+k < nz-1):
                  if (energy[ix+i, iy+j, iz+k] <= MINVAL):
                    couter = couter + 1
      
      print genc, couter, xstart, xstop, ystart, ystop, zstart, zstop

      genc = genc + 1

print "Done "
