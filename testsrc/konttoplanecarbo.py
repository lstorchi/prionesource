import re
import os
import sys
import sets
import math
import glob
import numpy
import subprocess


CUBESTEP = 7
MINVAL = -0.1

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

def read_kontfile (kontname):

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

  return energy, botx, boty, botz, topx, topy, topz, dx, dy, dz, nx, ny, nz

###############################################################################

kontname1 = ""
kontname2 = ""

if (len(sys.argv)) == 3: 
  kontname1 = sys.argv[1]
  kontname2 = sys.argv[2]
else:
  print "usage :", sys.argv[0] , " filename1.kont filename2.kont"
  exit(1)


#print "reading \"", kontname1, "\""
energy1, botx1, boty1, botz1, topx1, topy1, topz1, dx1, dy1, dz1, nx1, ny1, nz1 = \
    read_kontfile(kontname1)

#print "reading \"", kontname2, "\""
energy2, botx2, boty2, botz2, topx2, topy2, topz2, dx2, dy2, dz2, nx2, ny2, nz2 = \
    read_kontfile(kontname2)

if ((botx1 != botx2) or (boty1 != boty2) or (botz1 != botz2) or  \
    (topx1 != topx2) or (topy1 != topy2) or (topz1 != topz2) or \
    (dx1 != dx2) or (dy1 != dy2) or (dz1 != dz2) or \
    (nx1 != nx2) or (ny1 != ny2) or (nz1 != nz2)):
  print "non compatible kont files"
  exit(1)

dx = dx1
dy = dy1
dz = dz1

botx = botx1
boty = boty1
botz = botz1
topx = topx1
topy = topy1
topz = topz1

nx = nx1
ny = ny1
nz = nz1

#print "  Box: "
#print "           x: ", botx, " ", topx
#print "           y: ", boty, " ", topy
#print "           z: ", botz, " ", topz
#print "           dx: ", dx
#print "           dy: ", dy
#print "           dz: ", dz

#print "Done..."

#print "Start computing ... "

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

      num = 0.0
      denum1 = 0.0
      denum2 = 0.0
      tot1 = 0
      tot2 = 0
      for i in range(0,CUBESTEP):
        if (ix+i < nx-1):
          for j in range(0,CUBESTEP):
            if (iy+j < ny-1):
              for k in range(0,CUBESTEP):
                if (iz+k < nz-1):
                  #sum1 = 0.0
                  #sum2 = 0.0

                  #if (energy1[ix+i, iy+j, iz+k] <= MINVAL):
                  sum1 = energy1[ix+i, iy+j, iz+k]

                  #if (energy2[ix+i, iy+j, iz+k] <= MINVAL):
                  sum2 = energy2[ix+i, iy+j, iz+k]

                  if (energy1[ix+i, iy+j, iz+k] < MINVAL):
                    tot1 = tot1 + 1

                  if (energy2[ix+i, iy+j, iz+k] < MINVAL):
                    tot2 = tot2 + 1
                 
                  num = num + sum1*sum2
                  denum1 = denum1 + sum1*sum1
                  denum2 = denum2 + sum2*sum2

      if (denum1 == denum2):
        if (denum1 == 0.0):
          carboidx = 1.0
      elif (denum1 == 0.0) or (denum2 == 0.0):
        carboidx = -1.0
      else:
        carboidx = num/math.sqrt(denum1 * denum2)

      xs = xstart + (math.fabs(xstart - xstop)/2.0)
      ys = ystart + (math.fabs(ystart - ystop)/2.0)
      zs = zstart + (math.fabs(zstart - zstop)/2.0)
      #print genc, carboidx, xstart, xstop, ystart, ystop, zstart, zstop
      print genc, tot1 - tot2, xs, ys, zs, xstart, xstop, ystart, ystop, zstart, zstop

      genc = genc + 1

#print "Done "
