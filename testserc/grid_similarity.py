import re
import os
import sys
import sets
import math
import glob
import numpy
import subprocess

sys.path.append("./common")
import gridfield

###############################################################################

def remove_equal(input):
  output = []
  for x in input:
    if x not in output:
      output.append(x)
                        
  return output

###############################################################################

CUBESTEP = 7
MINVAL = -0.1
STEPVAL = 1.4
DELTAVAL = 10.0

probename = "DRY"

filename1 = ""
weightfile1 = ""
filename2 = ""
weightfile2 = ""

if (len(sys.argv)) == 5:
  filename1 = sys.argv[1]
  weightfile1 = sys.argv[2]
  filename2 = sys.argv[3]
  weightfile2 = sys.argv[4]
else:
  print "usage :", sys.argv[0] , " filename1.mol2 weight1 filename2.mol2 weight2"
  exit(1)

kontname1 = os.path.splitext(filename1)[0]
kontname2 = os.path.splitext(filename2)[0]

kontname1 = kontname1+".kont"
kontname2 = kontname2+".kont"

energy1, xmin1, ymin1, zmin1 = \
        gridfield.compute_grid_avg_field (filename1, weightfile1, \
        STEPVAL, DELTAVAL, probename)

gridfield.energytofile (energy1, kontname1, xmin1, ymin1, zmin1, STEPVAL)

energy2, xmin2, ymin2, zmin2 = \
        gridfield.compute_grid_avg_field (filename2, weightfile2, \
        STEPVAL, DELTAVAL, probename)

gridfield.energytofile (energy2, kontname2, xmin2, ymin2, zmin2, STEPVAL)

energy1, botx1, boty1, botz1, topx1, topy1, topz1, dx1, dy1, dz1, nx1, ny1, nz1 = \
    gridfield.read_kontfile(kontname1)

energy2, botx2, boty2, botz2, topx2, topy2, topz2, dx2, dy2, dz2, nx2, ny2, nz2 = \
    gridfield.read_kontfile(kontname2)

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
      print genc, tot1 - tot2, xs, ys, zs, xstart, xstop, ystart, ystop, zstart, zstop

      genc = genc + 1
