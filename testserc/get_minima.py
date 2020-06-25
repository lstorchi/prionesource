import re
import os
import sys
import sets
import math
import glob
import numpy
import subprocess

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

import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

def detect_local_minima(arr):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value 
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_minima = local_min - eroded_background
    return numpy.where(detected_minima)       

###############################################################################

def look_near_values (energy, ix, iy, iz, delta, dix, diy, diz):

  nx = energy.shape[0]
  ny = energy.shape[1]
  nz = energy.shape[2]

  minvalue = energy[ix, iy, iz]
  
  minidx = ix - dix
  if minidx < 0:
    minidx = 0
  maxidx = ix + dix
  if maxidx > nx:
    maxidx = nx

  minidy = iy - diy
  if minidy < 0:
    minidy = 0
  maxidy = iy + diy
  if maxidy > ny:
    maxidy = ny

  minidz = iz - diz
  if minidz < 0:
    minidz = 0
  maxidz = iz + diz
  if maxidz > nz:
    maxidz = nz

  listidx = []
  for i in range(minidx, maxidx):
    for j in range(minidy, maxidy):
      for k in range(minidz, maxidz):
        if math.fabs(energy[i, j, k] - minvalue) <= delta:
          listidx.append([i, j, k])

  return listidx

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

print "  field sizes: ", len(xsets), len(ysets), len(zsets)
print "           dx: ", dx
print "           dy: ", dy
print "           dz: ", dz

print "Done..."

print "Start filter e minima location ..."

local_minima_locations = detect_local_minima(energy) 
  
print(energy[local_minima_locations])
print(local_minima_locations)

print " Start writing minima "

for i in range(len(local_minima_locations[0])):
  ix = local_minima_locations[0][i]
  iy = local_minima_locations[1][i]
  iz = local_minima_locations[2][i]
  
  x = botx + ix * dx
  y = boty + iy * dy
  z = botz + iz * dz

  if (energy[ix, iy, iz] <= -3.0):
    print str(x) +  " " + str(y) + " " + str(z) + \
        " " + str(energy[ix, iy, iz])

print " Done "
