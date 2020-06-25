import re
import sys
import math
import mmap
import numpy 

import scipy.spatial

###############################################################################

def mapcount(filename):
  f = open(filename, "r+")
  buf = mmap.mmap(f.fileno(), 0)
  lines = 0
  readline = buf.readline
                      
  while readline():
    lines += 1

  return lines

###############################################################################

STEPVAL = 1.4
DELTAVAL = 20.0
coulombconst = 1.0

filename1 = ""
filename2 = ""

if (len(sys.argv)) == 3:
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]
else:
  print "usage :", sys.argv[0] , " filename1 filename2"
  exit(1)

atomnum1 = mapcount(sys.argv[1])
atomnum2 = mapcount(sys.argv[2])

f1 = open(sys.argv[1])
f2 = open(sys.argv[2])

# read data of first molecule
f1.readline()
# read header
mol1coord = numpy.zeros((atomnum1, 3))
mol1charges = numpy.zeros((atomnum1, 1))
print "Molecule: " , sys.argv[1], atomnum1
idx = 0
for line in f1:
  values = re.split('\s+', line)
  mol1coord[idx,0] = float(values[2])
  mol1coord[idx,1] = float(values[3])
  mol1coord[idx,2] = float(values[4])
  mol1charges[idx,0] = float(values[5])
  idx = idx + 1

# read data of second molecule
f2.readline()
# read header
mol2coord = numpy.zeros((atomnum2, 3))
mol2charges = numpy.zeros((atomnum2, 1))
print "Molecule: " , sys.argv[2], atomnum2
idx = 0
for line in f2:
  values = re.split('\s+', line)
  mol2coord[idx,0] = float(values[2])
  mol2coord[idx,1] = float(values[3])
  mol2coord[idx,2] = float(values[4])
  mol2charges[idx,0] = float(values[5])
  idx = idx + 1

f1.close()
f2.close()

# generate grid 
xmin = min(min(mol2coord[:,0]),min(mol1coord[:,0]))
ymin = min(min(mol2coord[:,1]),min(mol1coord[:,1]))
zmin = min(min(mol2coord[:,2]),min(mol1coord[:,2]))
xmax = max(max(mol2coord[:,0]),max(mol1coord[:,0]))
ymax = max(max(mol2coord[:,1]),max(mol1coord[:,1]))
zmax = max(max(mol2coord[:,2]),max(mol1coord[:,2]))
xmin = xmin - DELTAVAL
ymin = ymin - DELTAVAL
zmin = zmin - DELTAVAL
xmax = xmax + DELTAVAL
ymax = ymax + DELTAVAL
zmax = zmax + DELTAVAL

print "Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax

xnstep = int( ((xmax - xmin) / STEPVAL)+0.5)
ynstep = int( ((ymax - ymin) / STEPVAL)+0.5)
znstep = int( ((zmax - zmin) / STEPVAL)+0.5)


idx = 0
refpoint = numpy.zeros((1, 3))
refpoint[0,0] = xmin - STEPVAL
for ix in range(0,xnstep):
  refpoint[0,0] = refpoint[0,0] + STEPVAL
  refpoint[0,1] = ymin - STEPVAL
  
  num = 0.0
  denum1 = 0.0 
  denum2 = 0.0
  
  for iy in range(0,ynstep):
    refpoint[0,1] = refpoint[0,1] + STEPVAL
    refpoint[0,2] = zmin - STEPVAL
    for iz in range(0,znstep):
      refpoint[0,2] = refpoint[0,2] + STEPVAL
      
      dist1 = scipy.spatial.distance.cdist(mol1coord,refpoint)
      dist2 = scipy.spatial.distance.cdist(mol2coord,refpoint)
  
      sum1 = 0.0
      sum2 = 0.0
  
      if (dist1.min() > 1.0):
        ep1 = coulombconst * (mol1charges/dist1)
        sum1 = numpy.sum(ep1) 
  
      if (dist2.min() > 1.0):
        ep2 = coulombconst * (mol2charges/dist2)
        sum2 = numpy.sum(ep2)

      num = num + sum1*sum2
      denum1 = denum1 + sum1*sum1 
      denum2 = denum2 + sum2*sum2 
  
  carboidx = num/math.sqrt(denum1 * denum2)
  
  print refpoint[0,0], carboidx 

  idx = idx + 1
