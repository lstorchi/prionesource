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
weifile1 = ""
filename2 = ""
weifile2 = ""

if (len(sys.argv)) == 5:
  filename1 = sys.argv[1]
  weifile1 = sys.argv[2]
  filename2 = sys.argv[3]
  weifile2 = sys.argv[4]
else:
  print "usage :", sys.argv[0] , " listfilename1 weifile listfilename2 weifile"
  exit(1)

# weights file, this first nomalization is maybe useless
weights1 = []
weights2 = []
weightsfp1 = open(weifile1)
weightsfp2 = open(weifile2)
weightsfp1.readline()
idx = 0
sum = 0.0
for line in weightsfp1:
  values = re.split('\s+', line)
  sum = sum + float(values[6])
  weights1.append(float(values[6]))
  idx = idx + 1
for i in range(len(weights1)):
  weights1[i] = weights1[i] / sum
sum = 0.0
weightsfp2.readline()
idx = 0
for line in weightsfp2:
  values = re.split('\s+', line)
  sum = sum + float(values[6])
  weights2.append(float(values[6]))
  idx = idx + 1
for i in range(len(weights2)):
  weights2[i] = weights2[i] / sum
weightsfp1.close()
weightsfp2.close()

mol1num = mapcount(filename1)
mol2num = mapcount(filename2)

fp1 = open(filename1)
f1 = open(fp1.readline().rstrip())
atomnum1 = mapcount(fp1.readline().rstrip()) + 1
mol1coord = numpy.zeros((atomnum1, 3))
mol1charges = numpy.zeros((atomnum1, 1))
idx = 0
f1.readline()
for line in f1:
  values = re.split('\s+', line)
  mol1coord[idx,0] = float(values[2])
  mol1coord[idx,1] = float(values[3])
  mol1coord[idx,2] = float(values[4])
  mol1charges[idx,0] = float(values[5])
  idx = idx + 1

fp2 = open(filename2)
f2 = open(fp2.readline().rstrip())
atomnum2 = mapcount(fp2.readline().rstrip()) + 1
mol2coord = numpy.zeros((atomnum2, 3))
mol2charges = numpy.zeros((atomnum2, 1))
idx = 0
f2.readline()
for line in f2:
  values = re.split('\s+', line)
  mol2coord[idx,0] = float(values[2])
  mol2coord[idx,1] = float(values[3])
  mol2coord[idx,2] = float(values[4])
  mol2charges[idx,0] = float(values[5])
  idx = idx + 1

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
 
f1.close()
f2.close()
fp1.close()
fp2.close()

xrefpoints = numpy.zeros(xnstep)
carboidxs = numpy.zeros((mol1num*mol2num, xnstep))
weights = numpy.zeros(mol1num*mol2num)

# compute weight
generalidx = 0
for mol1idx in range(0,len(weights1)):
  for mol2idx in range(0,len(weights2)):
    # compute weights
    weights[generalidx] = weights1[mol1idx] + weights2[mol2idx]
    weights[generalidx] = weights[generalidx] / 2.0
                      
    generalidx = generalidx + 1

globalidx = 0
fp1 = open(filename1)
fp2 = open(filename2)
allfile1 = []
for fname1 in fp1:
  allfile1.append(fname1)
fp1.close()
allfile2 = []
for fname2 in fp2:
  allfile2.append(fname2)
fp2.close()

for fname1 in allfile1:
  for fname2 in allfile2:

    atomnum1 = mapcount(fname1.rstrip()) + 2
    atomnum2 = mapcount(fname2.rstrip()) + 2

    f1 = open(fname1.rstrip())
    f2 = open(fname2.rstrip())

    # read data of first molecule
    f1.readline()
    # read header
    mol1coord = numpy.zeros((atomnum1, 3))
    mol1charges = numpy.zeros((atomnum1, 1))
    print "Molecule: " , fname1.rstrip(), atomnum1
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
    print "Molecule: " , fname2.rstrip(), atomnum2
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
      xrefpoints[idx] = refpoint[0,0]
      carboidxs[globalidx,idx] = carboidx
    
      idx = idx + 1

    globalidx = globalidx + 1
    print "Done: ", globalidx, " of " , mol1num*mol2num

print "Weighted Mead and stdev"
waverage = numpy.average(carboidxs, 0, weights)
wvariance = numpy.average((carboidxs-waverage)**2, 0, weights)
idx = 0
for std in wvariance:
  print xrefpoints[idx], " ", waverage[idx] , " ", std
  idx = idx + 1

print "full matrix"
print waverage
print wvariance
