terrom scipy.cluster.vq import *

import Pycluster 

import numpy 
import scipy.spatial
import pybel
import sys
import math

STEPVAL = 5.0
DELTAVAL = 10.0
coulombconst = 1.0

filename1 = ""

if (len(sys.argv)) == 2:
  filename1 = sys.argv[1]
else:
  print "usage :", sys.argv[0] , " filename1.mol2"
  exit(1)

# read first molecule in the files
mol1 = pybel.readfile("mol2", filename1).next()

# read data of first molecule
atomnum1 = len(mol1.atoms)
mol1coord = numpy.zeros((atomnum1, 3))
mol1charges = numpy.zeros((atomnum1, 1))
print "Molecule: " , mol1.title
idx = 0
for atom in mol1:
  mol1coord[idx,0] = atom.coords[0]
  mol1coord[idx,1] = atom.coords[1]
  mol1coord[idx,2] = atom.coords[2]
  mol1charges[idx,0] = atom.partialcharge
  idx = idx + 1

# generate grid 
xmin = min(mol1coord[:,0])
ymin = min(mol1coord[:,1])
zmin = min(mol1coord[:,2])
xmax = max(mol1coord[:,0])
ymax = max(mol1coord[:,1])
zmax = max(mol1coord[:,2])
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

# read all molecules 
mol1list = list(pybel.readfile("mol2", filename1))

mol1field = numpy.zeros((len(mol1list), xnstep*ynstep*znstep))

for molidx in range(0,len(mol1list)):
  # read coord first molecule
  idx = 0
  for atom in mol1list[molidx]:
    mol1coord[idx,0] = atom.coords[0]
    mol1coord[idx,1] = atom.coords[1]
    mol1coord[idx,2] = atom.coords[2]
    idx = idx + 1

  refpoint = numpy.zeros((1, 3))
  refpoint[0,0] = xmin - STEPVAL
  for ix in range(0,xnstep):
    refpoint[0,0] = refpoint[0,0] + STEPVAL
    refpoint[0,1] = ymin - STEPVAL
    
    
    for iy in range(0,ynstep):
      refpoint[0,1] = refpoint[0,1] + STEPVAL
      refpoint[0,2] = zmin - STEPVAL
      for iz in range(0,znstep):
        refpoint[0,2] = refpoint[0,2] + STEPVAL
        
        #print refpoint[0,0], refpoint[0,1], refpoint[0,2], 1.0
        dist1 = scipy.spatial.distance.cdist(mol1coord,refpoint)
    
        sum1 = 0.0
    
        if (dist1.min() > 1.0):
          # per essere certo che non sono interno posso valuate se ho 
          # atomi a destra sinistra basso alto ... 
          ep1 = coulombconst * (mol1charges/dist1)
          sum1 = numpy.sum(ep1) 
        #per valurae ad esempio se non prendo i punti interni 
        #else:
        #  print refpoint[0,0], refpoint[0,1], refpoint[0,2], sum1
    
        counter = ix*ynstep*znstep + iy*znstep + iz
        mol1field[molidx, counter] = mol1field[molidx, counter] + sum1

  print "Done ", molidx+1, " of ", len(mol1list)

#print mol1field

#print "Start whitening" 
#whiten(mol1field)

print "Start k-means on the set"
labels, error, nfound = Pycluster.kcluster(mol1field,5)
#res, idx = kmeans2(mol1field,2)

print labels, error, nfound
