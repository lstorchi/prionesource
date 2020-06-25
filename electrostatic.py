import scipy.spatial
import numpy 
import pybel
import math
import sys
import re

CONVERTER = 1.889716
STEPVAL = 1.0*CONVERTER
DELTAVAL = 10.0*CONVERTER
coulombconst = 1.0

#####################################################################

def write_to_kont (mol1field, fname, xnstep, ynstep, znstep, step):

  opf = open(fname, "w")

  print "Writing final kont... "

  counter = 1
  for i in range(0, znstep):
    z = zmin + i * (step)
    for j in range(0, xnstep):
      x = xmin + j * (step)
      for k in range(0, ynstep):
        y = ymin + k * (step)
        opf.write(str(counter) + " " + str(x) +  " " + \
            str(y) + " " + str(z) + "\n")
        counter = counter + 1
  opf.write("Probe: XXX\n")
  for i in range(0, znstep):
    for j in range(0, xnstep):
      for k in range(0, ynstep):
        opf.write(str(mol1field[j, k, i]) + "\n")

  opf.close()

#####################################################################

def write_to_cube (mol1, mol1field, fname, xnstep, ynstep, znstep,\
    step, xmin, ymin, zmin):

  opf = open(fname, "w")

  print "Writing final cube... "

  zero = 0.0
  opf.write("El field\n")
  opf.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (len(mol1.atoms), xmin, ymin, zmin))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (xnstep, step, zero, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (ynstep, zero, step, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (znstep, zero, zero, step))

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

#####################################################################


filename1 = ""
weightsname1 = ""

if (len(sys.argv)) == 3:
  filename1 = sys.argv[1]
  weightsname1 = sys.argv[2]
else:
  print "usage :", sys.argv[0] , " filename1.mol2 weights1.txt"
  exit(1)

# read weights file, this first nomalization is maybe useless
weights1 = []
weightsfp1 = open(weightsname1)
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

weightsfp1.close()

# read first molecule in the files
mol1 = pybel.readfile("mol2", filename1).next()

# read data of first molecule
atomnum1 = len(mol1.atoms)
mol1coord = numpy.zeros((atomnum1, 3))
mol1charges = numpy.zeros((atomnum1, 1))
print "Molecule: " , mol1.title
idx = 0
for atom in mol1:
  mol1coord[idx,0] = CONVERTER*atom.coords[0]
  mol1coord[idx,1] = CONVERTER*atom.coords[1]
  mol1coord[idx,2] = CONVERTER*atom.coords[2]
  mol1charges[idx,0] = atom.partialcharge
  idx = idx + 1

# generate grid 
xmin = min(min(mol1coord[:,0]),min(mol1coord[:,0]))
ymin = min(min(mol1coord[:,1]),min(mol1coord[:,1]))
zmin = min(min(mol1coord[:,2]),min(mol1coord[:,2]))
xmax = max(max(mol1coord[:,0]),max(mol1coord[:,0]))
ymax = max(max(mol1coord[:,1]),max(mol1coord[:,1]))
zmax = max(max(mol1coord[:,2]),max(mol1coord[:,2]))
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

# check amount 
if (len(mol1list) != len(weights1)):
  print "Error different size"
  exit()

mol1field = numpy.zeros((xnstep, ynstep, znstep))
probecharge = -1.0

for mol1idx in range(0,len(weights1)):
  print "Step ", mol1idx+1, " of ", len(weights1)
  # read coord first molecule
  idx = 0
  atomnum1 = len(mol1list[mol1idx].atoms)
  mol1coord = numpy.zeros((atomnum1, 3))
  mol1charges = numpy.zeros((atomnum1, 1))
  for atom in mol1list[mol1idx]:
    mol1coord[idx,0] = CONVERTER*atom.coords[0]
    mol1coord[idx,1] = CONVERTER*atom.coords[1]
    mol1coord[idx,2] = CONVERTER*atom.coords[2]
    mol1charges[idx,0] = atom.partialcharge
    idx = idx + 1
    
  idx = 0
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
          ep1 = coulombconst * ((probecharge*mol1charges)/dist1)
          sum1 = numpy.sum(ep1)
        #per valurtare ad esempio se non prendo i punti interni 
        #else:
        #  print refpoint[0,0], refpoint[0,1], refpoint[0,2], sum1
    
        mol1field[ix,iy,iz] = mol1field[ix,iy,iz] + sum1*weights1[mol1idx]
    

#write_to_kont (mol1field, "final.kont", xnstep, ynstep, znstep, STEPVAL)

outfilename = filename1
outfilename = outfilename.replace(".mol2", "_elk.cube")

print "write result in ", outfilename
write_to_cube (mol1, mol1field, outfilename, xnstep, ynstep, znstep,\
    STEPVAL, xmin, ymin, zmin)
