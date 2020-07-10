import scipy.spatial
import subprocess
import numpy 
import pybel
import sets
import math
import sys
import re
import os

STEPVAL = 1.4
DELTAVAL = 20.0

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

  return energy

###############################################################################

filename1 = ""
filename2 = ""

if (len(sys.argv)) == 3:
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]
else:
  print "usage :", sys.argv[0] , " filename1.mol2 filename2.mol2"
  exit(1)

# generate grid 
xmin = float("inf")
ymin = float("inf")
zmin = float("inf")
xmax = float("-inf")
ymax = float("-inf")
zmax = float("-inf")

mol1list = list(pybel.readfile("mol2", filename1))
mol2list = list(pybel.readfile("mol2", filename2))

for conf1 in mol1list:
  for a in conf1.atoms:
    x, y, z = a.coords
    
    if x < xmin:
      xmin = x
    if y < ymin:
      ymin = y
    if z < zmin:
      zmin = z

    if x > xmax:
      xmax = x
    if y > ymax:
      ymax = y
    if z > zmax:
      zmax = z

for conf2 in mol2list:
  for a in conf2.atoms:
    x, y, z = a.coords
    
    if x < xmin:
      xmin = x
    if y < ymin:
      ymin = y
    if z < zmin:
      zmin = z

    if x > xmax:
      xmax = x
    if y > ymax:
      ymax = y
    if z > zmax:
      zmax = z

xmin = xmin - DELTAVAL
xmax = xmax + DELTAVAL

ymin = ymin - DELTAVAL
ymax = ymax + DELTAVAL

zmin = zmin - DELTAVAL
zmax = zmax + DELTAVAL

print "Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax

xnsteps = int((xmax - xmin) / 1.0) + 1
carboidxs = numpy.zeros((len(mol1list)*len(mol2list), xnsteps))

globalindex = 1
for conf1 in mol1list:
  output = pybel.Outputfile("pdb", str(globalindex)+".pdb")
  output.write(conf1)
  output.close()

  toexe = "./fixpdb --unkn-residue-to-grid-types --kout-out="+ \
      str(globalindex)+".kout "+str(globalindex)+".pdb"
  subprocess.call(toexe, shell=True)
  
  kontname = str(globalindex)+".kont"

  fg = open('grid.in','w')
  fg.write("LONT togrid.lont\n")
  fg.write("KONT "+kontname+"\n")
  fg.write("INPT "+str(globalindex)+".kout\n")
  fg.write("NPLA 1.0\n")
  fg.write("TOPX "+str(xmax)+"\n")
  fg.write("TOPY "+str(ymax)+"\n")
  fg.write("TOPZ "+str(zmax)+"\n")
  fg.write("BOTX "+str(xmin)+"\n")
  fg.write("BOTY "+str(ymin)+"\n")
  fg.write("BOTZ "+str(zmin)+"\n")
  fg.write("OH2\n")
  fg.write("IEND\n")
  fg.close()
                                                                                                       
  subprocess.call("./grid grid.in", shell=True)

  os.remove(str(globalindex)+".pdb")
  os.remove(str(globalindex)+".kout")
  os.remove("grid.in")
  os.remove("togrid.lont")

  # read kont file
  energy1 = readkontfile(kontname)

  os.remove(kontname)

  #print "Dealing with: ", kontname

  for conf2 in mol2list:

    globalindex = globalindex + 1

    output = pybel.Outputfile("pdb", str(globalindex)+".pdb")
    output.write(conf2)
    output.close()
    
    toexe = "./fixpdb --unkn-residue-to-grid-types --kout-out="+ \
        str(globalindex)+".kout "+str(globalindex)+".pdb"
    subprocess.call(toexe, shell=True)
    
    kontname_o = str(globalindex)+".kont"
    
    fg = open('grid.in','w')
    fg.write("LONT togrid.lont\n")
    fg.write("KONT "+kontname_o+"\n")
    fg.write("INPT "+str(globalindex)+".kout\n")
    fg.write("NPLA 1.0\n")
    fg.write("TOPX "+str(xmax)+"\n")
    fg.write("TOPY "+str(ymax)+"\n")
    fg.write("TOPZ "+str(zmax)+"\n")
    fg.write("BOTX "+str(xmin)+"\n")
    fg.write("BOTY "+str(ymin)+"\n")
    fg.write("BOTZ "+str(zmin)+"\n")
    fg.write("OH2\n")
    fg.write("IEND\n")
    fg.close()
                                                                                                         
    subprocess.call("./grid grid.in", shell=True)
    
    os.remove(str(globalindex)+".pdb")
    os.remove(str(globalindex)+".kout")
    os.remove("grid.in")
    os.remove("togrid.lont")
    
    # read kont file
    energy2 = readkontfile(kontname_o)
    
    os.remove(kontname_o)
    
    #print "Dealing with: ", kontname

    print "Starting to compare ", kontname , " and ", kontname_o

    # carboindexcoparison

    nx = energy1.shape[0]
    ny = energy1.shape[1]
    nz = energy1.shape[2]

    if (energy2.shape[0] != nx):
      print "Error in shape"
      exit(1)

    if (energy2.shape[1] != ny):
      print "Error in shape"
      exit(1)

    if (energy2.shape[2] != nz):
      print "Error in shape"
      exit(1)

    for i in range(0, nx):
      num = 0.0
      denum1 = 0.0
      denum2 = 0.0

      for j in range(0, ny):
        for k in range(0, nz):

          sum1 = energy1[i,j,k]
          sum2 = energy2[i,j,k]

          num = num + sum1*sum2
          denum1 = denum1 + sum1*sum1 
          denum2 = denum2 + sum2*sum2 
      
      carboidx = num/math.sqrt(denum1 * denum2)

      if (denum1 == denum2):
        if (denum1 == 0.0):
          carboidx = 1.0

      if (denum1 == 0.0) or (denum2 == 0.0):
        carboidx = -1.0

      carboidxs[globalindex-2,i] = carboidx

      print i, " " , carboidx


print "Weighted Mead and stdev"
average = carboidxs.mean(0)
variance = carboidxs.std(0)
idx = 0
xref = xmin
for std in variance:
  print xref, " ", average[idx] , " ", std
  idx = idx + 1
  xref = xref + 1.0
