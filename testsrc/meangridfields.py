import scipy.spatial
import subprocess
import openbabel
import numpy 
import pybel
import sets
import math
import sys
import re
import os

CONVERTER = 1.889716
STEPVAL = 1.0
DELTAVAL = 10.0

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

    if probename in l:
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

def write_to_cube (mol1, mol1field, fname, xnstep, ynstep, znstep,\
    step, xmin, ymin, zmin):

  opf = open(fname, "w")

  print "Writing final cube... "

  zero = 0.0
  opf.write("El field\n")
  opf.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (len(mol1.atoms), CONVERTER*xmin, \
      CONVERTER*ymin, CONVERTER*zmin))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (xnstep, CONVERTER*step, zero, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (ynstep, zero, CONVERTER*step, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (znstep, zero, zero, CONVERTER*step))

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


###############################################################################


filename1 = ""
weightfile = ""
probename = ""

if (len(sys.argv)) == 4:
  filename1 = sys.argv[1]
  weightfile = sys.argv[2]
  probename = sys.argv[3]
else:
  print "usage :", sys.argv[0] , " filename.mol2 weight probename"
  exit(1)

# generate grid 
xmin = float("inf")
ymin = float("inf")
zmin = float("inf")
xmax = float("-inf")
ymax = float("-inf")
zmax = float("-inf")

readlist = list(pybel.readfile("mol2", filename1))


mol1list = readlist

"""

#minidx = 0
#maxidx = 1794
#addidx = 0

minidx = 1794
maxidx = 4000
addidx = 1

mol1list = []

molidx = 1
for mol in readlist:
    obmol = openbabel.OBMol()

    for i in range(len(mol.atoms)):
        obatom = mol.atoms[i].OBAtom
        idx = obatom.GetIdx()
        if idx < maxidx and idx >= minidx:
            newidx = idx - minidx
            newobatom = openbabel.OBAtom()
            newobatom.SetAtomicNum(obatom.GetAtomicNum())
            #newobatom.SetVector(obatom.GetVector())
            newobatom.SetFormalCharge(obatom.GetFormalCharge())
            newobatom.SetFormalCharge(obatom.GetFormalCharge())
            newobatom.SetPartialCharge(obatom.GetPartialCharge())
            if (obatom.IsChiral()):
                newobatom.SetChiral()
            if (obatom.IsAromatic()):
                newobatom.SetAromatic ()
            if (obatom.IsClockwise()):
                newobatom.SetClockwise ()
            if (obatom.IsAntiClockwise()):
                newobatom.SetAntiClockwise ()
            if (obatom.IsPositiveStereo ()):
                newobatom.SetPositiveStereo ()
            if (obatom.IsNegativeStereo ()):
                newobatom.SetNegativeStereo ()
            newobatom.SetIsotope(obatom.GetIsotope())
            newobatom.SetType(obatom.GetType())
            newobatom.SetVector(obatom.x(), obatom.y(), obatom.z())
            newobatom.SetResidue(obatom.GetResidue())
            #newobatom.SetIdx(newidx)
            newobatom.SetId(newidx)
            obmol.AddAtom(newobatom)

    #for atom in openbabel.OBMolAtomIter(obmol):
    #    print atom.GetIdx(), atom.GetId(), atom.GetAtomicNum() 

    bondidx = 0
    for bond in openbabel.OBMolBondIter(mol.OBMol):
        startidx = bond.GetBeginAtom().GetIdx()
        endidx = bond.GetEndAtom().GetIdx()
        bondord = bond.GetBO()

        if startidx < maxidx and startidx >= minidx and \
                endidx < maxidx and endidx >= minidx :
            #print "renum: ", startidx-minidx, endidx-minidx, bondord
            obmol.AddBond(startidx-minidx+addidx, endidx-minidx+addidx, bondord, 0, \
                    bondidx)
            bondidx = bondidx + 1

    pybelmol = pybel.Molecule(obmol)

    pybelmol.write("mol2", str(molidx)+".mol2", True)
    mol1list.append(pybelmol)

    molidx = molidx + 1


"""

for conf1 in mol1list:
  for a in conf1.atoms:
    x, y, z = a.coords

    #x = x * CONVERTER
    #y = y * CONVERTER
    #z = z * CONVERTER
    
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

weights1 = []
weightsfp1 = open(weightfile)
weightsfp1.readline()
idx = 0
sum = 0.0
for line in weightsfp1:
  values = re.split('\s+', line)
  weights1.append(float(values[6]))

  sum = sum + float(values[6])

  idx = idx + 1

for i in range(len(weights1)):
  weights1[i] = weights1[i] / sum

weightsfp1.close()


if (len(mol1list) != len(weights1)):
  print "Dimension error ", len(mol1list) , " vs " , \
    len(weights1)
  exit(1)

energy = numpy.empty([1,1,1], float)
globalindex = 0
for conf1 in mol1list:
  outfl = filename1
  outfl = outfl.replace(".mol2", "_")
  basicfname = outfl+str(globalindex)

  output = pybel.Outputfile("pdb", basicfname+".pdb")
  output.write(conf1)
  output.close()

  toexe = "./fixpdb --remove-all-H2O --unkn-residue-to-grid-types --kout-out="+ \
      basicfname+".kout "+basicfname+".pdb"
  subprocess.call(toexe, shell=True)

  kontname = basicfname+".kont"

  limits = str(xmin)+";"+str(xmax)+";"+str(ymin)+";"+ \
          str(ymax)+";"+str(zmin)+";"+str(zmax)

  torun = "./kouttokont -g \""+limits+"\" "+basicfname+".kout "+probename + \
          " > "+kontname

  subprocess.call(torun, shell=True)


  os.remove(basicfname+".pdb")
  os.remove(basicfname+".kout")

  # read kont file
  energy1 = readkontfile(kontname)

  print "nx: ", energy1.shape[0], " ny: ", energy1.shape[1], \
      " nz: ", energy1.shape[2]

  os.remove(kontname)

  print "Dealing with: ", kontname, " w: ", weights1[globalindex]

  if  globalindex == 0:
    nx = energy1.shape[0]
    ny = energy1.shape[1]
    nz = energy1.shape[2]
    energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
    energy = numpy.zeros([nx,ny,nz], float)

  energy = energy + weights1[globalindex] * energy1

  globalindex = globalindex + 1

#energy = energy / float(globalindex)

mol1 = pybel.readfile("mol2", filename1).next()

nx = energy.shape[0]
ny = energy.shape[1]
nz = energy.shape[2]

outfilename = filename1
outfilename = outfilename.replace(".mol2", "_"+probename+".cube")

write_to_cube (mol1, energy, outfilename, nx, ny, nz,\
    1.0/STEPVAL, xmin, ymin, zmin)
