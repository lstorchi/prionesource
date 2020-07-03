import scipy.spatial
import subprocess
import openbabel
import numpy 
import pybel
import sets
import math
import re
import os
import os.path

###############################################################################

def ifextrm (filename):

    if os.path.isfile(filename):
        os.remove(filename)

    return 

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
    print("Maybe invalid kont file")
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

def energytofile (energy, name, botx, boty, botz, STEPVAL):

  ifextrm("./"+name)

  opf = open(name, "w")

  nx = energy.shape[0]
  ny = energy.shape[1]
  nz = energy.shape[2]

  print("Writing final kont... ")
  print("nx: ", nx, " ny: ", ny, " nz: ", nz)

  counter = 1
  for i in range(0, nz):
    z = botz + i * (1.0/STEPVAL)
    for j in range(0, nx):
      x = botx + j * (1.0/STEPVAL)
      for k in range(0, ny):
        y = boty + k * (1.0/STEPVAL)
        opf.write(str(counter) + " " + str(x) +  " " + \
            str(y) + " " + str(z) + "\n")
        counter = counter + 1
  opf.write("Probe: XXX\n")
  for i in range(0, nz):
    for j in range(0, nx):
      for k in range(0, ny):
        opf.write(str(energy[j, k, i]) + "\n")

  opf.close()

###############################################################################

def compute_grid_avg_field (filename1, weightfile, STEPVAL, DELTAVAL, \
        probename):

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
  
  minidx = 0
  maxidx = 1794
  addidx = 0
  
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
  
  # set values:
  #xmin = -43.000
  #xmax =  36.000
  
  #ymin = -32.000
  #ymax =  38.000
  
  #zmin = -29.000
  #zmax =  30.000
  
  print("Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax)
  
  xnsteps = int((xmax - xmin) / 1.0) + 1
  
  weights1 = []
  weightsfp1 = open(weightfile)
  weightsfp1.readline()
  idx = 0
  for line in weightsfp1:
    values = re.split('\s+', line)
    weights1.append(float(values[6]))
    idx = idx + 1
  weightsfp1.close()
  
  if (len(mol1list) != len(weights1)):
    print("Dimension error ", len(mol1list) , " vs " , \
      len(weights1))
    exit(1)
  
  energy = numpy.empty([1,1,1], float)
  globalindex = 0
  for conf1 in mol1list:
    ifextrm ("./"+str(globalindex)+".pdb")

    output = pybel.Outputfile("pdb", str(globalindex)+".pdb")
    output.write(conf1)
    output.close()
  
    toexe = "./fixpdb --remove-all-H2O --unkn-residue-to-grid-types --kout-out="+ \
        str(globalindex)+".kout "+str(globalindex)+".pdb"
    subprocess.call(toexe, shell=True)
    
    kontname = str(globalindex)+".kont"
  
    fg = open('grid.in','w')
    fg.write("LONT togrid.lont\n")
    fg.write("KONT "+kontname+"\n")
    fg.write("INPT "+str(globalindex)+".kout\n")
    fg.write("NPLA "+str(STEPVAL)+"\n")
    fg.write("TOPX "+str(xmax)+"\n")
    fg.write("TOPY "+str(ymax)+"\n")
    fg.write("TOPZ "+str(zmax)+"\n")
    fg.write("BOTX "+str(xmin)+"\n")
    fg.write("BOTY "+str(ymin)+"\n")
    fg.write("BOTZ "+str(zmin)+"\n")
    fg.write(probename+"\n")
    fg.write("IEND\n")
    fg.close()
                                                                                                         
    subprocess.call("./grid grid.in", shell=True)
  
    ifextrm ("./"+str(globalindex)+".pdb")
    ifextrm ("./"+str(globalindex)+".kout")
    ifextrm ("./grid.in")
    ifextrm ("./togrid.lont")
  
    # read kont file
    energy1 = readkontfile(kontname)
  
    print("nx: ", energy1.shape[0], " ny: ", energy1.shape[1], \
        " nz: ", energy1.shape[2])
  
    ifextrm ("./"+kontname)
  
    print("Dealing with: ", kontname, " w: ", weights1[globalindex])
  
    if  globalindex == 0:
      nx = energy1.shape[0]
      ny = energy1.shape[1]
      nz = energy1.shape[2]
      energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
      energy = numpy.zeros([nx,ny,nz], float)
  
    energy = energy + weights1[globalindex] * energy1
  
    globalindex = globalindex + 1
  
  #energy = energy / float(globalindex)
  #energytofile (energy, "mean.kont", xmin, ymin, zmin)
  
  return energy, xmin, ymin, zmin

###############################################################################

def read_kontfile (kontname):

  lineamnt = bufcount(kontname)
 
  dim = (lineamnt - 1)/2
 
  if ((dim * 2) + 1) != lineamnt :
    print("Maybe invalid kont file")
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