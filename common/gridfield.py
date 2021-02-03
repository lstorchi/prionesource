import scipy.spatial
import subprocess
import os.path
import numpy 
import math
import re
import os

import carbo

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

  xsets = set()
  ysets = set()
  zsets = set()
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
        #if e != 0.0:
        #  print (ix, iy, iz, e)

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
        probename, verbose=True):

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

  if verbose:  
    print("Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax)
  
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

    if verbose:  
      print("nx: ", energy1.shape[0], " ny: ", energy1.shape[1], \
         " nz: ", energy1.shape[2])
  
    ifextrm ("./"+kontname)

    if verbose:  
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

def compute_grid_box (filename, delta):

  # generate grid 
  xmin = float("inf")
  ymin = float("inf")
  zmin = float("inf")
  xmax = float("-inf")
  ymax = float("-inf")
  zmax = float("-inf")

  fp = open(filename, "r")

  mollist = []
  for l in fp:

    l = l.replace("\n","")
    
    m = carbo.pdbatomextractor("./"+l)

    mollist.extend(m)
  
  fp.close()

  for conf in mollist:
    for a in conf:
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
  
  xmin = xmin - delta
  xmax = xmax + delta
  
  ymin = ymin - delta
  ymax = ymax + delta
  
  zmin = zmin - delta
  zmax = zmax + delta

  return xmin, xmax, ymin, ymax, zmin, zmax

###############################################################################

def compute_grid_field (filename, box, \
        probename, step, verbose=True, savekont=False):

  xmin = box[0]
  ymin = box[2]
  zmin = box[4]
  xmax = box[1]
  ymax = box[3]
  zmax = box[5]

  if verbose:
    print("Grid will be used: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"%(\
      xmin, ymin, zmin, xmax, ymax, zmax))
  
  toexe = "./fixpdb --remove-all-H2O --unkn-residue-to-grid-types --kout-out="+ \
       filename[:-4]+".kout "+ filename
  results  = subprocess.run(toexe, shell=True, check=True, \
      stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
      universal_newlines=True)
    
  kontname = filename[:-4]+".kont"
  
  fg = open('grid.in','w')
  fg.write("LONT togrid.lont\n")
  fg.write("KONT "+kontname+"\n")
  fg.write("INPT "+filename[:-4]+".kout\n")
  fg.write("NPLA "+str(1.0/step)+"\n")
  fg.write("TOPX "+str(xmax)+"\n")
  fg.write("TOPY "+str(ymax)+"\n")
  fg.write("TOPZ "+str(zmax)+"\n")
  fg.write("BOTX "+str(xmin)+"\n")
  fg.write("BOTY "+str(ymin)+"\n")
  fg.write("BOTZ "+str(zmin)+"\n")
  fg.write(probename+"\n")
  fg.write("IEND\n")
  fg.close()
                                                                                                         
  results  = subprocess.run("./grid grid.in", shell=True, check=True, \
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
    universal_newlines=True)

  ifextrm ("./"+filename[:-4]+".kout")
  ifextrm ("./grid.in")
  ifextrm ("./togrid.lont")
  
  # read kont file
  energy = readkontfile(kontname)

  energy, energy_coords, \
    _ , _ , _ , _ , _ , _ , _ , _ , _ , _ , _ , _  \
       = read_kontfile_and_coords(kontname)
  
  if verbose:
    print("nx: ", energy.shape[0], " ny: ", energy.shape[1], \
      " nz: ", energy.shape[2])
  
  if not savekont:
    ifextrm ("./"+kontname)
  
  if verbose:
    print("Dealing with: ", kontname)
  
  return energy, energy_coords

###############################################################################

def compute_grid_mean_field (filename, step, delta, \
        probename, mol2pdb=True, verbose=True, savekont=False):

  # generate grid 
  xmin = float("inf")
  ymin = float("inf")
  zmin = float("inf")
  xmax = float("-inf")
  ymax = float("-inf")
  zmax = float("-inf")

  fp = open(filename, "r")

  sum = 0.0
  weights = []
  mollist = []
  for l in fp:
    sl = l.split()
    if (len(sl) != 2):
      print("Error in ", filename)
      exit(1)
    
    weights.append(float(sl[1]))
    sum += float(sl[1])

    m = None
    
    if mol2pdb:
      m = carbo.mol2atomextractor(sl[0])
    else:
      m = carbo.pdbatomextractor(sl[0])

    mollist.extend(m)
  
  fp.close()

  weights = [ v/sum for v in weights]

  for conf in mollist:
    for a in conf:
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
  
  xmin = xmin - delta
  xmax = xmax + delta
  
  ymin = ymin - delta
  ymax = ymax + delta
  
  zmin = zmin - delta
  zmax = zmax + delta

  # to set custom ranges
  #xmin = -42.0
  #ymin = -28.0
  #zmin = -42.0

  #xmax = 42.0
  #ymax = 39.0
  #zmax = 42.

  if verbose:
    print("Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax)
  
  if (len(mollist) != len(weights)):
    print("Dimension error ", len(mollist) , " vs " , \
      len(weights))
    exit(1)

  energy = numpy.empty([1,1,1], float)
  globalindex = 0
  fp = open(filename, "r")
  for conf in fp:

    sl = conf.split()
    
    ifextrm ("./"+str(globalindex)+".pdb")

    if mol2pdb:
      toexe = "obabel -imol2 " + sl[0] + " -opdb -O " + "./"+str(globalindex)+".pdb"
      results  = subprocess.run(toexe, shell=True, check=True, \
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
        universal_newlines=True)
    else:
      from shutil import copyfile
      copyfile (sl[0], str(globalindex)+".pdb")

  
    toexe = "./fixpdb --remove-all-H2O --unkn-residue-to-grid-types --kout-out="+ \
        str(globalindex)+".kout "+str(globalindex)+".pdb"
    results  = subprocess.run(toexe, shell=True, check=True, \
      stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
      universal_newlines=True)
    
    kontname = str(globalindex)+".kont"
  
    fg = open('grid.in','w')
    fg.write("LONT togrid.lont\n")
    fg.write("KONT "+kontname+"\n")
    fg.write("INPT "+str(globalindex)+".kout\n")
    fg.write("NPLA "+str(1.0/step)+"\n")
    fg.write("TOPX "+str(xmax)+"\n")
    fg.write("TOPY "+str(ymax)+"\n")
    fg.write("TOPZ "+str(zmax)+"\n")
    fg.write("BOTX "+str(xmin)+"\n")
    fg.write("BOTY "+str(ymin)+"\n")
    fg.write("BOTZ "+str(zmin)+"\n")
    fg.write(probename+"\n")
    fg.write("IEND\n")
    fg.close()
                                                                                                         
    results  = subprocess.run("./grid grid.in", shell=True, check=True, \
      stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
      universal_newlines=True)

    ifextrm ("./"+str(globalindex)+".pdb")
    ifextrm ("./"+str(globalindex)+".kout")
    ifextrm ("./grid.in")
    ifextrm ("./togrid.lont")
  
    # read kont file
    lenergy = readkontfile(kontname)
  
    if verbose:
      print("nx: ", lenergy.shape[0], " ny: ", lenergy.shape[1], \
        " nz: ", lenergy.shape[2])
  
    if savekont:
      newname = sl[0].replace(".pdb", "") + "_" + str(globalindex) + ".kont"
      os.rename("./"+kontname, "./" + newname )
    else:
      ifextrm ("./"+kontname)
  
    if verbose:
      print("Dealing with: ", kontname, " w: ", weights[globalindex])
  
    if  globalindex == 0:
      nx = lenergy.shape[0]
      ny = lenergy.shape[1]
      nz = lenergy.shape[2]
      energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
      energy = numpy.zeros([nx,ny,nz], float)

    energy += weights[globalindex] * lenergy
  
    globalindex = globalindex + 1
  
  #energy = energy / float(globalindex)
  #energytofile (energy, "mean.kont", xmin, ymin, zmin)
  
  fp.close()

  """
  for i in range(energy.shape[0]):
    for j in range(energy.shape[1]):
      for k in range(energy.shape[2]):
        e = energy[i, j, k]
        if e != 0.0:
          print(e)
  """

  return energy, xmin, ymin, zmin

###############################################################################

def read_kontfile_and_coords (kontname):

  lineamnt = bufcount(kontname)
 
  dim = (lineamnt - 1)/2
 
  if ((dim * 2) + 1) != lineamnt :
    print("Maybe invalid kont file")
    exit(1)
 
  fk = open(kontname)
 
  xsets = set()
  ysets = set()
  zsets = set()
  switchtofieldm = False

  energy = numpy.empty([1,1,1], float)
 
  nx = ny = nz = 0
  ix = iy = iz = 0
  for l in fk:

    if "Probe:" in l:
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
        n = ""
        x = ""
        y = ""
        z = ""

        if len(line.split(" ")) < 4:
            n = l[:7]
            x = l[8:16]
            y = l[17:24]
            z = l[25:]
        else:
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

  fk = open(kontname)
 
  coords = numpy.empty([nx,ny,nz], dtype=object)
 
  for iz in range(nz):
      for ix in range(nx):
          for iy in range(ny):
              l = fk.readline()
              p = re.compile(r'\s+')
              line = p.sub(' ', l)
              line = line.lstrip()
              line = line.rstrip()

              n = ""
              x = ""
              y = ""
              z = ""

              if len(line.split(" ")) < 4:
                  n = l[:7]
                  x = l[8:16]
                  y = l[17:24]
                  z = l[25:]
              else:
                  n, x, y, z = line.split(" ")

              coords[ix, iy, iz] = (float(x), float(y), float(z), int(n))

  fk.close()
 
 
  return energy, coords, \
          botx, boty, botz, topx, topy, topz, dx, dy, dz, nx, ny, nz

###############################################################################
