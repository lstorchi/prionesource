import scipy.spatial
import numpy 
import math
import re
import os

from scipy.io import FortranFile

IDXSLCT = 5

###############################################################

def get_distance_from_line (lp1, lp2, p0):
  
  # equation of ile sta + t * dir
  # dir = lp2 - lp1
  # sta = lp1

  # lp1 and lp2 define the line
  den = numpy.sqrt(numpy.sum(numpy.square(lp2 - lp1))) 
  cp = numpy.cross (p0 - lp1, p0 - lp2)
  num = numpy.sqrt(numpy.sum(numpy.square(cp))) 

  return num / den

###############################################################

class atom:

  def __init__ (self, id, name, x, y, z, c):
    self.id = id
    self.name = name
    self.coords = (x, y, z)
    self.partialcharge = c
    self.resname = ""
    self.atomname = ""
  
  def __repr__ (self):
    line = "%6d %6s %6s %10.5f %10.5f %10.5f %8.5f %3s\n"%( \
      self.id, self.name, self.resname, self.coords[0], \
      self.coords[1], self.coords[2], self.partialcharge, \
      self.atomname)

    return line
  
###############################################################

def get_eps (refpoint, mol1coord, mol1charges, \
  mol2coord, mol2charges, coulombconst):
            
  #print refpoint[0,0], refpoint[0,1], refpoint[0,2], 1.0
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

  return sum1, sum2

###############################################################

def read_delphi_file (filename):

  f = FortranFile(filename, 'r' )

  #character*20 uplbl
  #character*10 nxtlbl,character*60 toplbl
  #real*4 phi(65,65,65)
  #character*16 botlbl
  #real*4 scale,oldmid(3)

  uplbl = f.read_record('a20')
  nxtolbl = f.read_record('a70')
  epmap = f.read_reals(dtype='float32').reshape((65,65,65))
  botlbl = f.read_record('a16')
  scalemin = f.read_reals(dtype='float32')

  scale = scalemin[0]
  oldmid = scalemin[1:4]

  idx = 1
  for iz in range(65):
    IZ = iz + 1
    for ix in range(65):
      IX = ix + 1
      for iy in range(65):
        IY = iy + 1

        x = (IX - 33)/scale + oldmid[0]
        y = (IY - 33)/scale + oldmid[1]
        z = (IZ - 33)/scale + oldmid[2]

        #print("%6d %8.3f %8.3f %8.3f %10.5f"%(idx, x, y, z, \
        #  epmap[ix, iy, iz]))
        print("%6d %8.3f %8.3f %8.3f"%(idx, x, y, z))

        idx += 1

  print("EP")
  for iz in range(65):
    for ix in range(65):
      for iy in range(65):
        print("%8.3f"%(epmap[ix, iy, iz]))

  f.close()

###############################################################

def get_damped_eps (refpoint, mol1coord, mol1charges, \
  mol2coord, mol2charges, coulombconst):

  raise Exception("To be implemented")

  #print refpoint[0,0], refpoint[0,1], refpoint[0,2], 1.0
  dist1 = scipy.spatial.distance.cdist(mol1coord, refpoint)
  dist2 = scipy.spatial.distance.cdist(mol2coord, refpoint)
  
  avgatomcr = 3.0
  dampsfactors = []
  # compute dumping values
  for p2 in mol1coord:
    dists = [get_distance_from_line(refpoint, p2, p0) for p0 in mol1coord]
    df = sum(d < avgatomcr for d in dists)
    dampsfactors.append(df)
    print(df)

  # can compute once at the beginning need and stored per each point in the grid
  # compute end graficalli check 

  #print(dampsfactors)

  #d = get_distance_from_line (numpy.asarray( [4, 2, 4]) \
  #   ,  numpy.asarray([5, 5, 9]), numpy.asarray([1, 1, 2]))
  # chec 2.59 ...

  sum1 = 0.0
  sum2 = 0.0
        
  if (dist1.min() > 1.0):
    ep1 = coulombconst * (mol1charges/dist1)
    sum1 = numpy.sum(ep1) 
       
  if (dist2.min() > 1.0):
    ep2 = coulombconst * (mol2charges/dist2)
    sum2 = numpy.sum(ep2)

  return sum1, sum2

###############################################################

def mol2atomextractor (file=None, readresname= False):
  
  mols = []

  with open(file, 'r') as f:
      while not f.tell() == os.fstat(f.fileno()).st_size:
        line = f.readline()
        if line.startswith("@<TRIPOS>ATOM"):
            line = f.readline()
            mol = []
            while not line.startswith("@<TRIPOS>BOND"):
              sline = line.split()
              if len(sline) != 9 and len(sline) != 10 :
                raise Exception("Error in "+line+ " line is "+str(len(sline)))
              a = atom(int(sline[0]), sline[1], float(sline[2]), \
                  float(sline[3]),float(sline[4]), float(sline[8]) )

              if readresname:
                a.resname = sline[7][0:3]
                a.atomname = sline[5].split(".", 1)[0]

              line = f.readline()
              mol.append(a)
            mols.append(mol)

  return(mols)

###############################################################

def carbo_similarity (filename1, weightsname1, filename2, weightsname2, \
        STEPVAL, DELTAVAL, coulombconst, verbose = False, damped = False):

  # read weights file, this first alization is maybe useless
  weights1 = []
  weights2 = []
  weightsfp1 = open(weightsname1)
  weightsfp2 = open(weightsname2)
  weightsfp1.readline()
  idx = 0
  sum = 0.0
  for line in weightsfp1:
    p = re.compile(r'\s+')
    line = p.sub(' ', line)
    line = line.lstrip()
    line = line.rstrip()
                      
    values = re.split(' ', line)

    if len(values) < IDXSLCT+1:
      raise Exception("Error in weight file size")
  
    sum = sum + float(values[IDXSLCT])

    weights1.append(float(values[IDXSLCT]))
    idx = idx + 1

  for i in range(len(weights1)):
    weights1[i] = weights1[i] / sum
  
  sum = 0.0
  weightsfp2.readline()
  idx = 0
  
  for line in weightsfp2:
    values = re.split('\s+', line)

    if len(values) < IDXSLCT+1:
      raise Exception("Error in weight file size")

    sum = sum + float(values[IDXSLCT])
  
    weights2.append(float(values[IDXSLCT]))
    idx = idx + 1
  
  for i in range(len(weights2)):
    weights2[i] = weights2[i] / sum
  
  weightsfp1.close()
  weightsfp2.close()

  mols1 = None
  mols2 = None

  try:
    mols1 = mol2atomextractor (filename1)
    mols2 = mol2atomextractor (filename1)
  except Exception as exp:
    raise Exception(exp)

  if (len(mols1) < 1) or (len(mols2) < 1):
    raise Exception("Error in the number of molecules per file")

  # read data of first molecule
  mol1 = mols1[0]
  atomnum1 = len(mol1)
  mol1coord = numpy.zeros((atomnum1, 3))
  mol1charges = numpy.zeros((atomnum1, 1))
  if verbose:
    print("Molecule: " , filename1)
  idx = 0
  for atom in mol1:
    mol1coord[idx,0] = atom.coords[0]
    mol1coord[idx,1] = atom.coords[1]
    mol1coord[idx,2] = atom.coords[2]
    mol1charges[idx,0] = atom.partialcharge
    idx = idx + 1
  
  # read data of second molecule
  mol2 = mols2[1]
  atomnum2 = len(mol2)
  mol2coord = numpy.zeros((atomnum2, 3))
  mol2charges = numpy.zeros((atomnum2, 1))
  if verbose:
    print("Molecule: " , filename2)
  idx = 0
  for atom in mol2:
    mol2coord[idx,0] = atom.coords[0]
    mol2coord[idx,1] = atom.coords[1]
    mol2coord[idx,2] = atom.coords[2]
    mol2charges[idx,0] = atom.partialcharge
    idx = idx + 1
  
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
  
  if verbose:
    print("Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax)
  
  xnstep = int( ((xmax - xmin) / STEPVAL)+0.5)
  ynstep = int( ((ymax - ymin) / STEPVAL)+0.5)
  znstep = int( ((zmax - zmin) / STEPVAL)+0.5)
  
  # check amount 
  if (len(mols1) != len(weights1)):
    raise Exception("Error different size "+str(len(mol1list))+ \
            " vs "+str(len(weights1)))
  
  if (len(mols2) != len(weights2)):
    raise Exception("Error different size "+str(len(mol2list))+ \
            " vs "+str(len(weights2)))

  if verbose:
    print("Mol1 has %5d "%(len(mols1)) + " values ")
    print("Mol2 has %5d "%(len(mols2)) + " values ")

  if len(mols1) == 0 or len(mols2) == 0:
      raise Exception("Error no element")
  
  mol1field = numpy.zeros((xnstep, ynstep, znstep))
  mol2field = numpy.zeros((xnstep, ynstep, znstep))
  
  xrefpoints = numpy.zeros(xnstep)
  carboidxs = numpy.zeros((len(mols1)*len(mols2), xnstep))
  weights = numpy.zeros(len(mols1)*len(mols2))
  pweights = numpy.zeros(len(mols1)*len(mols2))
  
  # compute weight
  generalidx = 0
  for mol1idx in range(0,len(weights1)):
    for mol2idx in range(0,len(weights2)):
      # compute weights
      weights[generalidx] = weights1[mol1idx] + weights2[mol2idx]
      weights[generalidx] = weights[generalidx] / 2.0
      
      generalidx = generalidx + 1
  
  sum = 0.0
  generalidx = 0
  for mol1idx in range(0,len(weights1)):
    for mol2idx in range(0,len(weights2)):
      # compute weights
      pweights[generalidx] = weights1[mol1idx] * weights2[mol2idx]
      sum = sum + pweights[generalidx]
      
      generalidx = generalidx + 1
  
  for generalidx in range(0,len(pweights)):
    pweights[generalidx] = pweights[generalidx] / sum
  
  if verbose:
    print("")
    print("Start main computation")
    print("")

  generalidx = 0
  for mol1idx in range(0,len(weights1)):
    for mol2idx in range(0,len(weights2)):
  
      # read coord first molecule
      idx = 0
      atomnum1 = len(mols1[mol1idx])
      mol1coord = numpy.zeros((atomnum1, 3))
      mol1charges = numpy.zeros((atomnum1, 1))
      for atom in mols1[mol1idx]:
        mol1coord[idx,0] = atom.coords[0]
        mol1coord[idx,1] = atom.coords[1]
        mol1coord[idx,2] = atom.coords[2]
        mol1charges[idx,0] = atom.partialcharge
        idx = idx + 1
      
      # read data of second molecule
      idx = 0
      atomnum2 = len(mols2[mol2idx])
      mol2coord = numpy.zeros((atomnum2, 3))
      mol2charges = numpy.zeros((atomnum2, 1))
      for atom in mols2[mol2idx]:
        mol2coord[idx,0] = atom.coords[0]
        mol2coord[idx,1] = atom.coords[1]
        mol2coord[idx,2] = atom.coords[2]
        mol2charges[idx,0] = atom.partialcharge
        idx = idx + 1
     
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
            
            sum1 = sum2 = 0.0

            try:
              if damped:
                sum1, sum2 = get_damped_eps (refpoint, mol1coord, mol1charges, \
                  mol2coord, mol2charges, coulombconst)
              else:
                sum1, sum2 = get_eps (refpoint, mol1coord, mol1charges, \
                  mol2coord, mol2charges, coulombconst)
            except Exception as ex:
              raise (ex)

            mol1field[ix,iy,iz] = mol1field[ix,iy,iz] + sum1
            mol2field[ix,iy,iz] = mol2field[ix,iy,iz] + sum2
        
            num = num + sum1*sum2
            denum1 = denum1 + sum1*sum1 
            denum2 = denum2 + sum2*sum2 
        
        if denum1 == 0.0 and denum2 != 0.0:
            carboidx = -1.0
        elif denum2 == 0.0 and denum1 != 0.0:
            carboidx = -1.0
        elif denum2 == 0.0 and denum1 == 0.0:
            carboidx = 1.0
        else:
            carboidx = num/math.sqrt(denum1 * denum2)
        
        if verbose:
          print("%10.5f %5d %5d %10.5f"%(refpoint[0,0], mol1idx, mol2idx, carboidx)) 

        xrefpoints[idx] = refpoint[0,0]
        carboidxs[generalidx,idx] = carboidx
  
        idx = idx + 1

      generalidx = generalidx + 1
  
      if verbose:
        print("Done %5d of %5d "%(generalidx, (len(weights1)*len(weights2))))
        print("")
   
  return (carboidxs, xrefpoints, weights, pweights)


###############################################################

def returncarbodxs(set1, set2, verbose=False):

  g = next(iter(set1.values()))
  xrefpoints = numpy.zeros(g[1].grid.shape[0])
  carboidxs = numpy.zeros((len(set1)*len(set2), g[1].grid.shape[0]))

  generalidx = 0

  for v1 in set1:
    for v2 in set2:
      if verbose:
        print("Compare ", v1, v2)
      
      g1 = set1[v1][1]
      g2 = set2[v2][1]

      try:
        g1.check_compatible(g2)
      except TypeError as te:
        raise Exception(v1 + " and " + v2 + " are not compatible ")

      for i in range(g1.grid.shape[0]):
        x = g1.origin[0] + i*g1.delta[0]

        num = 0.0
        denum1 = 0.0 
        denum2 = 0.0
        
        for j in range(g1.grid.shape[1]):
          y = g1.origin[1] + j*g1.delta[2]
          for k in range(g1.grid.shape[2]):
            z = g1.origin[2] + k*g1.delta[2]

            num = num + g1.grid[i, j, k]*g2.grid[i, j, k]
            denum1 = denum1 + g1.grid[i, j, k]*g1.grid[i, j, k]
            denum2 = denum2 + g2.grid[i, j, k]*g2.grid[i, j, k]
        
        carboidx = 0.0
        if denum1 == 0.0 and denum2 != 0.0:
            carboidx = -1.0
        elif denum2 == 0.0 and denum1 != 0.0:
            carboidx = -1.0
        elif denum2 == 0.0 and denum1 == 0.0:
            carboidx = 1.0
        else:
            carboidx = num/math.sqrt(denum1 * denum2)

        if verbose:
          print("%10.5f %10.5f"%(x, carboidx)) 

        xrefpoints[i] = x
        carboidxs[generalidx,i] = carboidx

      generalidx += 1

  weights = numpy.zeros(len(set1)*len(set2))
  pweights = numpy.zeros(len(set1)*len(set2))
  
  # compute weight
  generalidx = 0
  sum = 0.0
  for v1 in set1:
    for v2 in set2:
      w1 = set1[v1][0]
      w2 = set2[v2][0]
      # compute weights
      weights[generalidx] = w1 + w2
      weights[generalidx] = weights[generalidx] / 2.0
      pweights[generalidx] = w1 * w2
      sum = sum + w1*w2
      generalidx = generalidx + 1
  
  pweights =  pweights / sum
  
  return (carboidxs, xrefpoints, weights, pweights)