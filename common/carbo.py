import scipy.spatial
import numpy 
import math
import re
import os

from rdkit import Chem
from rdkit.Chem import Draw,AllChem

IDXSLCT = 5

###############################################################

class atom:

  def __init__ (self, id, name, x, y, z, c):
    self.id = id
    self.name = name
    self.coords = (x, y, z)
    self.partialcharge = c
  
  def __repr__ (self):
    line = "%6d %6s %10.5f %10.5f %10.5f %8.5f \n"%( \
      self.id, self.name, self.coords[0], self.coords[1],
      self.coords[2], self.partialcharge )

    return line
  
###############################################################

def mol2atomextractor (file=None):
  
  mols = []

  with open(file, 'r') as f:
      while not f.tell() == os.fstat(f.fileno()).st_size:
        line = f.readline()
        if line.startswith("@<TRIPOS>ATOM"):
            line = f.readline()
            mol = []
            while not line.startswith("@<TRIPOS>BOND"):
              sline = line.split()
              if len(sline) != 9:
                print("Error in ", line)
                exit(1)
              a = atom(int(sline[0]), sline[1], float(sline[2]), \
                float(sline[3]),float(sline[4]), float(sline[8]) )
              line = f.readline()
              mol.append(a)
            mols.append(mol)

  return(mols)

###############################################################

def carbo_similarity (filename1, weightsname1, filename2, weightsname2, \
        STEPVAL, DELTAVAL, coulombconst):

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
        print("Error in weight file size")
        return 
  
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
        print("Error in weight file size")
        return 

    sum = sum + float(values[IDXSLCT])
  
    weights2.append(float(values[IDXSLCT]))
    idx = idx + 1
  
  for i in range(len(weights2)):
    weights2[i] = weights2[i] / sum
  
  weightsfp1.close()
  weightsfp2.close()
  
  mols1 = mol2atomextractor (filename1)
  mols2 = mol2atomextractor (filename1)

  # read data of first molecule
  mol1 = mols1[0]
  atomnum1 = len(mol1)
  mol1coord = numpy.zeros((atomnum1, 3))
  mol1charges = numpy.zeros((atomnum1, 1))
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
  
  print("Grid will be used: ", xmin, ymin, zmin, xmax, ymax, zmax)
  
  xnstep = int( ((xmax - xmin) / STEPVAL)+0.5)
  ynstep = int( ((ymax - ymin) / STEPVAL)+0.5)
  znstep = int( ((zmax - zmin) / STEPVAL)+0.5)
  
  # check amount 
  if (len(mols1) != len(weights1)):
    print("Error different size ", len(mol1list) , " vs ", len(weights1))
    exit()
  
  if (len(mols2) != len(weights2)):
    print("Error different size", len(mol2list) , " vs ", len(weights2))
    exit()

  print("Mol1 has ", len(mols1), " values ")
  print("Mol2 has ", len(mols2), " values ")

  if len(mols1) == 0 or len(mols2) == 0:
      print("Error no element")
      return 
  
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
            
            #print refpoint[0,0], refpoint[0,1], refpoint[0,2], 1.0
            dist1 = scipy.spatial.distance.cdist(mol1coord,refpoint)
            dist2 = scipy.spatial.distance.cdist(mol2coord,refpoint)
        
            sum1 = 0.0
            sum2 = 0.0
        
            if (dist1.min() > 1.0):
              # per essere certo che non sono interno posso valuate se ho 
              # atomi a destra sinistra basso alto ... 
              ep1 = coulombconst * (mol1charges/dist1)
              sum1 = numpy.sum(ep1) 
            #per valurae ad esempio se non prendo i punti interni 
            #else:
            #  print refpoint[0,0], refpoint[0,1], refpoint[0,2], sum1
        
            if (dist2.min() > 1.0):
              ep2 = coulombconst * (mol2charges/dist2)
              sum2 = numpy.sum(ep2)
      
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
        
  
        print(refpoint[0,0], mol1idx, mol2idx, carboidx) 
        xrefpoints[idx] = refpoint[0,0]
        carboidxs[generalidx,idx] = carboidx
  
        idx = idx + 1
      generalidx = generalidx + 1
  
      print("Done ", generalidx, " of ", (len(weights1)*len(weights2)))
  
  stdev = carboidxs.std(0)
  meanmtx = carboidxs.mean(0)
  
  print("Mead and stdev")
  idx = 0
  for std in stdev:
    print(xrefpoints[idx], " ", meanmtx[idx] , " ", std)
    idx = idx + 1
  
  print("full matrix")
  print(carboidxs.mean(0))
  print(carboidxs.std(0))
  
  print("Weighted Mead and stdev")
  waverage = numpy.average(carboidxs, 0, weights)
  wvariance = numpy.average((carboidxs-waverage)**2, 0, weights)
  idx = 0
  for std in wvariance:
    print(xrefpoints[idx], " ", waverage[idx] , " ", std)
    idx = idx + 1
  
  print("full matrix")
  print(waverage)
  print(wvariance)
  
  print("PWeighted Mead and stdev")
  waverage = numpy.average(carboidxs, 0, pweights)
  wvariance = numpy.average((carboidxs-waverage)**2, 0, pweights)
  idx = 0
  for std in wvariance:
    print(xrefpoints[idx], " ", waverage[idx] , " ", std)
    idx = idx + 1
  
  print("full matrix")
  print(waverage)
  print(wvariance)
