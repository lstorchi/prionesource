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
    self.id = 0
    self.element = ""
    self.radii = 0.0
    self.residueid = 0
  
  def __repr__ (self):
    line = "%6d %6s %6s %10.5f %10.5f %10.5f %8.5f %3s\n"%( \
      self.id, self.name, self.resname, self.coords[0], \
      self.coords[1], self.coords[2], self.partialcharge, \
      self.atomname)

    return line
  
###############################################################

def get_seps (refpoint, molcoord, molcharges, coulombconst, \
  ddielectric):
            
  #print refpoint[0,0], refpoint[0,1], refpoint[0,2], 1.0
  dist = scipy.spatial.distance.cdist(molcoord,refpoint)

  sum = 0.0

  #print("TODO")
  if (dist.min() > 1.0):
    if ddielectric:
      # simplest distance-dependent dielectric constant has been implemented in lammps
      ep = coulombconst * (molcharges/(dist**2))
    else:
      ep = coulombconst * (molcharges/dist)

    sum = numpy.sum(ep) 

  return sum

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

def pdbatomextractor (file=None):

  from Bio.PDB.PDBParser import PDBParser
  
  parser=PDBParser(PERMISSIVE=1)
  
  structure = parser.get_structure("singelemodel", file)

  mols = []

  for model in structure:

    id = 1
    mol = []
    for chain in model:
      for residue in chain:
        for a in residue:
          coords = a.get_coord()
          #print(id, residue.get_resname(), a.get_name(), coords)
          na = atom(id, a.get_name(), coords[0], coords[1], coords[2], 0.0 )
          na.id = id
          na.resname = residue.get_resname()
          na.residueid = residue.get_id()[1]
          na.atomname = a.get_name()
          na.element = a.element
          #na.radii = mendeleev.element(a.element).vdw_radius

          mol.append(na)
          id += 1

    mols.append(mol)

  return mols

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

def returncarbodxs(set1, set2, verbose=False, axis="x"):

  g = next(iter(set1.values()))
 
  generalidx = 0

  aidx1 = 0
  aidx2 = 1
  aidx3 = 2

  if axis == "x":
    aidx1 = 0
    aidx2 = 1
    aidx3 = 2
  elif axis == "y":
    aidx1 = 1
    aidx2 = 0
    aidx3 = 2
  elif axis == "z":
    aidx1 = 2
    aidx2 = 0
    aidx3 = 1

  xrefpoints = numpy.zeros(g[1].grid.shape[aidx1])
  carboidxs = numpy.zeros((len(set1)*len(set2), g[1].grid.shape[aidx1]))

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

      x = y = z = 0.0
      ix = iy = iz = 0
      movcoord = 0.0

      for i in range(g1.grid.shape[aidx1]):
        if aidx1 == 0:
          ix = i
          x = g1.origin[aidx1] + i*g1.delta[aidx1]
          movcoord = x
        elif aidx1 == 1:
          iy = i
          y = g1.origin[aidx1] + i*g1.delta[aidx1]
          movcoord = y
        elif aidx1 == 2:
          iz = i
          z = g1.origin[aidx1] + i*g1.delta[aidx1]
          movcoord = z

        num = 0.0
        denum1 = 0.0 
        denum2 = 0.0

        for j in range(g1.grid.shape[aidx2]):
          if aidx2 == 0:
            ix = j
            x = g1.origin[aidx2] + j*g1.delta[aidx2]
          elif aidx2 == 1:
            iy = j
            y = g1.origin[aidx2] + j*g1.delta[aidx2]
          elif aidx2 == 2:
            iz = j
            z = g1.origin[aidx2] + j*g1.delta[aidx2]

          for k in range(g1.grid.shape[aidx3]):
            if aidx3 == 0:
              ix = k
              x = g1.origin[aidx3] + k*g1.delta[aidx3]
            elif aidx3 == 1:
              iy = k
              y = g1.origin[aidx3] + k*g1.delta[aidx3]
            elif aidx3 == 2:
              iz = k
              z = g1.origin[aidx3] + k*g1.delta[aidx3]

            
            num = num + g1.grid[ix, iy, iz]*g2.grid[ix, iy, iz]
            denum1 = denum1 + g1.grid[ix, iy, iz]*g1.grid[ix, iy, iz]
            denum2 = denum2 + g2.grid[ix, iy, iz]*g2.grid[ix, iy, iz]
        
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
          print("%10.5f %10.5f"%(movcoord, carboidx)) 

        xrefpoints[i] = movcoord
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

#################################################################################

def get_cfields (mols, STEPVAL, DELTAVAL, coulombconst, verbose = False, \
     ddielectric=False ): 

  # read data of first molecule
  mol = mols[0]
  atomnum = len(mol)
  molcoord = numpy.zeros((atomnum, 3))
  molcharges = numpy.zeros((atomnum, 1))
  for idx, atom in enumerate(mol):
    molcoord[idx,0] = atom.coords[0]
    molcoord[idx,1] = atom.coords[1]
    molcoord[idx,2] = atom.coords[2]
    molcharges[idx,0] = atom.partialcharge
  
  # generate grid 
  xmin = min(molcoord[:,0])
  ymin = min(molcoord[:,1])
  zmin = min(molcoord[:,2])
  xmax = max(molcoord[:,0])
  ymax = max(molcoord[:,1])
  zmax = max(molcoord[:,2])
  xmin = xmin - DELTAVAL
  ymin = ymin - DELTAVAL
  zmin = zmin - DELTAVAL
  xmax = xmax + DELTAVAL
  ymax = ymax + DELTAVAL
  zmax = zmax + DELTAVAL
  
  if verbose:
    print("Grid will be used: %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"%( \
      xmin, ymin, zmin, xmax, ymax, zmax))
  
  xnstep = int( ((xmax - xmin) / STEPVAL)+0.5)
  ynstep = int( ((ymax - ymin) / STEPVAL)+0.5)
  znstep = int( ((zmax - zmin) / STEPVAL)+0.5)
  
  molsfield = []
  
  if verbose:
    print("")
    print("Start main computation")
    print("")

  for molidx in range(len(mols)):
    if verbose:
      print("Computing mol ", molidx+1)
    # read coord first molecule
    atomnum = len(mols[molidx])
    molcoord = numpy.zeros((atomnum, 3))
    molcharges = numpy.zeros((atomnum, 1))
    totcharge = 0.0
    for idx, atom in enumerate(mols[molidx]):
      molcoord[idx,0] = atom.coords[0]
      molcoord[idx,1] = atom.coords[1]
      molcoord[idx,2] = atom.coords[2]
      molcharges[idx,0] = atom.partialcharge
      totcharge += atom.partialcharge
    if verbose:
      print("   totalcharge: ", totcharge)
    
    molfield = numpy.zeros((xnstep, ynstep, znstep))
    coords = []
    refpoint = numpy.zeros((1, 3))
    for ix in range(0,xnstep):
      refpoint[0,0] = xmin + ix*STEPVAL
      for iy in range(0,ynstep):
        refpoint[0,1] = ymin + iy*STEPVAL
        for iz in range(0,znstep):
          refpoint[0,2] = zmin + iz*STEPVAL
          molfield[ix,iy,iz] = get_seps (refpoint, molcoord, \
            molcharges, coulombconst, ddielectric)
          coords.append((refpoint[0][0], refpoint[0][1], refpoint[0][2])) 

    molsfield.append((coords, molfield))

  return molsfield