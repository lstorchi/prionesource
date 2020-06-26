import scipy.spatial
import numpy 
import pybel
import math
import sys
import re

from scipy import cluster

CONVERTER = 1.889716
STEPVAL = 1.0*CONVERTER
DELTAVAL = 10.0*CONVERTER
coulombconst = 1.0

#####################################################################

def mean_elekfield (probecharge, filename, weifilename, DELTAVAL, \
        STEPVAL):

  # read weights file, this first nomalization is maybe useless
  weights = []
  weightsfp = open(weifilename)
  weightsfp.readline()
  idx = 0
  sum = 0.0
  for line in weightsfp:
    values = re.split('\s+', line)
  
    sum = sum + float(values[6])
  
    weights.append(float(values[6]))
    idx = idx + 1
  
  for i in range(len(weights)):
    weights[i] = weights[i] / sum
  
  weightsfp.close()
  
  # read first molecule in the files
  mol = next(pybel.readfile("mol2", filename))
  
  # read data of first molecule
  atomnum = len(mol.atoms)
  molcoord = numpy.zeros((atomnum, 3))
  molcharges = numpy.zeros((atomnum, 1))
  idx = 0
  for atom in mol:
    molcoord[idx,0] = atom.coords[0]
    molcoord[idx,1] = atom.coords[1]
    molcoord[idx,2] = atom.coords[2]
    molcharges[idx,0] = atom.partialcharge
    idx = idx + 1
  
  # generate grid 
  xmin = min(min(molcoord[:,0]),min(molcoord[:,0]))
  ymin = min(min(molcoord[:,1]),min(molcoord[:,1]))
  zmin = min(min(molcoord[:,2]),min(molcoord[:,2]))
  xmax = max(max(molcoord[:,0]),max(molcoord[:,0]))
  ymax = max(max(molcoord[:,1]),max(molcoord[:,1]))
  zmax = max(max(molcoord[:,2]),max(molcoord[:,2]))
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
  
  # read all molecules 
  mollist = list(pybel.readfile("mol2", filename))
  
  # check amount 
  if (len(mollist) != len(weights)):
    print("Error different size")
    exit()
  
  molfield = numpy.zeros((xnstep, ynstep, znstep))
  
  for molidx in range(0,len(weights)):
    print("Step ", molidx+1, " of ", len(weights))
    # read coord first molecule
    idx = 0
    atomnum = len(mollist[molidx].atoms)
    molcoord = numpy.zeros((atomnum, 3))
    molcharges = numpy.zeros((atomnum, 1))
    for atom in mollist[molidx]:
      molcoord[idx,0] = atom.coords[0]
      molcoord[idx,1] = atom.coords[1]
      molcoord[idx,2] = atom.coords[2]
      molcharges[idx,0] = atom.partialcharge
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
          
          dist = scipy.spatial.distance.cdist(molcoord,refpoint)
          sum = 0.0
      
          if (dist.min() > 1.0):
            ep1 = coulombconst * ((probecharge*molcharges)/dist)
            sum = numpy.sum(ep1)
      
          molfield[ix,iy,iz] = molfield[ix,iy,iz] + sum*weights[molidx]

  xset = []
  yset = []
  zset = []
  eset = []

  ixyz_to_xyzval_map = {}
  xyzval_to_ixyz_map = {}

  refpoint[0,0] = xmin - STEPVAL
  for ix in range(0,xnstep):
    refpoint[0,0] = refpoint[0,0] + STEPVAL
    refpoint[0,1] = ymin - STEPVAL
    xset.append(refpoint[0,0])
    
    for iy in range(0,ynstep):
      refpoint[0,1] = refpoint[0,1] + STEPVAL
      refpoint[0,2] = zmin - STEPVAL
      if (ix == 0):
        yset.append(refpoint[0,1])

      for iz in range(0,znstep):
        refpoint[0,2] = refpoint[0,2] + STEPVAL
        if (ix == 0) and (iy == 0):
          zset.append(refpoint[0,2])

        ixyzstr = str(ix)+'_'+str(iy)+'_'+str(iz)
        xstr = "{:.3f}".format(refpoint[0,0])
        ystr = "{:.3f}".format(refpoint[0,1])
        zstr = "{:.3f}".format(refpoint[0,2])
        xyzstr = xstr+'_'+ystr+'_'+zstr
        ixyz_to_xyzval_map.update({ixyzstr: xyzstr})
        xyzval_to_ixyz_map.update({xyzstr: ixyzstr})
   
        eset.append(molfield[ix,iy,iz])

  return molfield, xset, yset, zset, eset, ixyz_to_xyzval_map, \
          xyzval_to_ixyz_map

###############################################################################

def elekfield (probecharge, mol1coord, mol1charges, xmin, ymin, \
        zmin, xmax, ymax, zmax, STEPVAL):

  xset = []
  yset = []
  zset = []
  eset = []
    
  xnstep = int((xmax - xmin)/STEPVAL) + 1
  ynstep = int((ymax - ymin)/STEPVAL) + 1
  znstep = int((zmax - zmin)/STEPVAL) + 1

  print(xnstep, ynstep, znstep)

  mol1field = numpy.zeros((xnstep, ynstep, znstep))

  ixyz_to_xyzval_map = {}
  xyzval_to_ixyz_map = {}

  idx = 0
  refpoint = numpy.zeros((1, 3))
  refpoint[0,0] = xmin - STEPVAL
  for ix in range(0,xnstep):
    refpoint[0,0] = refpoint[0,0] + STEPVAL
    refpoint[0,1] = ymin - STEPVAL
    xset.append(refpoint[0,0])
    
    for iy in range(0,ynstep):
      refpoint[0,1] = refpoint[0,1] + STEPVAL
      refpoint[0,2] = zmin - STEPVAL
      if (ix == 0):
        yset.append(refpoint[0,1])

      for iz in range(0,znstep):
        refpoint[0,2] = refpoint[0,2] + STEPVAL
        if (ix == 0) and (iy == 0):
          zset.append(refpoint[0,2])

        ixyzstr = str(ix)+'_'+str(iy)+'_'+str(iz)
        xstr = "{:.3f}".format(refpoint[0,0])
        ystr = "{:.3f}".format(refpoint[0,1])
        zstr = "{:.3f}".format(refpoint[0,2])
        xyzstr = xstr+'_'+ystr+'_'+zstr
        ixyz_to_xyzval_map.update({ixyzstr: xyzstr})
        xyzval_to_ixyz_map.update({xyzstr: ixyzstr})
        
        #print refpoint[0,0], refpoint[0,1], refpoint[0,2], 1.0
        dist1 = scipy.spatial.distance.cdist(mol1coord,refpoint)
        sum1 = 0.0
    
        if (dist1.min() > 1.0):
          # per essere certo che non sono interno posso valuate se ho 
          # atomi a destra sinistra basso alto ... 
          ep1 = coulombconst * ((probecharge*mol1charges)/dist1)
          sum1 = numpy.sum(ep1)
    
        mol1field[ix,iy,iz] = mol1field[ix,iy,iz] + sum1

  for ix in range(0,xnstep):
    for iy in range(0,ynstep):
      for iz in range(0,znstep):
        eset.append(mol1field[ix,iy,iz])

  return mol1field, xset, yset, zset, eset, ixyz_to_xyzval_map, \
          xyzval_to_ixyz_map

#####################################################################

def computelek_and_get_centroids (filename, probecharge, STEPVAL, \
        MINDIM, NUMOFCLUST, DELTAVAL):

  mol = next(pybel.readfile("mol2", filename))

  atomnum = len(mol.atoms)
  mol1coord = numpy.zeros((atomnum, 3))
  mol1charges = numpy.zeros((atomnum, 1))
  idx = 0
  for atom in mol:
    mol1coord[idx,0] = atom.coords[0]
    mol1coord[idx,1] = atom.coords[1]
    mol1coord[idx,2] = atom.coords[2]
    mol1charges[idx,0] = atom.partialcharge
    idx = idx + 1

  # generate grid 
  botx = min(min(mol1coord[:,0]),min(mol1coord[:,0])) - DELTAVAL
  boty = min(min(mol1coord[:,1]),min(mol1coord[:,1])) - DELTAVAL
  botz = min(min(mol1coord[:,2]),min(mol1coord[:,2])) - DELTAVAL
  topx = max(max(mol1coord[:,0]),max(mol1coord[:,0])) + DELTAVAL
  topy = max(max(mol1coord[:,1]),max(mol1coord[:,1])) + DELTAVAL
  topz = max(max(mol1coord[:,2]),max(mol1coord[:,2])) + DELTAVAL

  print("Computing elek field ...")

  energy, xset, yset, zset, evalset, ixyz_to_xyzval_map, \
          xyzval_to_ixyz_map = elekfield (probecharge, mol1coord, \
          mol1charges, botx, boty, botz, topx, topy, topz, STEPVAL)

  print("Done")

  dx = STEPVAL
  dy = STEPVAL
  dz = STEPVAL

  mineval = []
  sortedeval = numpy.sort(evalset)

  #y, x = numpy.histogram(evalset, MINDIM)
  
  #for i in range(0,len(y)):
  #  print x[i], y[i]

  #exit(1)

  for i in range(0, min(MINDIM, len(sortedeval))):
    mineval.append(sortedeval[i])
  
  #print mineval
  
  dimcube = 5

  nx = len(xset)
  ny = len(yset)
  nz = len(zset)

  #print nx, ny, nz
  
  min_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
  
  for ix in range(0,nx):
    for iy in range(0,ny):
      for iz in range(0,nz):
        min_energy[ix, iy, iz] = 0.0
  
  xmin = []
  ymin = []
  zmin = []
  minvals = []

  print("Computing minima...")
  
  for ix in range(dimcube,nx-dimcube):
    for iy in range(dimcube,ny-dimcube):
      for iz in range(dimcube,nz-dimcube):
  
        x = xset[ix]
        y = yset[iy]
        z = zset[iz]
  
        eref = energy[ix, iy, iz]
  
        notminima = False
        
        if (eref in mineval):
  
          for ix_near in range(-dimcube,dimcube):
            for iy_near in range(-dimcube,dimcube):
              for iz_near in range(-dimcube,dimcube):
          
                if not ((ix_near == 0) and (iy_near == 0) and (iz_near == 0)):
          
                  x_near = x + (ix_near * dx)
                  y_near = y + (iy_near * dy)
                  z_near = z + (iz_near * dz)
          
                  xstr = "{:.3f}".format(x)
                  ystr = "{:.3f}".format(y)
                  zstr = "{:.3f}".format(z)
                  xyzstr = xstr+'_'+ystr+'_'+zstr
                  ixyz = xyzval_to_ixyz_map[xyzstr].split("_")
                  near_ix = int(ixyz[0])
                  near_iy = int(ixyz[1])
                  near_iz = int(ixyz[2])
  
                  e = energy[near_ix, near_iy, near_iz]
                  
                  if (e < eref):
                    notminima = True
          
          if (notminima):
            min_energy[ix, iy, iz] = 0.0
          else:
            #print "{:.3f}".format(x), " {:.3f}".format(y), " {:.3f}".format(z), " {:.3f}".format(eref)
            min_energy[ix, iy, iz] = eref
            minvals.append(eref)
            xmin.append(x)
            ymin.append(y)
            zmin.append(z)
  
  print("Done ")
  
  print("Start clustering ...")
  
  pointstocluster = numpy.zeros((len(minvals), 3))
  
  #print len(minvals)
  #print ""
  for i in range(0,len(minvals)):
    pointstocluster[i,0] = xmin[i]
    pointstocluster[i,1] = ymin[i]
    pointstocluster[i,2] = zmin[i]
    #print "H   ", xmin[i], ymin[i], zmin[i]
  
  #centroids, selected = cluster.vq.kmeans2 (pointstocluster, NUMOFCLUST, \
  #        iter=200, thresh=1e-9)

  inpoint = cluster.vq.whiten (pointstocluster)
  centroids, selected = cluster.vq.kmeans (inpoint, NUMOFCLUST, \
          iter=100, thresh=1e-9)

  '''
  clustermass = []
  for i in range(0, NUMOFCLUST):
    clustermass.append(0)
    for j in range(0,len(selected)):
      if (selected[j] == i):
        clustermass[i] = clustermass[i] + 1

  print clustermass
  ''' 

  print("Done.") 

  #print centroids
  #print len(centroids)
  #print ""
  #for i in range(0,len(centroids)):
  #  print "O   ", centroids[i][0], centroids[i][1], centroids[i][2]

  #exit(1)
  
  '''
  centroidvals = []
  
  for j in range(0, NUMOFCLUST):
    centroidvals.append(0.0)
  
  for j in range(0,len(selected)):
    idx = selected[j]
    centroidvals[idx] += minvals[j] 
  
  centroids_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)

  return centroids, centroids_energy
  '''

  return centroids

  #return pointstocluster

###############################################################################

def compute_meanelek_and_get_centroids (filename, weifilename, probecharge, \
        STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL):

  print("Computing mean elek field ...")

  energy, xset, yset, zset, evalset, ixyz_to_xyzval_map, \
          xyzval_to_ixyz_map = mean_elekfield (probecharge, filename, \
          weifilename, DELTAVAL, STEPVAL)

  print("Done")

  dx = STEPVAL
  dy = STEPVAL
  dz = STEPVAL

  mineval = []
  sortedeval = numpy.sort(evalset)

  for i in range(0, min(MINDIM, len(sortedeval))):
    mineval.append(sortedeval[i])
  
  dimcube = 1

  nx = len(xset)
  ny = len(yset)
  nz = len(zset)

  min_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
  
  for ix in range(0,nx):
    for iy in range(0,ny):
      for iz in range(0,nz):
        min_energy[ix, iy, iz] = 0.0
  
  xmin = []
  ymin = []
  zmin = []
  minvals = []

  print("Computing minima...")
  
  for ix in range(dimcube,nx-dimcube):
    for iy in range(dimcube,ny-dimcube):
      for iz in range(dimcube,nz-dimcube):
  
        x = xset[ix]
        y = yset[iy]
        z = zset[iz]
  
        eref = energy[ix, iy, iz]
  
        notminima = False
        
        if (eref in mineval):
  
          for ix_near in range(-dimcube,dimcube):
            for iy_near in range(-dimcube,dimcube):
              for iz_near in range(-dimcube,dimcube):
          
                if not ((ix_near == 0) and (iy_near == 0) and (iz_near == 0)):
          
                  x_near = x + (ix_near * dx)
                  y_near = y + (iy_near * dy)
                  z_near = z + (iz_near * dz)
          
                  xstr = "{:.3f}".format(x)
                  ystr = "{:.3f}".format(y)
                  zstr = "{:.3f}".format(z)
                  xyzstr = xstr+'_'+ystr+'_'+zstr
                  ixyz = xyzval_to_ixyz_map[xyzstr].split("_")
                  near_ix = int(ixyz[0])
                  near_iy = int(ixyz[1])
                  near_iz = int(ixyz[2])
  
                  e = energy[near_ix, near_iy, near_iz]
                  
                  if (e < eref):
                    notminima = True
          
          if (notminima):
            min_energy[ix, iy, iz] = 0.0
          else:
            #print "{:.3f}".format(x), " {:.3f}".format(y), " {:.3f}".format(z), " {:.3f}".format(eref)
            min_energy[ix, iy, iz] = eref
            minvals.append(eref)
            xmin.append(x)
            ymin.append(y)
            zmin.append(z)
  
  print("Done ")
  
  print("Start clustering ...")
  
  pointstocluster = numpy.zeros((len(minvals), 3))
  
  for i in range(0,len(minvals)):
    pointstocluster[i,0] = xmin[i]
    pointstocluster[i,1] = ymin[i]
    pointstocluster[i,2] = zmin[i]
    print("H   ", xmin[i], ymin[i], zmin[i])

  #return pointstocluster
  # se faccio questo i centroid non sono piu' ciorretti nelle coordinate 
  #inpoint = cluster.vq.whiten (pointstocluster)

  if (len(minvals) <= NUMOFCLUST):
    print("Num of cluster is >= of minvals")
    if (len(minvals) == NUMOFCLUST):
      return pointstocluster
    exit(1)

  centroids, selected = cluster.vq.kmeans2 (pointstocluster, NUMOFCLUST)
          #iter=100, thresh=1e-9)

  print("Done.") 

  return centroids

#####################################################################
