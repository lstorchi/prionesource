import re
import os
import sys
import sets
import math
import glob
import numpy
import pybel
import subprocess

import kabsh_minima

from scipy import cluster

CONVERTER = 1.889716

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

    return lines

###############################################################################

def remove_equal(input):
  output = []
  for x in input:
    if x not in output:
      output.append(x)
                        
  return output

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

def read_kont_and_get_centroids (kontname, MINDIM, NUMOFCLUST):

  lineamnt = bufcount(kontname)
  
  dim = (lineamnt - 1)/2
  
  if ((dim * 2) + 1) != lineamnt :
    print "Maybe invalid kont file"
    exit(1)

  print "Reading file ..."
  
  fk = open(kontname)
  
  xsets = sets.Set()
  ysets = sets.Set()
  zsets = sets.Set()
  
  evalset = []
  
  xvals = []
  yvals = []
  zvals = []
  
  ixyz_to_xyzval_map = {}
  xyzval_to_ixyz_map = {}
  
  switchtofieldm = False
  
  energy = numpy.empty([1,1,1], float)
  nx = ny = nz = 0
  ix = iy = iz = 0
  counter = 0
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
       
        if not ( e in evalset):
          evalset.append(e)
  
        ixyzstr = str(ix)+'_'+str(iy)+'_'+str(iz)
        xstr = "{:.3f}".format(xvals[counter])
        ystr = "{:.3f}".format(yvals[counter])
        zstr = "{:.3f}".format(zvals[counter])
        xyzstr = xstr+'_'+ystr+'_'+zstr
        #print ixyzstr , ' ==> ', xyzstr
        counter = counter + 1
        ixyz_to_xyzval_map.update({ixyzstr: xyzstr})
        xyzval_to_ixyz_map.update({xyzstr: ixyzstr})
  
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
  
        xvals.append(float(x))
        yvals.append(float(y))
        zvals.append(float(z))
  
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

  mineval = []
  sortedeval = numpy.sort(evalset)
  
  for i in range(0, min(MINDIM, len(sortedeval))):
    mineval.append(sortedeval[i])
  
  #print mineval
  
  print "Done..."
  
  print "Start computing ... "
  
  dimcube = 5
  
  min_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
  
  for ix in range(0,nx):
    for iy in range(0,ny):
      for iz in range(0,nz):
        min_energy[ix, iy, iz] = 0.0
  
  xmin = []
  ymin = []
  zmin = []
  minvals = []
  
  for ix in range(dimcube,nx-dimcube):
    for iy in range(dimcube,ny-dimcube):
      for iz in range(dimcube,nz-dimcube):
  
        x = botx + ix * dx
        y = boty + iy * dy
        z = botz + iz * dz
  
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
  
  print "Done "
  
  print "Start clustering ..."
  
  pointstocluster = numpy.zeros((len(minvals), 3))
  
  for i in range(0,len(minvals)):
    pointstocluster[i,0] = xmin[i]
    pointstocluster[i,1] = ymin[i]
    pointstocluster[i,2] = zmin[i]
  
  centroids, selected = cluster.vq.kmeans2 (pointstocluster, NUMOFCLUST)
  
  centroidvals = []
  
  for j in range(0, NUMOFCLUST):
    centroidvals.append(0.0)
  
  for j in range(0,len(selected)):
    idx = selected[j]
    centroidvals[idx] += minvals[j] 
  
  centroids_energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)

  return centroids, centroids_energy

###############################################################################

kontname1 = ""
kontname2 = ""
filenamemol21 = ""
filenamemol22 = ""
NUMOFCLUST = 2
MINDIM = 20

print "Per allineare proteine diverse in funzione dei cluster o dei minimi dei campi "

if (len(sys.argv)) == 7: 
  kontname1 = sys.argv[1]
  filenamemol21 = sys.argv[2]
  kontname2 = sys.argv[3]
  filenamemol22 = sys.argv[4]
 
  MINDIM = int(sys.argv[5])
  NUMOFCLUST = int(sys.argv[6])
else:
  print "usage :", sys.argv[0] , " filename1.kont filename1.mol2 filename2.kont filename2.mol2 num_of_min num_of_cluster"
  exit(1)

print "Reading and compute min \"", kontname1, "\""

mol1 = pybel.readfile("mol2", filenamemol21).next()
centroids1, centroids_energy1 = \
    read_kont_and_get_centroids (kontname1, MINDIM, NUMOFCLUST)

print "Done "

print "Reading and compute min \"", kontname2, "\""

mol2 = pybel.readfile("mol2", filenamemol22).next()
centroids2, centroids_energy2 = \
    read_kont_and_get_centroids (kontname2, MINDIM, NUMOFCLUST)

print "Done "

print "Kabsh minim..."

rmatrix, translate1, translate2 = \
    kabsh_minima.return_rotation_matrix(centroids1, centroids2)

#print rmatrix, translate1, translate2
# da rivedere 

ormatrix = pybel.ob.matrix3x3()
for i in range(3):
  for j in range(3):
    ormatrix.Set(i, j, rmatrix[i,j]) 
myrm = pybel.ob.doubleArray(9)
ormatrix.GetArray(myrm)

mol1.OBMol.Translate(pybel.ob.vector3(translate1[0], translate1[1], translate1[2]))
mol1.OBMol.Rotate(myrm)
mol2.OBMol.Translate(pybel.ob.vector3(translate2[0], translate2[1], translate2[2]));


###apply rotation to atoms coordinates R * v
