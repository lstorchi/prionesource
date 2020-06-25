import re
import os
import sys
import sets
import math
import glob
import numpy
import pybel
import scipy
import subprocess

sys.path.append("./common")
import elektrofield
import kabsch_minima

CONVERTER = 1.889716
STEPVAL = 2.0
DELTAVAL = 5.0

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

def order_via_eucldist (centroids1, centroids2):

  dists = scipy.spatial.distance.cdist(centroids1, centroids2)
  
  idx1toidx2 = {}
  mindist = []
  idx = 0
  print "idx1 idx2 dist.min."
  for d in dists:
    mv = min(d)
    mindist.append(mv)
    i = numpy.argwhere(d == mv)
    if len(i) > 1:
      print "ERROR internal error"
      exit(1)
      if (len(i[0])) > 1:
        print "ERROR not a single min val"
        exit(1)
    idx1toidx2[idx] = i[0][0]
    print idx, " ", i[0][0], " dist ", mv
    idx = idx + 1
  
  mindist_sort = numpy.sort(mindist)
  print " "
  print "Min. dists"
  print mindist_sort
  print " "
  
  idxslct1 = []
  idxslct2 = []
  
  print "Selected pairs"
  print "idx1 idx2 dist"
  for d in mindist_sort:
    i = numpy.argwhere(mindist == d)
    if len(i) > 1:
      print "ERROR internal error"
      exit(1)
      if (len(i[0])) > 1:
        print "ERROR not a single min val"
        exit(1)
   
    counter = 0
    for k, v in idx1toidx2.iteritems():
      if counter == i[0][0]:
        if (k not in idxslct1) and (v not in idxslct2):
          print k, " ", v, " ", d
          idxslct1.append(k)
          idxslct2.append(v)
          break
      counter = counter + 1
  
  centroids1_sort = numpy.empty ((len(idxslct1), 3))
  centroids2_sort = numpy.empty ((len(idxslct2), 3))
  
  i = 0
  for idx in idxslct1:
    centroids1_sort[i] = centroids1[idx]
    i = i + 1
  
  i = 0
  for idx in idxslct2:
    centroids2_sort[i] = centroids2[idx]
    i = i + 1

  return centroids1_sort, centroids2_sort

###############################################################################

def order_via_distfromcenter (centroids1, centroids2):

  translate1 = kabsch_minima.centroid(centroids1)
  translate2 = kabsch_minima.centroid(centroids2)

  c1 = numpy.empty_like (centroids1)

  for i in range(0,len(centroids1)):
    for j in range(0,3):
      c1[i][j] = centroids1[i][j] - translate1[j]

  dist1 = []
  for i in range(0,len(c1)):
    sum = 0.0
    for j in range(0,3):
      sum = sum + c1[i][j]**2

    d = math.sqrt(sum)
    dist1.append(d)

  c2 = numpy.empty_like (centroids2)

  for i in range(0,len(centroids2)):
    for j in range(0,3):
      c2[i][j] = centroids2[i][j] - translate2[j]

  dist2 = []
  for i in range(0,len(c2)):
    sum = 0.0
    for j in range(0,3):
      sum = sum + c1[i][j]**2

    d = math.sqrt(sum)
    dist2.append(d)

  centroids1_sort = numpy.empty_like (centroids1)
  centroids1_sort[:] = centroids1[:]

  for i in range(len(dist1)):
    for j in range(len(dist1)-1-i):
      if dist1[j] > dist1[j+1]:
        dist1[j], dist1[j+1] = dist1[j+1], dist1[j] 

        swapper = []
        for k in range(0,3):
            swapper.append(centroids1_sort[j][k])
        centroids1_sort[j][:] = centroids1_sort[j+1][:]
        for k in range(0,3):
          centroids1_sort[j+1][k] = swapper[k]


  #print dist1
  #print centroids1_sort

  centroids2_sort = numpy.empty_like (centroids2)
  centroids2_sort[:] = centroids2[:]

  for i in range(len(dist2)):
    for j in range(len(dist2)-1-i):
      if dist2[j] > dist2[j+1]:
        dist2[j], dist2[j+1] = dist2[j+1], dist2[j] 

        swapper = []
        for k in range(0,3):
            swapper.append(centroids2_sort[j][k])
        centroids2_sort[j][:] = centroids2_sort[j+1][:]
        for k in range(0,3):
          centroids2_sort[j+1][k] = swapper[k]

  #print dist2
  #print centroids2_sort

  return centroids1_sort, centroids2_sort
  
###############################################################################

filenamemol21 = ""
filenamemol22 = ""
# per determinare questi valori ho usato 
# python fileds_minima_superposition.py h_irf1_1.mol2 h_irf1_1.mol2 50 5 
# volendo chiaramente sempre ottenere N centroidi tali da  avere sempre almeno 
# lo stesso nel secondo run (essendo le due strutture sempre la stessa) 
# da capire ma c'e' una certa randomicicita' che porta ad una certa instabilita'
NUMOFCLUST = 5
MINDIM = 50

#print "Per allineare proteine diverse in funzione dei cluster o dei minimi dei campi "

if (len(sys.argv)) == 3: 
  filenamemol21 = sys.argv[1]
  filenamemol22 = sys.argv[2]
else:
  print "usage :", sys.argv[0] , " filename1.mol2 filename2.mol2"
  exit(1)

fileout1 = filenamemol21[:-5]+"_rt.mol2"
fileout2 = filenamemol22[:-5]+"_rt.mol2"

print "Compute fields and minima \"", filenamemol21, "\""

centroids1_1 = elektrofield.computelek_and_get_centroids (\
        filenamemol21, -1.0, STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL)
centroids1_2 = elektrofield.computelek_and_get_centroids (\
        filenamemol21, 1.0, STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL)
centroids1 = numpy.append(centroids1_1, centroids1_2, axis=0)

#centroids1 = elektrofield.computelek_and_get_centroids (\
#        filenamemol21, -1.0, STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL)

print centroids1

print "Done "

print "Compute fields and minima \"", filenamemol22, "\""

centroids2_1 = elektrofield.computelek_and_get_centroids (\
        filenamemol22, -1.0, STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL)
centroids2_2 = elektrofield.computelek_and_get_centroids (\
        filenamemol22, 1.0, STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL)
centroids2 = numpy.append(centroids2_1, centroids2_2, axis=0)

#centroids2 = elektrofield.computelek_and_get_centroids (\
#        filenamemol22, -1.0, STEPVAL, MINDIM, NUMOFCLUST, DELTAVAL)

print centroids2

print "Done "

print "Order and select centroids... "

centroids1_sort, centroids2_sort = order_via_eucldist (centroids1, centroids2)
#centroids1_sort, centroids2_sort = order_via_distfromcenter (centroids1, centroids2)

#centroids1_sort = numpy.empty_like (centroids1)
#centroids1_sort[:] = centroids1[:]

#centroids2_sort = numpy.empty_like (centroids2)
#centroids2_sort[:] = centroids2[:]

print "Selected and sorted centroids1 ..."
print centroids1_sort
print "Selected and sorted centroids2 ..."
print centroids2_sort
print "Check dists "
dists = scipy.spatial.distance.cdist(centroids1_sort, centroids2_sort)
print dists

print "RMSD: ", kabsch_minima.rmsd(centroids1_sort, centroids2_sort)

o_file = open("centroids1.xyz","w")
o_file.write(str(len(centroids1_sort))+"\n")
o_file.write("centroids 1 \n")
for i in range(0,len(centroids1_sort)):
  o_file.write("O   "+str(centroids1_sort[i][0])+" "+\
          str(centroids1_sort[i][1])+" "+\
          str(centroids1_sort[i][2])+"\n")
o_file.close()

o_file = open("centroids2.xyz","w")
o_file.write(str(len(centroids2_sort))+"\n")
o_file.write("centroids 2 \n")
for i in range(0,len(centroids2_sort)):
  o_file.write("O   "+str(centroids2_sort[i][0])+" "+\
          str(centroids2_sort[i][1])+" "+\
          str(centroids2_sort[i][2])+"\n")
o_file.close()

print "Done"

print "Kabsch minim..."

#translate1 = kabsch_minima.centroid(centroids1_sort)
#translate2 = kabsch_minima.centroid(centroids2_sort)
rmatrix, translate2, translate1 = \
    kabsch_minima.return_rotation_matrix(centroids2_sort, centroids1_sort)

ormatrix = pybel.ob.matrix3x3()
for i in range(3):
  for j in range(3):
    ormatrix.Set(i, j, rmatrix[i,j]) 
myrm = pybel.ob.doubleArray(9)
ormatrix.GetArray(myrm)

mol1 = pybel.readfile("mol2", filenamemol21).next()
mol1.OBMol.Translate(pybel.ob.vector3(-translate1[0], -translate1[1], -translate1[2]));
output = pybel.Outputfile("mol2", fileout1, overwrite=True)
output.write(mol1)
output.close()

mol2 = pybel.readfile("mol2", filenamemol22).next()
mol2.OBMol.Translate(pybel.ob.vector3(-translate2[0], -translate2[1], -translate2[2]));
mol2.OBMol.Rotate(myrm)
output = pybel.Outputfile("mol2", fileout2, overwrite=True)
output.write(mol2)
output.close()

centroids = numpy.empty_like (centroids1_sort)

for i in range(0,len(centroids1_sort)):
  for j in range(0,3):
    centroids[i][j] = centroids1_sort[i][j] - translate1[j]

centroids1_sort[:] = centroids

centroids = numpy.empty_like (centroids2_sort)

for i in range(0,len(centroids2_sort)):
  
  for j in range(0,3):
    sum = 0.0
    for k in range(0,3):
      c = centroids2_sort[i][k] - translate2[k]
      sum = sum + rmatrix[j, k] * c

    centroids[i][j] = sum

centroids2_sort[:] = centroids

print "Check dists "
dists = scipy.spatial.distance.cdist(centroids1_sort, centroids2_sort)
print dists

print "RMSD: ", kabsch_minima.rmsd(centroids1_sort, centroids2_sort)

o_file = open("centroids1_rt.xyz","w")
o_file.write(str(len(centroids1_sort))+"\n")
o_file.write("centroids 1 \n")
for i in range(0,len(centroids1_sort)):
  o_file.write("O   "+str(centroids1_sort[i][0])+" "+\
          str(centroids1_sort[i][1])+" "+\
          str(centroids1_sort[i][2])+"\n")
o_file.close()

o_file = open("centroids2_rt.xyz","w")
o_file.write(str(len(centroids2_sort))+"\n")
o_file.write("centroids 2 \n")
for i in range(0,len(centroids2_sort)):
  o_file.write("O   "+str(centroids2_sort[i][0])+" "+\
          str(centroids2_sort[i][1])+" "+\
          str(centroids2_sort[i][2])+"\n")
o_file.close()
