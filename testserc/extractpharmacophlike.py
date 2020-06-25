import re
import os
import sys
import sets
import math
import glob
import numpy
import pybel
import subprocess

from scipy import cluster
import scipy.spatial.distance

#####################################################################

def is_a_number(s):

    try:
        float(s)
    except ValueError:
        return False
                                    
    return True

#####################################################################

numoff = 0

if (len(sys.argv)) > 2: 
  numoff = int(sys.argv[1])
else:
  print "usage :", sys.argv[0] , " NUM centroid1 cntroid2 ... "
  exit(1)

list_of_centroid = []
coords = []

for struct in range(1, numoff+1):
    #print "reading \"", sys.argv[1+struct], "\""

    fp = open(sys.argv[1+struct], "r") 

    start = False

    list = []
    c = []
    
    for l in fp:
      p = re.compile(r'\s+')
      line = p.sub(' ', l)
      line = line.lstrip()
      line = line.rstrip()

      if start:
          if line == "Done":
              start = False
          else:
              sline = line.split(" ")

              if (is_a_number(sline[0])) and (len(sline) == 7):
                  list.append([float(s) for s in sline])
                  c.append([float(sline[j]) for j in range(3)]) 

      if line == "Controid":
          start = True

    list_of_centroid.append(list)
    coords.append(c)

    fp.close()

alldists = {}

#for struct1 in range(0, numoff):
#    for c in coords[struct1]:
#        print c
#    print ""

for struct1 in range(0, numoff):
    for struct2 in range(struct1, numoff):
        c1 = numpy.asarray(coords[struct1])
        c2 = numpy.asarray(coords[struct2])
        
        dists = numpy.empty((c1.shape[0], c2.shape[0]), dtype=float)

        for i in range(c1.shape[0]):
            for j in range(c2.shape[0]):
                x1 = c1[i][0]
                y1 = c1[i][1]
                z1 = c1[i][2]

                x2 = c2[j][0]
                y2 = c2[j][1]
                z2 = c2[j][2]

                d = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

                dists[i][j] = d

        #dists = scipy.spatial.distance.cdist(c1, c2)

        #for i in range(dists.shape[0]):
        #    for j in range(i+1,dists.shape[1]):
        #        dists[j][i] = dists[i][j]

        #print struct1, struct2 
        #print len(coords[i]), len(coords[j])
        #print dists.shape
        #print dists
        #print ""

        alldists[str(struct1)+"_"+str(struct2)] = dists
        if (struct1 != struct2):
            alldists[str(struct2)+"_"+str(struct1)] = numpy.transpose(dists)

R = 4
MULT = 1.2
MINDISTVAL = 4.0

pointtoprint = []
for struct1 in range(numoff):
    pointtoprint.append(set())

final_list_of_p = []
for struct1 in range(numoff):

    for c1 in range(len(list_of_centroid[struct1])):
        howmanynear = 0
        r1 = list_of_centroid[struct1][c1][R]
        pointsfound = []

        for struct2 in range(numoff):
            atleastonenear = False

            if struct1 != struct2:
                #print struct1, " vs ", struct2
                dists = alldists[str(struct1)+"_"+str(struct2)]
                for c2 in range(len(list_of_centroid[struct2])):
                
                    r2 = list_of_centroid[struct2][c2][R]

                    mindist = ((r1+r2) * MULT)
                    #print mindist
                    mindist = MINDISTVAL

                    if dists[c1][c2] <= mindist:
                        atleastonenear = True
                        pointsfound.append(list_of_centroid[struct2][c2])
                        break

                        #print "struct ( ", struct1, \
                        #        " , ", struct2, ") point (", c1, " , ", c2 , ") "\
                        #        " dist ", dists[c1][c2], " ", r1+r2 
                        #print dists[c1][c2] , (r1+r2), dists[c2][c1]
                
            if atleastonenear:
                howmanynear = howmanynear + 1
                
        if howmanynear == (numoff - 1):
            if not (list_of_centroid[struct1][c1] in final_list_of_p):
                final_list_of_p.append(list_of_centroid[struct1][c1])
            for c in pointsfound:
                if not (c in final_list_of_p):
                    final_list_of_p.append(c)

print len(final_list_of_p)
print ""
for c in final_list_of_p:
    print " H %10.5f %10.5f %10.5f"%(c[0], c[1], c[2])

#for c in list_of_centroid:
#    for l in c:
#        print l   
#    print ""

