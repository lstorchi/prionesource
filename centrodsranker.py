import os
import sys
import uuid
import math
import numpy
import argparse
import subprocess

import gridfieldcentroids

sys.path.append("./common")
import gridfield

###############################################################

def split_PDBfile_by_chains(filename, chainlist) :

    for c in chainlist:
        toexe = "pdb_selchain -" + c + " " + filename, 
        results  = subprocess.run(toexe, shell=True, check=True, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
            universal_newlines=True)

        f = open(c+".pdb", "w")
        f.write(str(results.stdout))
        f.close()

###############################################################

def dist (c1, c2):

    d2 = (c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2  

    return math.sqrt(d2)

###############################################################


if __name__ == "__main__":

    probe = "DRY"
    STEPVAL = 2.0
    DELTAVAL = 10.0
    MINDIM = 30
    NUMOFCLUST = 4

    listaname = unique_filename = str(uuid.uuid4())
    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", \
        help="input the PDB file", \
            required=True, default="", type=str)
    parser.add_argument("-p","--probe", help="the probe to be used [default="+probe+"]", \
      required=False, default=probe, type=str)
    parser.add_argument("-s","--stepval", help="input stepval value [defaul="+str(STEPVAL)+"]", \
      required=False, default=STEPVAL, type=float)
    parser.add_argument("-d","--deltaval", help="input deltaval value [defaul="+str(DELTAVAL)+"]", \
      required=False, default=DELTAVAL, type=float)
    parser.add_argument("-m","--numofdim", help="input number of minima to be used [defaul="+ \
      str(MINDIM)+"]", required=False, default=MINDIM, type=int)
    parser.add_argument("-c","--numofcluster", help="input number of clusters to be used [defaul=" + \
      str(NUMOFCLUST)+"]", required=False, default=NUMOFCLUST, type=int)

    args = parser.parse_args()

    probe = args.probe
    MINDIM = args.numofdim
    NUMOFCLUST = args.numofcluster
    STEPVAL = args.stepval
    DELTAVAL = args.stepval

    chainlist = ["A", "B"]
    split_PDBfile_by_chains (args.file, chainlist)


    valuefp = []

    for c in chainlist:
        f = open(listaname, "w")
        f.write(c + ".pdb 1.0\n")
        f.close()

        energy, xmin, ymin, zmin = gridfield.compute_grid_mean_field (listaname, \
            STEPVAL, DELTAVAL, probe, False)

        os.remove(listaname)

        centroids, centroidvals, rmins, rmaxs, ravgs = gridfieldcentroids.get_centroids(energy, \
            STEPVAL, NUMOFCLUST, MINDIM, 0.0, 0.0, 0.0)

        valuefp.append((centroids, centroidvals, rmins, rmaxs, ravgs))
        os.remove(c+".pdb")

    print("Controids")
    print("X          Y          Z          EMIN       RMIN       RMAX       RAVG")
    for v in valuefp:
        centroids = v[0]
        centroidvals = v[1]
        rmin = v[2]
        rmax = v[3]
        ravg = v[4]

        for j in range(len(centroids)):
          print("%10.5f"%centroids[j][0], \
              "%10.5f"%centroids[j][1], \
              "%10.5f"%centroids[j][2], \
              "%10.5f"%centroidvals[j], \
              "%10.5f"%rmins[j], \
              "%10.5f"%rmaxs[j], \
              "%10.5f"%ravgs[j])
        print("")

    for i in range(len(valuefp)):
        cxyz1 = valuefp[i][0]
        cval1 = valuefp[i][1]
        rmin1 = valuefp[i][2]
        rmax1 = valuefp[i][3]
        ravg1 = valuefp[i][4]
        for j in range(i+1, len(valuefp)):
            cxyz2 = valuefp[j][0]
            cval2 = valuefp[j][1]
            rmin2 = valuefp[j][2]
            rmax2 = valuefp[j][3]
            ravg2 = valuefp[j][4]

            summa = 0.0
            isin = 0
            for idx1, c1 in enumerate(cxyz1):
                for idx2, c2 in enumerate(cxyz2):
                    s = dist(c1, c2) * cval1[idx1] * cval2[idx2]
                    summa += s
                    if dist(c1, c2) < (rmax1[idx1] + rmax2[idx2]):
                        isin +=  1
                    print("%10.5f"%(s))

            print ("Pair %3d %3d has Summa %10.5f  InDist %4d "%(i, j, summa, isin))