import os
import sys
import uuid
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

if __name__ == "__main__":

    listaname = unique_filename = str(uuid.uuid4())
    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", \
        help="input the PDB file", \
            required=True, default="", type=str)

    args = parser.parse_args()

    chainlist = ["A", "B"]
    split_PDBfile_by_chains (args.file, chainlist)

    probe = "DRY"
    STEPVAL = 2.0
    DELTAVAL = 10.0
    MINDIM = 30
    NUMOFCLUST = 4

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

    print("Controid")
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