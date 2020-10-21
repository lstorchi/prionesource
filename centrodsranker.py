import os
import sys
import uuid
import numpy
import argparse
import subprocess

import gridfieldcentroids

sys.path.append("./common")
import gridfield

def split_PDBfile_by_chains(filename, chainlist) :

    for c in chainlist:
        toexe = "pdb_selchain -" + c + " " + filename, 
        results  = subprocess.run(toexe, shell=True, check=True, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
            universal_newlines=True)

        f = open(c+".pdb", "w")
        f.write(str(results.stdout))
        f.close()

        toexe = "obabel -ipdb " + c + ".pdb -omol2 -O " + "./"+ c +".mol2"
        results  = subprocess.run(toexe, shell=True, check=True, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
            universal_newlines=True)
        os.remove(c+".pdb")
    
if __name__ == "__main__":

    listaname = unique_filename = str(uuid.uuid4())
    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input the PDB file", \
          required=True, default="", type=str)

    args = parser.parse_args()

    chainlist = ["A", "B"]
    split_PDBfile_by_chains (args.file, chainlist)

    probe = "DRY"
    STEPVAL = 2.0
    DELTAVAL = 10.0
    MINDIM = 30
    NUMOFCLUST = 4

    for c in chainlist:
        f = open(listaname, "w")
        f.write(c + ".mol2 1.0\n")
        f.close()

        energy, xmin, ymin, zmin = gridfield.compute_grid_mean_field (listaname, \
            STEPVAL, DELTAVAL, probe)

        gridfieldcentroids.get_centroids(energy, STEPVAL, NUMOFCLUST, MINDIM, 0.0, 0.0, 0.0)

        os.remove(c+".mol2")
        os.remove(listaname)


    