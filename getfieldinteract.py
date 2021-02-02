import os
import sys
import uuid
import math
import numpy
import argparse
import mendeleev
import subprocess

import gridfieldcentroids

sys.path.append("./common")
import gridfield
import carbo



###############################################################

def split_PDBfile_by_chains(filename, chainlist) :

    chainsfile = []

    for c in chainlist:
        toexe = "pdb_selchain -" + c + " " + filename, 
        results  = subprocess.run(toexe, shell=True, check=True, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
            universal_newlines=True)

        basicname = filename[:-4]
        cname = basicname + "_"+c+".pdb"
        f = open(cname, "w")
        f.write(str(results.stdout))
        f.close()
        
        chainsfile.append(cname)

    return chainsfile

###############################################################

def dist (c1, c2):

    d2 = (c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2  

    return math.sqrt(d2)


###############################################################

def compare (energy_coords, coords, ELIMIT, eradii):

  VDW = 1.0

  counter = 0
  counter_multiple = 0
  peratom_counter = []
  peratom_counter_multiple = []

  for iy in range(energy_coords.shape[1]):
    for ix in range(energy_coords.shape[0]):
        for iz in range(energy_coords.shape[2]):
            x, y, z, n = energy_coords[ix, iy, iz] 
            e = energy[ix, iy, iz]
            partialconter = 0
            distfromatom = []
            isnearatom = []

            for ai in range(len(coords)):
                ax = coords[ai].coords[0]
                ay = coords[ai].coords[1]
                az = coords[ai].coords[2]

                dist = (ax-x)**2 + (ay-y)**2 + (az-z)**2 
                distfromatom.append(dist)
                isnearatom.append(0)

                #radii = mendeleev.element(coords[ai].element).vdw_radius
                radii = eradii[coords[ai].element]

                if dist < VDW*radii and e <= ELIMIT:
                    partialconter += 1
                    isnearatom[ai] = 1
                    peratom_counter_multiple[ai] += 1

            counter_multiple += partialconter

            if partialconter > 1:
                #print partialconter
                partialconter = 1

                mindistai = -1
                mindist = 0.0
                for ai in range(len(coords)):
                    if isnearatom[ai] != 0:
                        if mindistai < 0:
                            mindist = distfromatom[ai]
                            mindistai = ai
                        else:
                            if distfromatom[ai] < mindist:
                                mindist = distfromatom[ai]
                                mindistai = ai

                if mindistai >= 0:
                    peratom_counter[mindistai] += 1
                else:
                    print("Error")
                    exit(1)

            elif partialconter  == 1:
                for ai in range(len(coords)):
                    if isnearatom[ai] != 0:
                        peratom_counter[ai] += 1

            counter += partialconter

###############################################################

def get_comp_values(first, energy1, energycoord1, \
    second, energy2, energycoord2, ELIMIT):

    mols1 = carbo.pdbatomextractor (first)
    mols2 = carbo.pdbatomextractor (second)

    if len(mols1) == 1 and len(mols2) == 1:
        coords1 = mols1[0]
        coords2 = mols2[0]

        eradii = {}

        print("Getting elements radii")
        for a in coords1:
            if a.element not in eradii:
                eradii[a.element] = mendeleev.element(a.element).vdw_radius
        for a in coords2:
            if a.element not in eradii:
                eradii[a.element] = mendeleev.element(a.element).vdw_radius
        print("Done")

        compare (energycoord1, coords2, ELIMIT, eradii)
        compare (energycoord2, coords1, ELIMIT, eradii)
    else:
        print("Moltiple models in PDB ?")

###############################################################


if __name__ == "__main__":

    probe = "DRY"
    STEPVAL = 2.0
    DELTAVAL = 10.0
    MINDIM = 30
    NUMOFCLUST = 4
    verbose = True

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", \
        help="input filelist each file is splitted into A and B chain", \
            required=True, default="", type=str)
    parser.add_argument("-p","--probe", help="the probe to be used [default="+probe+"]", \
      required=False, default=probe, type=str)
    parser.add_argument("-s","--stepval", help="input stepval value [defaul="+str(STEPVAL)+"]", \
      required=False, default=STEPVAL, type=float)
    parser.add_argument("-e","--eminilimit", help="input the emin minimun value", \
        required=True, default=0.0, type=float)
    parser.add_argument("-k", "--savekont", help="save kont file during the computation", \
        default=False, action="store_true")
    parser.add_argument("-d","--deltaval", help="input deltaval value [defaul="+str(DELTAVAL)+"]", \
        required=False, default=DELTAVAL, type=float)
    args = parser.parse_args()

    probe = args.probe
    STEPVAL = args.stepval
    DELTAVAL = args.stepval

    # compute box 
    xmin, xmax, ymin, ymax, zmin, zmax = \
        gridfield.compute_grid_box (args.file, DELTAVAL)

    print("Grid to use: ")
    print("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"%(\
        xmin, xmax, ymin, ymax, zmin, zmax))
    
    fp = open(args.file, "r")

    namestofields = {}
    chainslist = []

    for filename in fp:
        chainlist = ["A", "B"]
        filename = filename.replace("\n", "")
        chainsfile = split_PDBfile_by_chains (filename, chainlist)

        chainslist.append(chainsfile)

        for cname in chainsfile:
            energy, energy_coords = gridfield.compute_grid_field (cname, \
                (xmin, xmax, ymin, ymax, zmin, zmax),
                args.probe, STEPVAL, verbose)

            namestofields[cname] = (energy, energy_coords)
            

    for pair in chainslist:
        first = pair[0]
        second = pair[1]

        print(first, " VS ", second)

        get_comp_values(first, namestofields[first][0], namestofields[first][1], \
            second, namestofields[second][0], namestofields[second][1], args.eminilimit)

        if os.path.isfile(first):
            os.remove(first)

        if os.path.isfile(second):
            os.remove(second)


        
            
        