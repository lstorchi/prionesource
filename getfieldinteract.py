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

  peratom_counter_energy = []
  peratom_counter_multiple_energy = []

  for ai in range(len(coords)):
    peratom_counter.append(0)
    peratom_counter_multiple.append(0)
    peratom_counter_energy.append(0.0)
    peratom_counter_multiple_energy.append(0.0)

  for iy in range(energy_coords.shape[1]):
    print("%5d of %5d"%(iy, energy_coords.shape[1]), flush=True)
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
                    peratom_counter_multiple_energy[ai] += e

            counter_multiple += partialconter

            if partialconter > 1:
                partialconter = 1

                mindist = min(distfromatom)
                mindistai = distfromatom.index(mindist)

                if mindistai >= 0:
                    peratom_counter[mindistai] += 1
                    peratom_counter_energy[mindistai] += e
                else:
                    print("Error")
                    exit(1)
            elif partialconter == 1:
                notzero = [i for i, v in enumerate(isnearatom) if v != 0]
                if len(notzero) == 1:
                    peratom_counter[notzero[0]] += 1
                    peratom_counter_energy[notzero[0]] += e
                else:
                    print("Error internal ", notzero)
                    exit(1)
 
            counter += partialconter

  print ("     ATOMNAME , RESIDUENAME , ATOMINDEX , " + \ 
    "RESIDUEID , POINTINATOM , TOTALENERGY , MultiPOINTINATOM , " + \
        "MultiTOTALENERGY")
  for ai in range(len(coords)):
      if peratom_counter[ai] != 0 or \
          peratom_counter_multiple[ai] != 0:

              print ("ATOM %5s , %5s , %8d , %5d , %5d , %10.5f , %5d , %10.5f"%(
                  coords[ai].atomname, \
                  coords[ai].resname , \
                  coords[ai].id , \
                  coords[ai].residueid , \
                  peratom_counter[ai], \
                  peratom_counter_energy[ai], \
                  peratom_counter_multiple[ai],
                  peratom_counter_multiple_energy[ai]),flush=True )

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
                eradii[a.element] = mendeleev.element(a.element).vdw_radius/100.0
        for a in coords2:
            if a.element not in eradii:
                eradii[a.element] = mendeleev.element(a.element).vdw_radius/100.0
        print("Done")

        print("  Field1 vs Coord2")
        compare (energycoord1, coords2, ELIMIT, eradii)
        print("  Field2 vs Coord1")
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
    parser.add_argument("-g","--grid", help="Specify a grid \"xmin,xmax,ymin,ymax,zmin,zmax\"", \
        required=False, default="", type=str)
    args = parser.parse_args()

    probe = args.probe
    STEPVAL = args.stepval
    DELTAVAL = args.stepval

    # compute box 
    xmin = 0.0
    xmax = 0.0
    ymin = 0.0
    ymax = 0.0
    zmin = 0.0
    zmax = 0.0

    if args.grid == "":
        xmin, xmax, ymin, ymax, zmin, zmax = \
            gridfield.compute_grid_box (args.file, DELTAVAL)
    else:
        boxs = args.grid.split(",")
        if len(boxs) == 6:
            xmin = float(boxs[0])
            xmax = float(boxs[1])
            ymin = float(boxs[2])
            ymax = float(boxs[3])
            zmin = float(boxs[4])
            zmax = float(boxs[5])
        else:
            print("Error in specified grid")
            exit(1)

    print("Grid to use: ")
    print("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"%(\
        xmin, xmax, ymin, ymax, zmin, zmax))
    print("  Center: %10.5f %10.5f %10.5f"%( \
        ((xmax-xmin)/2)+xmin, \
        ((ymax-ymin)/2)+ymin, \
        ((zmax-zmin)/2)+zmin))
    
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
                args.probe, STEPVAL, verbose, args.savekont)

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