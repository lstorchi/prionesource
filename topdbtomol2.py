import os
import sys
import glob
import numpy
import argparse
import subprocess

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--pdbfile", help="input the pdb file", 
            required=True, default="", type=str)
    parser.add_argument("-t","--topfile", help="input the top file", \
            required=True, default="", type=str)

    args = parser.parse_args()

    pdbfile = args.pdbfile
    topfile = args.topfile

    basename = pdbfile[0:-4]

    
    print("Producing mol2 files")
    result = subprocess.run("obabel -ipdb "+ pdbfile + \
            " -omol2 -O " + basename+".mol2", shell=True, check=True, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,  \
                    universal_newlines=True)

    fp = open(topfile)

    atoms = []
    stateinatom = False
    for l in fp:
        if l.find("[ atoms ]") == 0:
            stateinatom = True
        elif l.find("[ bonds ]") == 0:
            stateinatom = False
        
        if stateinatom:
            sl = l.split(";")
            values = sl[0].split()
            if len(values) == 8:
                atoms.append(values)

    fp.close()

    fp = open(basename+".mol2")

    stateinatom = False
    for l in fp:
        if l.find("@<TRIPOS>ATOM") == 0:
            stateinatom = True
            print(l, end="")
        elif l.find("@<TRIPOS>BOND") == 0:
            stateinatom = False

        if stateinatom:
            values = l.split()
            if len(values) == 9:
                idx = int(values[0])
                if idx == int(atoms[idx-1][0]) and \
                    values[7].find(atoms[idx-1][3][0:2]) == 0 :
                    print("%7d  %8s  %8.4f  %8.4f  %8.4f %6s  %3d  %10s  %8.4f"%( \
                        idx, values[1], float(values[2]), float(values[3]), float(values[4]), \
                            values[5], int(values[6]), values[7], float(atoms[idx-1][6])))
                else:
                    print("Error in id")
                    print(l)
                    print(atoms[idx-1])
                    exit(0)
            else:
                print(l, end="")
        else:
            print(l, end="")


        
