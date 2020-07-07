import os
import sys
import argparse
import subprocess

import mendeleev

sys.path.append("./common")
import carbo

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input the mol2 file", \
            required=True, default="", type=str)

    args = parser.parse_args()

    basename = os.path.splitext(args.file)[0]

    mols = carbo.mol2atomextractor(args.file, True)
    
    result = subprocess.run("obabel -imol2 "+  basename+ ".mol2 " + \
            "-opqr -O " + basename+".pqr", shell=True, check=True, \
            stdout=subprocess.PIPE, universal_newlines=True)

    exit()

    fpcrg = open(basename+".crg", "w")
    fpsiz = open(basename+".siz", "w")
    fpcrg.write("atom__resnumbc_charge_ \n")
    fpsiz.write("atom__res_radius_\n")
    for atoms in mols:
        for a in atoms:
                name = "{:<5}".format(a.name)
                element = mendeleev.element(a.atomname)
                fpcrg.write(name + " %5s %+8.3f\n"%(a.resname, a.partialcharge))
                fpsiz.write(name + " %5s %8.3f\n"%(a.resname, element.vdw_radius/100.0))

    fpcrg.close()
    fpsiz.close()
