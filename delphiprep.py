import os
import sys
import argparse
import subprocess

sys.path.append("./common")
import carbo

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input the mol2 file", \
            required=True, default="", type=str)

    args = parser.parse_args()

    basename = os.path.splitext(args.file)[0]

    atoms = carbo.mol2atomextractor(args.file)
    
    result = subprocess.run("obabel -imol2 "+  basename+ ".mol2 " + \
            "-opdb -O " + basename+".pdb", shell=True, check=True, \
            stdout=subprocess.PIPE, universal_newlines=True)

    print (atoms)



      
