import os
import argparse
import subprocess

def split_PDBfile_by_chains(filename, chainlist) :

    for c in chainlist:
        toexe = "pdb_selchain -" + c + " " + filename, 
        results  = subprocess.run(toexe, shell=True, check=True, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
            universal_newlines=True)

        f = open(c+".pdb", "w")
        f.write(str(results.stdout))
        f.close()
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input the PDB file", \
          required=True, default="", type=str)

    args = parser.parse_args()

    chainlist = ["A", "B"]
    split_PDBfile_by_chains (args.file, chainlist)

    for c in chainlist:

        os.remove(c+".pdb")


    