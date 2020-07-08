import sys
import numpy
import argparse

import matplotlib.pyplot as plt
from gridData import Grid

sys.path.append("./common")
import carbo

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f1","--file1", help="input the first files list and weight values", \
      required=True, default="", type=str)
    parser.add_argument("-f2","--file2", help="input the second files list and weight values", \
      required=True, default="", type=str)
    parser.add_argument("-v","--verbose", help="Verbose mode on", \
      required=False, default=False, action="store_true")
    parser.add_argument("--show", help="Show plot", \
      required=False, default=False, action="store_true")
    
    args = parser.parse_args()

    set1 = {}
    with open(args.file1) as fp:
      for line in fp:
        dxfname = line.split()[0]
        w = float(line.split()[1])
        g = Grid(dxfname)

        set1[dxfname] = (w, g)

    set2 = {}
    with open(args.file2) as fp:
      for line in fp:
        dxfname = line.split()[0]
        w = float(line.split()[1])
        g = Grid(dxfname)

        set2[dxfname] = (w, g)
  
    try:
      carbo.returncarbodxs(set1, set2, args.verbose)
    except Exception as exp:
      print(exp)
      
