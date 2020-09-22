import sys
import math
import numpy
import argparse

import matplotlib.pyplot as plt
from gridData import Grid

sys.path.append("./common")
import carbo

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input the dx file", \
      required=True, default="", type=str)
    
    args = parser.parse_args()

    g = Grid(args.file)

    maxv = numpy.max(g.grid)
    minv = numpy.min(g.grid)
    print("Max %10.4f Min %10.4f"%(maxv, minv))

    for i in range(g.grid.shape[0]):
      for j in range(g.grid.shape[1]):
        for k in range(g.grid.shape[2]):
          if g.grid[i, j, k] > 0.0:
            g.grid[i, j, k] /= maxv
          elif g.grid[i, j, k] < 0.0:
            g.grid[i, j, k] /= math.fabs(minv)
    
    g.export(args.file)