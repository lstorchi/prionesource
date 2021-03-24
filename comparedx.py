import numpy as np
import argparse

from gridData import Grid


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f1","--file1", help="input the first dx file", \
            required=True, default="", type=str)
    parser.add_argument("-f2","--file2", help="input the second dx file", \
            required=True, default="", type=str)

    args = parser.parse_args()

    g1 = Grid(args.file1)
    g2in = Grid(args.file2)
    g2 = g2in.resample(g1)

    norig = g1.origin

    g1min = np.min(g1.grid) 
    g2min = np.min(g2.grid) 

    g1max = np.max(g1.grid) 
    g2max = np.max(g2.grid)

    print(g1min, g1max, g2min, g2max) 

    grid1_scaled = np.array([(x - g1min) / (g1max - g1min) for x in g1.grid])
    grid2_scaled = np.array([(x - g2min) / (g2max - g2min) for x in g2.grid])

    g1min = np.min(grid1_scaled) 
    g2min = np.min(grid1_scaled) 

    g1max = np.max(grid1_scaled) 
    g2max = np.max(grid2_scaled)

    print(g1min, g1max, g2min, g2max) 

    ng1 = Grid(grid1_scaled, origin=norig , \
           delta=g1.delta)

    ng2 = Grid(grid2_scaled, origin=norig , \
           delta=g2.delta)

    ng1.export("first.dx")
    ng2.export("second.dx")

    diffg = ng1 - ng2

    print(np.min(diffg.grid), np.max(diffg.grid))

    diffg.export("diff.dx")