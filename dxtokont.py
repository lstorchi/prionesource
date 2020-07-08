from gridData import Grid
import numpy
import argparse

def to_kont(g):

  idx = 1
  for k in range(g.grid.shape[2]):
    z = g.origin[2] + k*g.delta[2]
    for i in range(g.grid.shape[0]):
      x = g.origin[0] + i*g.delta[0]
      for j in range(g.grid.shape[1]):
        y = g.origin[1] + j*g.delta[2]

        print("%6d %8.3f %8.3f %8.3f"%(idx, x, y, z))
        idx += 1

  print("Fake")

  for k in range(g.grid.shape[2]):
    for i in range(g.grid.shape[0]):
      for j in range(g.grid.shape[1]):

        print("%8.3f"%(g.grid[i, j, k]))

def to_txt(g):

  for i in range(g.grid.shape[0]):
    x = g.origin[0] + i*g.delta[0]
    for j in range(g.grid.shape[1]):
      y = g.origin[1] + j*g.delta[2]
      for k in range(g.grid.shape[2]):
        z = g.origin[2] + k*g.delta[2]

        print("%8.3f %8.3f %8.3f %8.3f"%(x, y, z, g.grid[i, j, k]))

if __name__ ==  "__main__":   
  
  parser = argparse.ArgumentParser()
  parser.add_argument("-f","--file", help="input the first DX file", \
    required=True, default="", type=str)

  args = parser.parse_args()

  g = Grid(args.file)

  to_kont(g)
