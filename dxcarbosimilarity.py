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
      carboidxs, xrefpoints, weights, pweights = \
        carbo.returncarbodxs(set1, set2, args.verbose)

      stdev = carboidxs.std(0)
      meanmtx = carboidxs.mean(0)
  
      print("")
      print("Final Mean and stdev")
      for idx, std in enumerate(stdev):
        print("%+11.8e %+11.8e %+11.8e"%(xrefpoints[idx], meanmtx[idx] , std))

      plt.errorbar(xrefpoints, meanmtx, stdev,  linestyle='None', \
        marker='^', label="Mean and stdev")
      plt.plot(xrefpoints, meanmtx, linestyle='--')
  
      #print("Full matrix")
      #print(carboidxs.mean(0))
      #print(carboidxs.std(0))
  
      print("Weighted Mead and stdev")
      waverage = numpy.average(carboidxs, 0, weights)
      wvariance = numpy.average((carboidxs-waverage)**2, 0, weights)
      for idx, std in enumerate(wvariance):
        print("%+11.8e %+11.8e %+11.8e"%(xrefpoints[idx], waverage[idx] , std))

      plt.errorbar(xrefpoints, waverage, wvariance,  linestyle='None', \
        marker='^', label="Weighted Mean and stdev")
      plt.plot(xrefpoints, waverage, linestyle='--')
  
      #print("full matrix")
      #print(waverage)
      #print(wvariance)
  
      print("PWeighted Mead and stdev")
      waverage = numpy.average(carboidxs, 0, pweights)
      wvariance = numpy.average((carboidxs-waverage)**2, 0, pweights)
      for idx, std in enumerate(wvariance):
        print("%+11.8e %+11.8e %+11.8e"%(xrefpoints[idx], waverage[idx] , std))
      
      plt.errorbar(xrefpoints, waverage, wvariance,  linestyle='None', \
        marker='^', label="PWeighted Mean and stdev")
      plt.plot(xrefpoints, waverage, linestyle='--')
      
      #print("full matrix")
      #print(waverage)
      #Sprint(wvariance)
      
      if args.show:
        plt.legend(loc="lower left")
        plt.show()

    except Exception as exp:
      print(exp)