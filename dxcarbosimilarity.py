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
    parser.add_argument("--axis", help="Select axis to be used for the C similarity [default=x]", \
      required=False, default="x", type=str)
    
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
        carbo.returncarbodxs(set1, set2, args.verbose, args.axis)
        
      stdev = carboidxs.std(0)
      meanmtx = carboidxs.mean(0)
  
      plt.errorbar(xrefpoints, meanmtx, stdev,  linestyle='None', \
        marker='^', label="Mean and stdev")
      plt.plot(xrefpoints, meanmtx, linestyle='--')
  
      #print("Full matrix")
      #print(carboidxs.mean(0))
      #print(carboidxs.std(0))
  
      waverage = numpy.average(carboidxs, 0, weights)
      wvariance = numpy.average((carboidxs-waverage)**2, 0, weights)
      plt.errorbar(xrefpoints, waverage, wvariance,  linestyle='None', \
        marker='^', label="Weighted Mean and stdev")
      plt.plot(xrefpoints, waverage, linestyle='--')
  
      #print("full matrix")
      #print(waverage)
      #print(wvariance)
  
      pwaverage = numpy.average(carboidxs, 0, pweights)
      pwvariance = numpy.average((carboidxs-waverage)**2, 0, pweights)
     
      plt.errorbar(xrefpoints, pwaverage, pwvariance,  linestyle='None', \
        marker='^', label="PWeighted Mean and stdev")
      plt.plot(xrefpoints, pwaverage, linestyle='--')

      print("%13s %13s %13s %13s %13s %13s %13s"%("X", "SMean",  \
        "SStdev", "WMean", "WStdev", "PWMean", "PWStdev"))
      for idx, std in enumerate(stdev):
        print("%+8.6e %+8.6e %+8.6e %+8.6e %+8.6e %+8.6e %+8.6e"%(\
          xrefpoints[idx], meanmtx[idx] , std, waverage[idx], \
            wvariance[idx], pwaverage[idx], pwvariance[idx] ))
      
      #print("full matrix")
      #print(waverage)
      #Sprint(wvariance)
      
      if args.show:
        plt.legend(loc="lower left")
        plt.show()
      else:
        plt.legend(loc="lower left")
        #plt.figure(num=None, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
        plt.savefig('finalcarbo.png')

    except Exception as exp:
      print(exp)