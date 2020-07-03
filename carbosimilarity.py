import sys
import numpy
import argparse

import matplotlib.pyplot as plt

sys.path.append("./common")
import carbo

if __name__ == "__main__":
    STEPVAL = 1.4
    DELTAVAL = 20.0

    """
    coulombconst = 1.0
    carbo.read_delphi_file ("wt_90-231_aligned_2.phi")
    exit(1)
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f1","--file1", help="input the first mol2 file and weight filename  \n"+ \
      "   filenames are comma separated \"file1.mol,weight1.txt", required=True, default="", type=str)
    parser.add_argument("-f2","--file2", help="input the second mol2 file and weight filename  \n"+ \
      "   filenames are comma separated \"file2.mol,weight2.txt", required=True, default="", type=str)
    parser.add_argument("-s","--stepvalue", help="input step value defaul="+str(STEPVAL), \
      required=False, default=STEPVAL, type=float)
    parser.add_argument("-d","--deltaval", help="input delta value defaul="+str(DELTAVAL), \
      required=False, default=DELTAVAL, type=float)
    parser.add_argument("-v","--verbose", help="Verbose mode on", \
      required=False, default=False, action="store_true")
    parser.add_argument("--show", help="Show plot", \
      required=False, default=False, action="store_true")
    parser.add_argument("--damped", help="Use damped approache", \
      required=False, default=False, action="store_true")

    args = parser.parse_args()

    if len(args.file1.split(",")) != 2:
      print("Error in --file1 option two comma separated filenames expected")
      exit(1)

    if len(args.file2.split(",")) != 2:
      print("Error in --file2 option two comma separated filenames expected")
      exit(1)

    filename1 = args.file1.split(",")[0]
    weightsname1 = args.file1.split(",")[1]

    filename2 = args.file2.split(",")[0]
    weightsname2 = args.file2.split(",")[1] 
  
    try:
      carboidxs, xrefpoints, weights, pweights = carbo.carbo_similarity (\
        filename1, weightsname1, filename2, weightsname2, \
          args.stepvalue, args.deltaval, coulombconst, args.verbose, \
            args.damped)

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
      
