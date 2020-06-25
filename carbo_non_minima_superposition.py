import re
import os
import sys
import sets
import math
import glob
import numpy
import pybel
import scipy
import subprocess

sys.path.append("./common")
import carbo

if (len(sys.argv)) == 5: 
  filenamemol21 = sys.argv[1]
  weifilename1 = sys.argv[2] 
  filenamemol22 = sys.argv[3]
  weifilename2 = sys.argv[4]
else:
  print "usage :", sys.argv[0] , " filename1.mol2 weight1.txt filename2.mol2 weight2.txt"
  exit(1)

STEPVAL = 1.4
DELTAVAL = 20.0
coulombconst = 1.0

carbo.carbo_similarity (filenamemol21, weifilename1, filenamemol22, weifilename2, \
        STEPVAL, DELTAVAL, coulombconst)
