import sys

sys.path.append("./common")
import carbo

STEPVAL = 1.4
DELTAVAL = 20.0
coulombconst = 1.0

filename1 = ""
filename2 = ""
weightsname1 = ""
weightsname2 = ""

if (len(sys.argv)) == 5:
  filename1 = sys.argv[1]
  weightsname1 = sys.argv[2]
  filename2 = sys.argv[3]
  weightsname2 = sys.argv[4]
else:
  print "usage :", sys.argv[0] , " filename1.mol2 weights1.txt filename2.mol2 weights2.txt"
  exit(1)

carbo.carbo_similarity (filename1, weightsname1, filename2, weightsname2, \
        STEPVAL, DELTAVAL, coulombconst)

